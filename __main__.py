# -*- coding: utf-8 -*-

import atexit
import os
import sys
import json
import fnmatch
import time
import logging
import logging.handlers
import threading

from cli import lie_cli

# Package info
__module__ = 'eTOXlie'
__docformat__ = 'restructuredtext'
__version__ = '{major:d}.{minor:d}.{micro:d}'.format(major=1, minor=0, micro=0)
__author__ = 'Marc van Dijk'
__status__ = 'release 1'
__date__ = '20 june 2017'
__copyright__ = 'Copyright (c) 2017, VU University, Amsterdam, the Netherlands'
__rootpath__ = os.path.dirname(os.path.abspath(__file__))
__pyv__ = sys.version_info[0:2]

verbLevel={
    'debug':logging.DEBUG,
    'info':logging.INFO,
    'warn':logging.WARNING,
    'error':logging.ERROR,
    'critical':logging.CRITICAL
}

# Check if Python virtual environment is in sys.path
venv_path = os.getenv('VIRTUAL_ENV', '{0}/.etox_venv'.format(__rootpath__))
venv_active = False
for path in sys.path:
    packages = '{0}*/site-packages'.format(venv_path)
    if fnmatch.fnmatch(path, packages):
        venvpath = path
        venv_active = True
if not venv_active:
    print('Python virtual environment not active. Activate using the shell (bash) command:')
    print('source {0}/.etox_venv/bin/activate'.format(__rootpath__))
    sys.exit(1)

# Init base logger
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')

# Add eTOXlie application to Python path
if __rootpath__ not in sys.path:
    sys.path.append(__rootpath__)
    logging.debug('Add eTOXlie to Python Path')

# Parse settings
json_settings = os.path.join(__rootpath__, 'data/settings.json')
if not os.path.exists(json_settings):
    logging.error('eTOXlie settings file not found: {0}'.format(json_settings))
    exit(1)

settings = {}
with open(json_settings) as st:
    settings = json.load(st)

def bootstrap_app(args):
    
    # Update global settings with CLI arguments
    for key,value in args.__dict__.items():
        if key in settings:
            settings[key] = value
    
    # Update the eTOXlie project folder
    etoxlie_folder = settings.get('etoxlie_folder', '')
    if not len(etoxlie_folder):
        etoxlie_folder = os.path.abspath('{0}/eTOXLie_projects'.format(os.environ['HOME']))
    
    logging.debug('eTOXlie project directory: {0}'.format(etoxlie_folder))
    if not os.path.exists(etoxlie_folder):
        os.makedirs(etoxlie_folder)
        logging.info('eTOXlie project directory made at: {0}'.format(etoxlie_folder))
        
    settings['etoxlie_folder'] = etoxlie_folder
    
    # Adjust the logger
    if args.etoxlie_loglevel:
        logging.info('Set application log level to: {0}'.format(args.etoxlie_loglevel))
        logging.basicConfig(level=verbLevel[args.etoxlie_loglevel.lower()], format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
    
    # Check if essential external programs are available
    unavailable_exe = []
    for exe in ('ACPYPE','AMBERHOME','GROMACSHOME','GMXRC','PARADOCKS','PLANTS'):
        path = settings.get(exe, '')
        if not path:
            logging.warn('Path to program {0} not defined'.format(exe))
            unavailable_exe.append(exe)
            continue
        path = os.path.abspath(path)
        if not os.path.exists(path):
            logging.warn('Program {0} not available at path {1}'.format(exe,path))
            unavailable_exe.append(exe)
        else:
            logging.info('Program {0} available at {1}'.format(exe,path))
    
    if 'PLANTS' in unavailable_exe and 'PARADOCKS' in unavailable_exe:
        logging.error('Docking programs PLANTS and PARADOCKS not defined, unable to continue')
        exit(1)
    essential = set(unavailable_exe).intersection(('ACPYPE','AMBERHOME','GROMACSHOME','GMXRC'))
    if essential:
        logging.error('Essential programs not defined: {0}'.format(','.join(list(essential))))
    
    # Save modified settings to JSON file
    if args.save_settings:
        with open(json_settings, 'w') as st:
            json.dump(settings, st, indent=2, sort_keys=True)
            logging.info('Update eTOXlie settings.json file: {0}'.format(json_settings))
    
def eTOXlie_deamon(settings):
    
    logging.info('Start a new eTOXlie deamon')
    
    # Setup for remote execution
    remoteSettings = None
    sshConn = None
    localExe = settings.get('etoxlie_localExe', True)
    
    if not localExe:
        try:
            remoteSettings = {
                'submit' : settings['SUBMITJOB'],
                'check' : settings['CHECKJOB'],
                'cancel' : settings['CANCELJOB'],
                'jobHead' : settings['JOB_HEAD'],
                'wdir' : settings['JOB_WORKDIR'],
                'tdir' : settings['JOB_TEMPDIR']
            }
            
            sshConn=sshHandler.Connection(
                settings['CLUSTER'],
                settings['ACCOUNT'],
                password=settings.get('PASSWORD') or None
            )
            if not sshConn.connected:
                raise Exception(sshConn.info)
        except Exception, e:
            localExe=True
            sshConn=None
            logging.error('Error in configuring the remote execution: %s'%e)
    
    killer = jobHandler.WarmKiller()
    while True:
    
        # Check jobs. The function sort the jobs according to the submission date
        listJobs=jobHandler.collectJobs(settings['etoxlie_folder'], sort=True)
        
        # Cycle through jobs
        for job in listJobs:
            logging.debug("Job submission time: {0}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(job.dateSub))))
            logging.debug("Job status: %s"%job.status)
            logging.debug("Job filename: %s"%job.filename)
            jobUpdated=calculateInteractions(job,
                localExe=localExe,
                sshConn=sshConn,
                killer=killer,
                dock_soft=settings.get('etoxlie_docking_method', 'PLANTS'),
                timemd=.1,
                clean=False,
                remoteSettings=remoteSettings,
                modelFN='model.dat'
            )
            
            # kill if getting signal after each job
            if killer.kill_now:
                logging.info("%s Warm Killing"%time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time())))
                if not localExe:
                    sshConn.close()
                sys.exit(0)

        # In case no jobs are present
        if killer.kill_now:
            logging.info("%s Warm Killing"%time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time())))
            if not localExe:
                sshConn.close()
            sys.exit(0) 
        
        time.sleep(settings.get('etoxlie_deamon_updateinterval', 600))
    
    
if __name__ == '__main__':
    
    # Launch eTOXlie from command line.
    # Parse CLI arguments and bootstrap
    bootstrap_app(lie_cli(__rootpath__, defaults=settings))
    
    # Always log a start
    verbosity = verbLevel[settings['etoxlie_loglevel']]
    logging.log(verbosity, 'Starting of a new LIE daemon')
    logging.log(verbosity, 'Deamon update frequency in seconds: {0}'.format(settings['etoxlie_deamon_updateinterval']))
    logging.log(verbosity, 'Folder for temporary files: {0}'.format(settings['etoxlie_folder']))
    
    # Setup the file logger
    from etox.core import jobHandler, sshHandler
    from etox.core.main import calculateInteractions
    
    logfile = os.path.join(__rootpath__, 'data/logs', settings.get('etoxlie_logfile', 'eTOXlie.log'))
    lielogger = logging.getLogger()
    lielogger.setLevel(verbosity)
    handler = logging.handlers.RotatingFileHandler(
        logfile,
	    maxBytes=settings.get('etoxlie_logfile_maxsize', 10*1024*1024),
        backupCount=settings.get('etoxlie_logfile_maxbackup', 100))
    handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s: %(message)s'))
    lielogger.addHandler(handler)
    
    # Run the eTOXlie deamon on a new tread
    logging.debug('Start the eTOXlie worker deamon')
    deamon_thread = threading.Thread(target=eTOXlie_deamon, args=[settings])
    deamon_thread.daemon = True
    deamon_thread.start()
    time.sleep(1)
    
    # Launch Flask app
    if settings.get('run_gui', True):    

        from gui.app import app
        app.config.update(**settings)
        app.run(host=settings.get('flask_host'),
                port=settings.get('flask_port'),
                debug=settings.get('etoxlie_loglevel') == 'debug',
                use_reloader=settings.get('flask_reloader', False))