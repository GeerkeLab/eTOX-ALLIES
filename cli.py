# -*- coding: utf-8 -*-

"""
Command Line Interface to the eTOXlie application
"""

import argparse
import os


def lie_cli(root_path, prog="eTOXlie", defaults={}):
    """
    Command Line Interface to the eTOXlie application
    
    :param root_path: path to the application directory
    :type root_path:  str
    :param prog:      CLI program name
    :type prog:       str
    """
    
    # Create the top-level parser
    parser = argparse.ArgumentParser(prog=prog,
                                     description="eTOXlie, binding affinity prediction.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Set global application logger
    parser.add_argument('--loglevel',
                        type=str,
                        dest='etoxlie_loglevel',
                        default=defaults.get('loglevel','info'),
                        choices=['error', 'warn', 'info', 'debug','critical'],
                        help=("Global application log level"))
    parser.add_argument('--logfile',
                        type=str,
                        dest='etoxlie_logfile',
                        default=defaults.get('etoxlie_logfile'),
                        help=("Global application log file"))
    
    # Save modified settings to JSON file again
    parser.add_argument('--save-settings',
                        type=bool,
                        default=True,
                        help=("Save command line settings to eTOXlie settings.json file"))
    
    # parse cmd line args
    options = parser.parse_args()
    
    return options
