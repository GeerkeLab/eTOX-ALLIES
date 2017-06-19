import os
import sys
import paramiko
import select
import logging

class Connection:
    def __init__(self, address, account, password=None):
        
        self.address = address
        self.account = account
        self.password = password
        
        self.sshtunnel()
        logging.debug(self.info)


    def sshtunnel(self):
        maxcycles=100
        try:
            '''open ssh tunnel to the cluster using the authorized_key credentials'''
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            '''In case the server's key is unknown,
            we will be adding it automatically to the list of known hosts'''
            ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
            
            if self.password:
                '''SSH connect based on password'''
                ssh.connect(self.address, username=self.account, password=self.password)
            else:            
                '''Loads the user's local known host file.'''
                ssh.connect(self.address, username=self.account)
            
            sftp = ssh.open_sftp()
            self.home=sftp.getcwd()
            self.ssh=ssh
            self.sftp=sftp
            self.connected=True
            self.info='connected to %s@%s'%(self.account,self.address)      
        
        except Exception,e:
            msg="Error in open ssh connection: %s"%e
            logging.error(msg)
            self.connected=False
            self.ssh=None
            self.sftp=None
            self.info=msg


    def checkAlive(self):
        if self.ssh.get_transport().is_active():
            return
        else:
            self.sshtunnel()
            return


    def execCmd(self,cmd):
        self.checkAlive()
        if self.connected:
            stdOut=''
            transport = self.ssh.get_transport()
            channel=transport.open_session()
            channel.exec_command(cmd)
            while not channel.exit_status_ready():               
                # Only print data if there is data to read in the channel
                if channel.recv_ready():
                    outLine=channel.recv(1024)
                    #rl, wl, xl = select.select([channel], [], [], 0.0)
                    #if len(rl) > 0:
                        # Print data from stdout
                    stdOut+=outLine
                    logging.debug("OUT: %s"%outLine)
                if channel.recv_stderr_ready():
                #if ssh_stderr.channel.recv_ready():    
                    #rl, wl, xl = select.select([channel.recv_stderr(-1)], [], [], 0.0)
                    errLine=channel.recv_stderr(1024)
                    #if len(rl) > 0:
                    logging.debug("ERR: %s"%errLine)

            outLine=channel.recv(1024)
            stdOut+=outLine            
            exitStatus=channel.recv_exit_status()
            
            return (True,(exitStatus,stdOut))

        else:
            msg='No connection available'
            return (False,msg)

                
    def putFile(self,oriFile,destDir, overwrite=True):
        self.checkAlive()
        if self.connected:
            ## check absolute path or not:
            if not os.path.abspath(destDir)==destDir:
                self.sftp.chdir('.')
                RcurrDir=self.sftp.getcwd()
                destDir=os.path.join(RcurrDir,destDir)
            
            ## Check if directory exists
            try:
                self.sftp.chdir(destDir)
            except IOError:
                logging.debug('Destination directory does not exists. Trying to create it')
                fullpath=[]
                head=destDir
                while True:
                    head,tail=os.path.split(head)
                    if tail =='':
                        break
                    fullpath.append(tail)
                
                fullpath.reverse()
                head='/'
                for folder in fullpath:
                    head=os.path.join(head,folder)
                    try:
                        self.sftp.chdir(head)
                    except:
                        try:
                            self.sftp.mkdir(head)
                        except Exception,e:
                            msg='Error in creation the directory %s: %s'%(head,e)
                            return (False,msg)
                

            destFile=os.path.join(destDir,os.path.split(oriFile)[1])
            try:    # if file exists:
                self.sftp.stat(destFile)
                if overwrite==True:
                    logging.debug('Overwriting existing file: %s:%s'%(self.address,destFile))
                    self.sftp.remove(destFile)
                else:
                    msg='File %s:%s already exists. overwrite=True if you wanna it.'%(self.address,destFile)
                    return (False,msg)
            except:
                pass    # File does not exist

            try:
                self.sftp.put(oriFile,destFile)
                permOri=oct(os.stat(oriFile).st_mode)[-4:]
                self.sftp.chmod(destFile, int(permOri,8))       # Change permissions to file
            except Exception, e:
                msg="Error during copy of the file %s to %s: %s"%(oriFile,destDir,e)
                return (False,msg)

            msg='Successful transfer of file %s to %s'%(oriFile,destDir)
            return (True,msg)

        else:
            msg='Transfer of file %s to %s:%s failed because no connection was possible'%(oriFile,self.address,destDir)
            return (False,msg)


    def getFile(self,oriFile,destDir,overwrite=True):
        currDir=os.getcwd()
        destDir=os.path.abspath(destDir)
        self.checkAlive()
        if self.connected:        
            self.sftp.chdir(self.home)
            #Check if file exists:
            try:
                self.sftp.stat(oriFile)
            except:
                msg='File not found'
                return (False,msg)
            
            ## Check if directory exists
            try:
                os.chdir(destDir)
            except OSError:
                logging.debug('Destination directory does not exists. Trying to create it')
                fullpath=[]
                head=destDir
                while True:
                    head,tail=os.path.split(head)
                    if tail =='':
                        break
                    fullpath.append(tail)
                
                fullpath.reverse()
                head='/'
                for folder in fullpath:
                    head=os.path.join(head,folder)
                    try:
                        os.chdir(head)
                    except:
                        try:
                            os.mkdir(head)
                        except Exception,e:
                            msg='Error in creation the directory %s: %s'%(head,e)
                            os.chdir(currDir)
                            return (False,msg)
            
            os.chdir(currDir)
            destFile=os.path.join(destDir,os.path.split(oriFile)[1])

            if os.path.isfile(destFile):
                if overwrite:
                    logging.debug('Overwriting existing file: %s'%(destFile))
                    os.remove(destFile)
                else:
                    msg='File %s already exists. overwrite=True if you wanna it.'%(destFile)
                    return (False,msg)
            
            self.sftp.stat(oriFile)
            self.sftp.get(oriFile,destFile)
            permOri=oct(self.sftp.stat(oriFile).st_mode)[-4:]

            os.chmod(destFile, int(permOri,8))

            return (True,destFile)
            
        else:
            msg='Transfer of file %s:%s to %s failed because no connection was possible'%(self.address,oriFile,destDir)
            return (False,msg)


    def close(self):
        if self.ssh.get_transport().is_active():
            self.sftp.close()
            self.ssh.close()
        self.ssh=None
        self.sftp=None
        self.connected=False
        self.info='Disconnected'
 
 
    def fileExists(self,filein):
        self.checkAlive()
        if self.connected:        
            self.sftp.chdir(self.home)
            try:
                self.sftp.stat(filein)
                return True
            except:
                return False
        else:
            logging.error('Connection was not possible for checking file %s:%s'%(filein, self.address))
            return False