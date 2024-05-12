import re, os, time, binascii
import paramiko

CODE_INIT = 0
CODE_OK = 1
CODE_WAITING = 2
CODE_ERR = 3

CODE_EXISTS = 10
CODE_IS_RUNNING = 11
CODE_HAS_RESULT = 12

STEP_INPUT_CHECK = 0
STEP_INIT = 1

# Global variables
param_server = False
param_username = False
param_password = False
sshClient = False

ansi_escape = re.compile(r'''
    \x1B  # ESC
    (?:   # 7-bit C1 Fe (except CSI)
        [@-Z\\-_]
    |     # or [ for CSI, followed by a control sequence
        \[
        [0-?]*  # Parameter bytes
        [ -/]*  # Intermediate bytes
        [@-~]   # Final byte
    )
''', re.VERBOSE)

def _err(text, sshClient = None):
    import sys
    print("Error:", text)
    if sshClient: sshClient.close()
    sys.exit(1)

def _step(num, text):
    print("STEP", str(num) + ":", text)

def _log(step, code, suffix = ""):
    print("_LOG_: " + str(step) + "/" + str(code) + " [" + suffix + "]")

class SSHClient:
    def __init__(self, username, password, host, port=22, ignoreSFTP = False):
        self.username = username
        self.password = password
        self.host = host
        self.shell = None
        self.qsub = None
        self.ssh = None
        self.sftp = True if not ignoreSFTP else False
        self.connect()
        
    def isConnected(self):
        if not self.ssh or self.ssh.get_transport() is None:
            return False
        if not self.ssh.get_transport().is_active():
            return False
        try:
            transport = self.ssh.get_transport()
            transport.send_ignore()
            return True
        except EOFError as e:
            return False
        
    def reconnect(self):
        if self.isConnected():
            return True
        self.connect()

    def connect(self):
        print("Connecting...")
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        # Connect
        self.ssh.connect(
            self.host, 
            username=self.username, 
            password=self.password
        )
        # init shell
        self.shell = self.ssh.invoke_shell()
        self.shell_exec("clear")
        # Switch to bash
#        self.shell_exec("bash")
        if not self.kinit():
            self.close()
            _err("Cannot initialize kerberos token.")
        
        # Get qsub path
        for i in range(3):
            output = self.shell_exec("which qsub")
            if i < 2 and (not len(output) or not str(output[0]).startswith("/")):
                time.sleep(2)

        if not len(output) or not str(output[0]).startswith("/"):
            self.close()
            _err("Cannot get `qsub` remote path")
        self.qsub = output[0]
        if self.sftp:
            self.sftp = self.ssh.open_sftp()
        else:
            self.sftp = None

    def kinit(self):
        out = self.shell_exec("kinit")
        time.sleep(4)
        if not len(out):
            return True
        # Set password if required
        if "password" in out[0].lower():
            out = self.shell_exec(self.password, clear=False)
        return True

    def close(self):
        if self.sftp: self.sftp.close()
        if self.ssh: self.ssh.close()

    def shell_exec(self, command, clear=False, sleepTime=1):
        if clear:
            pass
#            self.shell.send("clear\n")
        # Execute command
        self.shell.send(command+"\n")
        # Read output
        out = ""
        # sleep is essential, recv_ready returns False without sleep
        time.sleep(sleepTime)
        while self.shell.recv_ready():
            out += str(self.shell.recv(2048))

        out = out.replace("\\r", "")
        out = out.replace("\\t", "")
        out = out.replace("\\x1b", "\x1b")
        out = out.replace("\\x00", "")
        out = re.sub(r'\s+', " ", out)
        out = out.split('\\n')

        if len(out) < 2:
            return []
        out = [ansi_escape.sub('', t) for t in out]
        if "password" in str(out[-1]).lower():
            return [str(out[-1])]
        if (self.username + "@") in str(out[-1]).lower():
            out.pop()

        if len(command) > 10:
            command = command[:10]

        if len(out) and str(command).lower() in str(out[0]).lower():
            out.pop(0)
        return out
    

def connect(server, username, password, ignoreSFTP = False):
    if(not server or not username or not password):
        return False
    
    global param_username, param_server, param_password, sshClient
    
    # If params changed or server is disconnected, reconnect
    if (server != param_server or username != param_username or password != param_password or not sshClient or not sshClient.isConnected()):
        if isinstance(sshClient, SSHClient) and sshClient.isConnected():
            sshClient.close()
        sshClient = SSHClient(username=username, password=password,host=server, port=22, ignoreSFTP=ignoreSFTP)
        param_server = server
        param_username = username
        param_password = password
    return True


############################
######## ENDPOINTS #########
############################
def getRunningJobs(params):
    from services.sshWorkers import getRunningJobs as GRJ
    global sshClient, param_username
    sshClient.reconnect()
    
    if(not sshClient or not sshClient.isConnected()):
        return {"status": False, "message": "SSHclient not connected."}
    
    return GRJ.run(
        sshClient=sshClient, 
        username=param_username, 
        include_finished=params["include_finished"] if "include_finished" in params else False)


##################################################
### Removes remote files - compute cosmo results #
### Files are zipped and returned as attachment ##
##################################################
def cosmo_download_results(parameters):
    # check parameters
    required = {
        "idFragment": None,
        "idIon": None,
        "cosmoType": None,
        "temp": None,
        "membraneId": None
    }
    
    for k in required:
        if not k in parameters:
            return {"status": False, "message": "Missing parameter '" + k + "'."}
        
    from services.sshWorkers import cosmoPipeline
    global sshClient, param_username
    sshClient.reconnect()
    
    if(not sshClient or not sshClient.isConnected()):
        return {"status": False, "message": "SSHclient not connected."}
    
    return cosmoPipeline.download_results(
        parameters["idFragment"], 
        parameters["idIon"], 
        sshClient=sshClient, 
        cosmoType=parameters["cosmoType"], 
        temp=parameters["temp"], 
        membraneId=parameters["membraneId"])
    
#############################################
### Removes remote folder structure #########
#############################################
def cosmo_clear_folder_structure(parameters):
    # check parameters
    required = {
        "idFragment": None
    }
    
    for k in required:
        if not k in parameters:
            return {"status": False, "message": "Missing parameter '" + k + "'."}
    
    from services.sshWorkers import cosmoPipeline
    global sshClient, param_username
    sshClient.reconnect()
    
    if(not sshClient or not sshClient.isConnected()):
        return {"status": False, "message": "SSHclient not connected."}
    
    return cosmoPipeline.remove_remote_folder_structure(
        parameters["idFragment"], 
        parameters["idIon"] if "idIon" in parameters else None,  
        sshClient=sshClient)
    
#############################################
### Uploads SDF FILES to remote server ######
#############################################
def cosmo_upload_sdf(parameters):
    # check parameters
    required = {
        "idFragment": None,
        "idIon": None,
    }
    
    for k in required:
        if not k in parameters:
            return {"status": False, "message": "Missing parameter '" + k + "'."}
    
    from services.sshWorkers import cosmoPipeline
    global sshClient, param_username
    sshClient.reconnect()
    
    if(not sshClient or not sshClient.isConnected()):
        return {"status": False, "message": "SSHclient not connected."}
    
    return cosmoPipeline.upload_sdf_to_remote_server(parameters["idFragment"], parameters["idIon"],  sshClient=sshClient)

    
###################################################
### Checks optimization status for each file ######
###################################################
def cosmo_check_optimization_status(parameters):
    # check parameters
    required = {
        "idFragment": None,
        "idIon": None,
    }
    
    for k in required:
        if not k in parameters:
            return {"status": False, "message": "Missing parameter '" + k + "'."}
    
    from services.sshWorkers import cosmoPipeline
    global sshClient, param_username
    sshClient.reconnect()
    
    if(not sshClient or not sshClient.isConnected()):
        return {"status": False, "message": "SSHclient not connected."}
    
    return cosmoPipeline.check_optimization_status(parameters["idFragment"], parameters["idIon"],  sshClient=sshClient)

    
###################################################
### Run CUBY4 on remote server for each file ######
###################################################
def cosmo_optimize_sdf(parameters):
    # check parameters
    required = {
        "idFragment": None,
        "idIon": None,
        "queue": None
    }
    
    for k in required:
        if not k in parameters:
            return {"status": False, "message": "Missing parameter '" + k + "'."}
    
    from services.sshWorkers import cosmoPipeline
    global sshClient, param_username
    sshClient.reconnect()
    
    if(not sshClient or not sshClient.isConnected()):
        return {"status": False, "message": "SSHclient not connected."}
    
    return cosmoPipeline.run_cuby4(
        parameters["idFragment"], 
        parameters["idIon"],  
        sshClient=sshClient, 
        queue=parameters["queue"], 
        reRun=parameters["reRun"] if "reRun" in parameters else False,
        limitHs=parameters["hours"] if "hours" in parameters else 20)

#############################################
### Run COSMO on remote server ##############
#############################################
def runCosmo(parameters):
    required = {
        "idFragment": None,
        "idIon": None,
        "membraneId": None,
        "temp": None,
        "cosmoType": None,
        "queue": None
    }
    
    for k in required:
        if not k in parameters:
            return {"status": False, "message": "Missing parameter '" + k + "'."}
    
    from services.sshWorkers import cosmoPipeline
    global sshClient, param_username
    sshClient.reconnect()
    
    if(not sshClient or not sshClient.isConnected()):
        return {"status": False, "message": "SSHclient not connected."}
    
    return cosmoPipeline.run_cosmo(
        parameters["idFragment"],
        parameters["idIon"],
        sshClient=sshClient,
        membraneId=parameters["membraneId"],
        temp=parameters["temp"],
        cosmoType=parameters["cosmoType"],
        queue=parameters["queue"],
        forceRun=parameters["forceRun"] if "forceRun" in parameters else False
    )


