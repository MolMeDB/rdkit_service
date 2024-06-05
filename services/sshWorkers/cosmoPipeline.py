import re, os, time
from services.sshWorkers.lib.file import File
import paramiko
from services.sshWorkers.lib.cuby4 import CUBY4
import requests

CODE_INIT = 0
CODE_OK = 1
CODE_WAITING = 2
CODE_ERR = 3

CODE_EXISTS = 10
CODE_IS_RUNNING = 11
CODE_HAS_RESULT = 12

STEP_INPUT_CHECK = 0
STEP_INIT = 1
STEP_CHECK_RES = 2
STEP_REMOTE_FOLDER = 3
STEP_UPLOAD = 4
STEP_CUBY_CHECK_STATUS = 5
STEP_CUBY_GENERATE = 6
STEP_CUBY_UPLOAD = 7
STEP_CUBY_RUN = 8
STEP_CUBY_RUN_JOB = 81
STEP_COSMO_PREPARE = 9
STEP_COSMO_RUN = 10
STEP_COSMO_DOWNLOAD = 11

sshClient = None

# MOLMEDB_SERVER = "http://192.168.0.183:1515/"
MOLMEDB_SERVER = "https://molmedb.upol.cz/"

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

# paramiko.util.log_to_file("paramiko.log")

def _err(text):
    print("Error:", text)
    raise Exception(text)
def _step(num, text):
    print("STEP", str(num) + ":", text)

def _log(step, code, suffix = ""):
    print("_LOG_: " + str(step) + "/" + str(code) + " [" + suffix + "]")
    
    
########################################
#### Prepare files #####################
########################################
def prepare_files(idFragment, idIon, sshClient, forceRun = False, reRun = False):
    # Check if folder exists
    if not idFragment:
        _log(STEP_INPUT_CHECK, CODE_ERR)
        _err("IdFragment not set.")
        
    if idIon:
        files = requests.get(MOLMEDB_SERVER + "api/conformers/files/" + str(idFragment) + "?idion=" + str(idIon)).json()
    else:
        files = requests.get(MOLMEDB_SERVER + "api/conformers/files/" + str(idFragment)).json()
    
    if len(files) == 0:
        _err("No SDF files found.")
        
    files = [SDFFile(
        fileName=f["name"], 
        folder=f["folder"], 
        content=f["content"], 
        charge=f["charge"], 
        sshClient=sshClient, 
        forceRun=forceRun, 
        reRun=reRun) 
             for f in files]

    # Update base paths
    print("Checking remote folder structure")
    files[0].checkRemoteFolderStructure()
    
    for f in files:
        f.remoteBasePath = files[0].remoteBasePath
    
    return files


########################################
## Force RUN - remote remote data ######
########################################
def remove_remote_folder_structure(idFragment, idIon, sshClient):
    # Check if folder exists
    files = prepare_files(idFragment, idIon, sshClient)

    files[0].removeRemoteFolderStructure()
    return {"status": 'ok'}

#########################################
### Download remote data if exists ######
#########################################
def download_results(idFragment, idIon, sshClient, cosmoType="perm", temp=25, queue="elixir", membraneId=None):
    if not membraneId:
        return {"status": False, "message": "Membrane not set."}
    else:
        membrane = requests.get(MOLMEDB_SERVER + "api/membranes/get/cosmo/" + str(membraneId)).json()
        
    if not membrane["name"]:
        return {"status": False, "message": "Membrane get - invalid response from MolMeDB server."}
    
    files = prepare_files(idFragment, idIon, sshClient)
    
    cosmo = COSMO(
        sshClient=sshClient, 
        forceRun = False, 
        type=cosmoType,
        temperature=temp,
        membrane=membrane,
        membraneName=membrane["name"],
        files=files,
        queue=queue)
    
    # Check status on remote server
    running, hasResult = cosmo.getState()
    
    if not hasResult:
        return {"status": False}
    
    if running:
        return {"status": "running"}
    
    dataPath = cosmo.downloadResults()
    
    if not os.path.exists(dataPath):
        return {"status": False, "message": "Failed to download and archive results."}
    
    return dataPath
    
#########################################
### Upload SDFs to remote server ########
#########################################
def upload_sdf_to_remote_server(idFragment, idIon, sshClient):
    files = prepare_files(idFragment, idIon, sshClient)
    skipped = uploaded = 0
    for f in files:
        r = f.uploadSdfFile()
        if r is True:
            uploaded += 1
        else:
            skipped += 1
            
    return {"status": "ok", "skipped": skipped, "uploaded": uploaded}
    
#########################################
### Upload SDFs to remote server ########
#########################################
def check_optimization_status(idFragment, idIon, sshClient):
    files = prepare_files(idFragment, idIon, sshClient)
    result = {}
    for f in files:
        f.checkOptimizationStatus()
        result[f.name] = {
            "isRunning": f.optimizationRunning,
            "isDone": f.optimizationResults,
            "hasHistory": f.optimizationHistoryRun
        }
            
    return {"status": "ok", "response": result}
  
#########################################
### Run CUBY4 on remote server #########
#########################################          
def run_cuby4(idFragment, idIon, sshClient, cpu=8, ram=32, limitHs=10, queue="elixir", reRun = False):
    files = prepare_files(idFragment, idIon, sshClient, reRun=reRun)
    for f in files:
        if reRun:
            f.clearOptAndCosmoFolders()
        f.checkOptimizationStatus()

    for f in files:
        if not f.optimizationResults and not f.optimizationRunning:
            f.getOptimizeInputs(ncpu=cpu, ram=ram, walltimeHs=limitHs,queue=queue)

    skipped = uploaded = running = 0
    # _log(STEP_CUBY_UPLOAD, CODE_INIT)
    for f in files: 
        if f.optimizationRunning:
            running += 1
        else:
            r = f.uploadCubyFiles()
        
        if r == "skip":
            skipped += 1
        elif not f.optimizationResults:
            uploaded += 1
            
    running = 0
    running_new = 0
    hasResult = 0
    errorFiles = []
    for f in files: 
        if f.optimizationHistoryRun and not f.optimizationResults:
            errorFiles.append(f.fileName)
            continue

        if not f.readyToOptimize:
            if f.optimizationResults:
                hasResult += 1
            else: 
                running += 1
            continue

        jobID = f.runOptimization()
        if jobID:
            running_new += 1
    
    return {
        "status": "ok",
        "running_total": running + running_new,
        "running_new": running_new,
        "uploaded": uploaded,
        "skipped": skipped,
        "hasResult": hasResult,
        "errorFiles": errorFiles  
    }

#########################################
####### Run COSMO #######################
#########################################
def run_cosmo(idFragment, idIon, sshClient, membraneId=None, temp=25, queue="elixir", cosmoType="perm", forceRun = False):
    if not membraneId:
        return {"status": False, "message": "Membrane not set."}
    else:
        membrane = requests.get(MOLMEDB_SERVER + "api/membranes/get/cosmo/" + str(membraneId)).json()
        
    if not membrane["name"]:
        return {"status": False, "message": "Membrane get - invalid response from MolMeDB server."}
        
        
    if cosmoType not in ["perm", "mic"]:
       return {"status": False, "message": "Invalid COSMO type."}
    
    files = prepare_files(idFragment, idIon, sshClient)
    
    cosmo = COSMO(
        sshClient=sshClient, 
        forceRun = forceRun, 
        type=cosmoType,
        temperature=temp,
        membrane=membrane,
        membraneName=membrane["name"],
        files=files,
        queue=queue)
    
    # Check, if SDF are optimized
    optimized = 0
    for f in files:
        isRunning, hasResults = f.checkOptimizationStatus()
        
        if isRunning:
            return {"status": "waiting", "message": "Cuby4 optimization is still running."}
        
        if hasResults:
            optimized += 1
    
    if not optimized or optimized < len(files)/2:
        return {"status": "waiting", "message": "Waiting for at least half of optimized SDF files."}
        
    # Prepare and copy cosmo files
    cosmoFiles = list()
    for f in files:
        if f.optimizationResults:
            f.copyCosmoFile()
            cosmoFiles.append(f)

    cosmo.prepareRun(cosmoFiles)
    # Get COSMO run state
    isRunning, hasResults = cosmo.getState()
    if isRunning:
        return {"status": "ok", "message": "COSMO is already running."}
    if hasResults:
        return {"status": "ok", "message": "COSMO results already exists."}
    # Upload run files
    cosmo.uploadRunFiles()
    # Run cosmo
    jobId = cosmo.run()
    if jobId:
        return {"status": "ok", "message": "COSMO job submitted with number " + str(jobId) + "."}
    

class COSMO:
    QUEUE_ELIXIR = "elixircz@pbs-m1.metacentrum.cz"
    QUEUE_METACENTRUM = "default@pbs-m1.metacentrum.cz	"
    QUEUE_CERIT = "default@pbs-m1.metacentrum.cz	"
    
    def __init__(self, sshClient, forceRun, type, temperature, membrane, membraneName, queue = "elixir", files = []): 
        self.files = files
        self.name = ""
        self.force = forceRun
        self.ssh = sshClient
        self.type = type # perm, mic
        self.temperature = temperature
        self.membraneObject = membrane
        self.membraneName = membraneName
        if "elixir" in str(queue):
            self.queue = self.QUEUE_ELIXIR
        elif "cerit" in str(queue):
            self.queue = self.QUEUE_CERIT
        else:
            self.queue = self.QUEUE_METACENTRUM

        self.cosmoINP = self.cosmoJob = None
        self.cosmoRunning = self.cosmoResults = False
        self.readyToRun = True

    def prepareRun(self, sdfFiles):
        self.files = sdfFiles
        # Prepare cosmo input
        self.cosmoINP = self.genCosmoInput()
        self.cosmoJob = self.genCosmoJob()

    # V2, OK
    def getFolderName(self):
        return self.type + "_" + str(self.membraneName).replace(" ", "-").replace("/", "_") + "_" + str(self.temperature).replace(".", ",")

    # V2, OK
    def getState(self):
        running, hasResult = False, False
        file = self.files[0]
        # Check, if job is potentialy running
        folder = self.getFolderName()
        in_path = str(file.getCosmoInputPath()).rstrip("/") + "/" + folder
        self.ssh.shell_exec("mkdir -p " + in_path)
        # Check if file already exists
        existing = self.ssh.sftp.listdir(in_path)
        if type(existing) is not list:
            existing = list()
        # Checking, if job is not running
        for f in existing:
            f=str(f)
            if f.endswith(".run"):
                running = True
                
        out_path = file.remoteBasePath + "04-COSMO_RESULTS/" + str(self.getFolderName()) + "/"
        self.ssh.shell_exec("mkdir -p " + out_path)
        existing = self.ssh.sftp.listdir(out_path)
        if type(existing) is not list:
            existing = list()
            
        existing = list(filter(lambda f: f.endswith(".tab") or f.endswith(".xml"), existing))

        if len(existing) > 1:
            hasResult = True

        self.cosmoRunning = running
        self.cosmoResults = hasResult
        return running, hasResult

    def uploadRunFiles(self):
        if self.cosmoResults and not self.force:
            self.readyToRun = False
            print(" --- Cosmo already computed. Skipping.")
            return

        if self.cosmoRunning and not self.force:
            self.readyToRun = False
            print(" --- Cosmo looks like running. Skipping.")
            return

        f = self.files[0]
        target = f.getCosmoInputPath()
        sftp = self.ssh.sftp
        folder = self.getFolderName()
        import tempfile

         # Check if file already exists
        self.ssh.shell_exec("mkdir -p " + target)
        existing = sftp.listdir(target)
        if self.membraneName not in existing:
            existing = list()
        else:
            existing = sftp.listdir(target.rstrip("/") + "/" + folder)
            if type(existing) is not list:
                existing = list()

        target = target.rstrip("/") + "/" + folder
        self.ssh.shell_exec("mkdir -p " + target)

        tmpJob = tempfile.NamedTemporaryFile(delete=False)
        tmpInp = tempfile.NamedTemporaryFile(delete=False)
        skipped = False
        try:
            with open(tmpJob.name, "w") as tY:
                tY.write(self.cosmoJob)
            with open(tmpInp.name, "w") as tY:
                tY.write(self.cosmoINP)
            # Upload
            if "cosmo.inp" not in existing or self.force:
                sftp.put(tmpInp.name, target + "/" + "cosmo.inp")
            else: 
                skipped = True
            if "cosmo.job" not in existing or self.force:
                sftp.put(tmpJob.name, target + "/" + "cosmo.job")
            else: 
                skipped = True
            # Copy micelle file
            mic = tempfile.NamedTemporaryFile(delete=False)
            with open(mic.name, "w") as tY:
                tY.write(self.membraneObject["file_content"])
            sftp.put(mic.name, target + "/" + "micelle.mic")
        finally:
            tmpJob.close()
            tmpInp.close()
            os.unlink(tmpJob.name)
            os.unlink(tmpInp.name)

        if skipped:
            return "skip"
        return True

    # V2, OK
    def downloadResults(self):
        if not self.cosmoResults:
            return

        f = self.files[0]
        sftp = self.ssh.sftp

        out_path = f.remoteBasePath + "04-COSMO_RESULTS/" + str(self.getFolderName()) + "/"
        existing = sftp.listdir(out_path)

        if type(existing) is not list or not len(existing):
            _err("Output folder doesnt contain any file...")

        existing = existing if type(existing) is list else list()

        localOutPath = f.folder + "COSMO/" + str(self.getFolderName()) + "/"
        if not os.path.exists(localOutPath):
            os.makedirs(localOutPath)

        for fp in existing:
            p = out_path + fp
            sftp.get(p, localOutPath + fp)
        
        # Make zip    
        import shutil
        shutil.make_archive(f.folder + "archive", "zip", localOutPath)
        
        return f.folder + "archive.zip"

    def run(self):
        if not self.readyToRun:
            return

        f = self.files[0]
        in_folder = f.getCosmoInputPath() + "/" + str(self.getFolderName())
        self.ssh.shell_exec("cd " + in_folder)

        output = self.ssh.shell_exec(self.ssh.qsub + " cosmo.job")
        if len(output) != 1 or not re.match(r"^\d+", output[0]):
            print(output)
            _err("Cannot run the job.", sshClient=self.ssh)
        jobID = re.sub(r'\..*$', "", output[0])
        # Create log to inform, that job was run
        self.ssh.shell_exec("echo 'running' > " + jobID + ".run")
        print(" --- OK. JobID:", jobID)
        return jobID
    
    def genCosmoJob(self):
        f = self.files[0]
        COSMO_INPUT = f.getCosmoInputPath() + "/" + str(self.getFolderName())
        OUTPATH = f.remoteBasePath + "04-COSMO_RESULTS/" + str(self.getFolderName()) + "/"

        return f"""#!/bin/bash
#PBS -q {self.queue}
#PBS -N MMDB_COSMO_{self.type}_{f.name}
#PBS -l select=1:ncpus=10:mem=5gb:scratch_local=5gb
#PBS -l walltime=10:00:00
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR || exit 1

COSMO=/storage/praha5-elixir/home/xjur2k/COSMOlogic/COSMOthermX18/COSMOtherm/BIN-LINUX/cosmotherm
INP={COSMO_INPUT}
OUT={OUTPATH}

cp $INP/cosmo.inp .
cp $INP/micelle.mic .

$COSMO cosmo.inp

rm -f $INP/*.run

mkdir -p $OUT

cp -r * $OUT || export CLEAN_SCRATCH=false
        """

    def genCosmoInput(self):
        content = """ctd = BP_TZVPD_FINE_18.ctd cdir = "/storage/praha5-elixir/home/xjur2k/COSMOlogic/COSMOthermX18/COSMOtherm/CTDATA-FILES" ldir = "/storage/praha5-elixir/home/xjur2k/COSMOlogic/COSMOthermX18/licensefiles"
rmic=micelle.mic"""
        if self.type != "mic": # Cosmo perm
            content += " unit notempty wtln ehfile"

        l = self.files.copy()
        first = l.pop(0)
        content += " accc \n! " + first.name + " conformer computation !\n"
        # Add files

        self.name = first.name

        content += "f = " + str(first.name) + ".ccf fdir=" + str(first.getCosmoInputPath()) + " Comp = " + str(first.name)

        if len(l):
            content += " ["
            for f in l:
                content += "\n" + "f = " + str(f.name) + ".ccf fdir=" + str(f.getCosmoInputPath())
            content += " ]"

        content += "\n"

        if self.type == "mic":
            content += "tc=" + str(self.temperature) + " x_pure=Micelle"
        else:
            content += "tc=" + str(self.temperature) + " micelle permeability centersig2 rmic=micelle.mic"

        return content

class SDFFile:
    # Constants
    SERVER_FOLDER_PREFIX = "~/.MolMeDB/COSMO/"
    LAST_FOLDER = "conformers"
    FILE_FOLDERS = [
        "01-INPUT",
        "02-OPTIMIZE",
        "03-COSMO_INPUT",
        "04-COSMO_RESULTS"
    ]
    
    SCRIPTPATH_PS = "[SCRIPT_PATH]"
    SDFPATH_PS = "[SDFPATH_PS]"
    LOGPATH_PS = "[LOGPATH_PS]"
    

    def __init__(self, fileName, folder, content, sshClient=None, charge = 0, forceRun=False, reRun = False):
        # self.path = str(path).strip()
        self.sshClient = sshClient
        self.content = content
        self.fileName = fileName
        self.folder = folder
        self.remoteBasePath = None
        self.forceRun = forceRun
        self.name = re.sub(r"\.sdf", "", self.fileName)
        self.cuby = CUBY4()
        self.readyToOptimize = False
        self.charge = charge
        self.reRun = reRun

        self.optimizeScripts = None

        self.optimizationRunning = self.optimizationResults = self.optimizationHistoryRun = False

        if not fileName or not content:
            _err("Target file not exists.")

    def getCosmoInputPath(self):
        return self.remoteBasePath + self.FILE_FOLDERS[2]

    # def getFolder(self, path):
    #     path = str(path)
    #     # Check if contains file with suffix
    #     lastPos = path.rfind("/")
    #     if path.rfind(".") > lastPos:
    #         file = path[lastPos+1:].strip()
    #         path = path[:lastPos] + "/"
    #     else:
    #         file = None
    #         path = path.strip("/") + "/"
    #     return path, file

    def removeRemoteFolderStructure(self):
        path = str(self.folder)
        # index = path.find(self.LAST_FOLDER)
        # if not index:
        #     _log(STEP_REMOTE_FOLDER, CODE_ERR)
        #     _err("Invalid input folder structure.")
        # path = path[index+(self.LAST_FOLDER.__len__())+1:]
        path = path.strip("/") + "/"
        if not path:
            _log(STEP_REMOTE_FOLDER, CODE_ERR)
            _err("Invalid input folder structure.")
        # Remove file structure
        path = self.SERVER_FOLDER_PREFIX + path
        self.sshClient.shell_exec('rm -rf ' + path)
        
    def clearOptAndCosmoFolders(self):
        # Check, if job is potentialy running
        out_folder = self.remoteBasePath + self.FILE_FOLDERS[1] + "/" + self.name
        self.sshClient.shell_exec("mkdir -p " + out_folder)
        # Check if file already exists
        existing = self.sshClient.sftp.listdir(out_folder)
        # Check, if results exists
        if type(existing) is not list:
            existing = list()
        # Checking, if job has results
        for f in existing:
            f=str(f)
            if f == "OUTPUT":
                # Check if content is valid
                l = self.sshClient.sftp.listdir(out_folder + "/OUTPUT")
                l = list(filter(lambda f: f.startswith("job_step"), l))
                l.sort()
                if not len(l):
                    break
                else:
                    l = l[-1] # Get last result folder
                    l = self.sshClient.sftp.listdir(out_folder + "/OUTPUT/" + l)
                    if "out.ccf" in l:
                        return True
        
        # No results? Clear
        self.sshClient.shell_exec("rm -rf " + out_folder)

    def checkRemoteFolderStructure(self):
        try:
            path = str(self.folder)
            # Make file structure
            path = self.SERVER_FOLDER_PREFIX + path
            cmd = ""
            for fold in self.FILE_FOLDERS:
                if len(cmd):
                    cmd += " && "
                cmd += "mkdir -p " + path + fold
            response = self.sshClient.shell_exec(cmd)
            if len(response) != 0:
                return False, response.join("\n")
            # Save remote base path for next usage
            print("Readlink")
            out = self.sshClient.shell_exec("cd " + path + " && pwd", sleepTime=2)

            if out and len(out) and str(out[-1]).startswith("/"):
                self.remoteBasePath = str(out[-1].rstrip("/")) + "/"
            else:
                print(path, out)
                _err("Cannot obtain full remote base path.")
            return True
        except Exception as e:
            print(e)
            _err("Exception occured during making file structure.")

    def uploadSdfFile(self):
        import tempfile
        try:
            if self.remoteBasePath == "":
                _err("Remote folder seems to not exist.")

            target = str(self.remoteBasePath) + self.FILE_FOLDERS[0] + "/"
            # Check if file already exists
            existing = self.sshClient.sftp.listdir(target)
            if existing and len(existing) and self.fileName in existing and not self.forceRun:
                return None
            print(" - Uploading ", self.name , "file.")
            ## Create temporary file and upload
            tmp = tempfile.NamedTemporaryFile(delete=False)
            with open(tmp.name, "w") as fl:
                fl.write(self.content)
            self.sshClient.sftp.put(tmp.name, target + self.fileName, confirm=False)
        except Exception as e:
            _log(STEP_UPLOAD, CODE_ERR)
            _err("Exception occured during uploading files..")
        return True

    def checkOptimizationStatus(self):
        running, hasResult, wasRun = False, False, False
        # Check, if job is potentialy running
        out_folder = self.remoteBasePath + self.FILE_FOLDERS[1] + "/" + self.name
        self.sshClient.shell_exec("mkdir -p " + out_folder)
        # Check if file already exists
        existing = self.sshClient.sftp.listdir(out_folder)
        if type(existing) is not list:
            existing = list()
        # Checking, if job is not running
        for f in existing:
            f=str(f)
            if f.endswith(".run"):
                running = True
            if f == "OUTPUT":
                wasRun = True
                # Check if content is valid
                l = self.sshClient.sftp.listdir(out_folder + "/OUTPUT")
                l = list(filter(lambda f: f.startswith("job_step"), l))
                l.sort()
                if not len(l):
                    hasResult = False
                else:
                    l = l[-1] # Get last result folder
                    l = self.sshClient.sftp.listdir(out_folder + "/OUTPUT/" + l)
                    if "out.ccf" in l:
                        hasResult = True

        self.optimizationRunning = running
        self.optimizationResults = hasResult
        self.optimizationHistoryRun = wasRun

        return running, hasResult


    def getOptimizeInputs(self, ncpu=8, ram=32, walltimeHs=10,queue="elixir"):
        self.cuby.generate(self, ncpu=ncpu, ram=ram, walltimeHs=walltimeHs, queue=queue, includeTurobmolePreComputation=False)

    def uploadCubyFiles(self):
        if self.optimizationResults and not self.forceRun:
            self.readyToOptimize = False
            print(" --- Job", self.name, "already computed. Skipping.")
            return

        if not self.cuby.filesContent:
            _err("Nothing to upload for file" + self.name + ".")

        import tempfile
        print(" - Uploading optimizing job files for", self.name)
        # Create folder for each structure
        out_folder = self.remoteBasePath + self.FILE_FOLDERS[1] + "/" + self.name
        self.sshClient.shell_exec("mkdir -p " + out_folder)
        sdf_folder = self.remoteBasePath + self.FILE_FOLDERS[0]
        log_foler = self.remoteBasePath.split("/COSMO/")[0] + "/COSMO/LOGS/RUNS/"
        
        # Check if file already exists
        existing = self.sshClient.sftp.listdir(out_folder)
        if type(existing) is not list:
            existing = list()

        # Checking, if job is not running
        if self.optimizationRunning and not self.forceRun and not self.reRun:
            print(" --- Job", self.name, "looks like running. Skipping.")
            self.readyToOptimize = False
            return "running"
        # Fill variables with paths
        contentYAML = str(self.cuby.getYaml())
        contentjob = str(self.cuby.getJob())
        contentjob = contentjob.replace(self.SCRIPTPATH_PS, out_folder)
        contentjob = contentjob.replace(self.SDFPATH_PS, sdf_folder)
        contentjob = contentjob.replace(self.LOGPATH_PS, log_foler)
        # Upload both to the remote server
        # Create temporary file and upload to server
        tmpYaml = tempfile.NamedTemporaryFile(delete=False)
        tmpJob = tempfile.NamedTemporaryFile(delete=False)
        skipped = False
        self.readyToOptimize = True
        try:
            with open(tmpYaml.name, "w") as tY:
                tY.write(contentYAML)
            with open(tmpJob.name, "w") as tY:
                tY.write(contentjob)
            # Upload
            if self.name + ".yaml" not in existing or self.forceRun or self.reRun:
                self.sshClient.sftp.put(tmpYaml.name, out_folder + "/" + self.name + ".yaml")
            else: 
                skipped = True
            if self.name + ".job" not in existing or self.forceRun or self.reRun:
                self.sshClient.sftp.put(tmpJob.name, out_folder + "/" + self.name + ".job")
            else: 
                skipped = True
        finally:
            tmpJob.close()
            tmpYaml.close()
            os.unlink(tmpJob.name)
            os.unlink(tmpYaml.name)

        self.sshClient.shell_exec("cd ~/")

        if skipped:
            return "skip"
        return True

    def runOptimization(self):
        try:
            print(" - Trying to run:", self.name)
            out_folder = self.remoteBasePath + self.FILE_FOLDERS[1] + "/" + self.name
            self.sshClient.shell_exec("cd " + out_folder)
            # Remove output folder if present
            self.sshClient.shell_exec("rm -r OUTPUT")
            # Run new job
            scriptFile = self.name + ".job"
            output = self.sshClient.shell_exec(self.sshClient.qsub + " " + scriptFile, sleepTime=3)
            if not len(output) or not re.match(r"^\d+", output[-1]):
                print(output)
                _err("Cannot run the job.")
            jobID = re.sub(r'\..*$', "", output[-1])
            # Create log to inform, that job was run
            self.sshClient.shell_exec("echo 'running' > " + jobID + ".run")
            self.optimizeJobID = jobID
            print(" --- OK. JobID:", jobID)
            return jobID
        except Exception as e: 
            _log(STEP_CUBY_RUN, CODE_ERR)
            print(e)
            _err("Exception occured.")

    def copyCosmoFile(self):
        # Check, if job is potentialy running
        out_folder = self.remoteBasePath + self.FILE_FOLDERS[1] + "/" + self.name + "/OUTPUT"
        existing = self.sshClient.sftp.listdir(out_folder)
        if type(existing) is not list:
            existing = list()
        # Check if content is valid
        l = self.sshClient.sftp.listdir(out_folder)
        l = list(filter(lambda f: f.startswith("job_step"), l))
        l.sort()
        if not len(l):
            _err("Cannot find cosmo result file for structure " + self.name + ".")
        else:
            l = l[-1] # Get last result folder
            l1 = self.sshClient.sftp.listdir(out_folder + "/" + l)
            if "out.ccf" not in l1:
                _err("Cannot find cosmo result file for structure " + self.name + ".")
            out = self.sshClient.shell_exec("cp " + out_folder + "/" + l + "/out.ccf " + self.remoteBasePath + self.FILE_FOLDERS[2] + "/" + self.name + ".ccf")
            if len(out) and str(out[0]).startswith(""):
                print(out)
                _err("Cannot copy result file.")



class SSHClient:
    def __init__(self, username, password, host, port=22):
        self.username = username
        self.password = password
        self.host = host
        self.shell = None
        self.qsub = None
        self.ssh = None
        self.sftp = None

        self.connect()

    def connect(self):
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        # Connect
        self.ssh.connect(
            self.host, 
            username=self.username, 
            password=self.password
        )
        import time
        # init shell
        self.shell = self.ssh.invoke_shell()
        self.shell_exec("clear")
        # Switch to bash
        self.shell_exec("bash")
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
        self.sftp = self.ssh.open_sftp()

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

    def shell_exec(self, command, clear=True):
        if clear:
            pass
            # self.shell.send("clear\n")
        # Execute command
        self.shell.send(command+"\n")
        # Read output
        out = ""
        # sleep is essential, recv_ready returns False without sleep
        time.sleep(1)
        while self.shell.recv_ready():
            out += str(self.shell.recv(2048))

        out = out.replace("\\r", "")
        out = out.replace("\\t", "")
        out = out.replace("\\x1b", "\x1b")
        out = out.split('\\n')
        if len(out) < 2:
            return []
        out = [ansi_escape.sub('', t) for t in out]
        if "password" in str(out[-1]).lower():
            return [str(out[-1])]
        out.pop()

        if len(command) > 10:
            command = command[:10]

        if len(out) and str(command).lower() in str(out[0]).lower():
            out.pop(0)
        return out

def get_conformer_folder(id_fragment, id_ion):
    group = int(id_fragment / 10000)
    group = group * 10000
    group = str(group) + "-" + str(group+10000)
    return "media/files/conformers/" + group + "/" + str(id_fragment) + "/" + str(id_ion) + "/"

