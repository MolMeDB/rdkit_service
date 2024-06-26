import re

QUEUE_ELIXIR = "elixircz@pbs-m1.metacentrum.cz"
QUEUE_METACENTRUM = "default@pbs-m1.metacentrum.cz"
QUEUE_CERIT = "default@pbs-m1.metacentrum.cz"
SCRIPTPATH_PS = "[SCRIPT_PATH]"
SDFPATH_PS = "[SDFPATH_PS]"
LOGPATH_PS = "[LOGPATH_PS]"
JOB_NAME_PREFIX = "MMDB_C_"

class CUBY4:
    def __init__(self): 
        self.filesContent = {}

    def returnTemplates(self):
        return self.filesContent

    def getYaml(self):
        return self.filesContent["yaml"]

    def getJob(self):
        return self.filesContent["job"]

    def generate(self, file, ncpu = 8, ram=32, walltimeHs = 10, queue = "elixir", includeTurobmolePreComputation=False):
        self.filesContent = {
            "yaml": self.generateYAML(file.fileName, charge=file.charge, ncpu=ncpu, ram=ram, includePreComp=includeTurobmolePreComputation),
            "job": self.generateJob(fileName=file.fileName, ncpu=ncpu, ram=ram, walltimeHs=walltimeHs, queue=queue) 
        }

    def generateYAML(self, fileName, charge=0, ncpu = 8, ram=32, includePreComp = False):
        if includePreComp:
            content = f"""# Calculation
job: multistep
steps: gopt, genergy, wopt, wcosmo

calculation_common:
  interface: turbomole
  method: dft
  dft_grid_custom: m3
  functional: b-p
  optimizer: lbfgs
  basisset: def-TZVP
  charge: {charge}
  geometry: {fileName}
  job_cleanup: no
  mem: {ram}000
  delete_large_files: yes
  cuby_threads: {ncpu}
  parallel: {ncpu}
  parallel_mode: shm

calculation_gopt:
  job: optimize
  opt_quality: 0.1
  maxcycles: 2000
  density_convergence: 7

calculation_genergy:
  job: energy
  density_convergence: 6
  basisset: def2-TZVPD

calculation_wopt:
  job: optimize
  maxcycles: 2000
  opt_quality: 0.1
  density_convergence: 7
  solvent_model: cosmo

calculation_wcosmo:
  job: energy
  density_convergence: 7
  basisset: def2-TZVPD
  solvent_model: cosmo
"""
        else:
            content = f"""# Calculation
job: multistep
steps: wopt, wcosmo

calculation_common:
  interface: turbomole
  method: dft
  dft_grid_custom: m3
  functional: b-p
  optimizer: lbfgs
  basisset: def-TZVP
  charge: {charge}
  geometry: {fileName}
  job_cleanup: no
  mem: {ram}000
  delete_large_files: yes
  cuby_threads: {ncpu}
  parallel: {ncpu}
  parallel_mode: shm

calculation_wopt:
  job: optimize
  maxcycles: 2000
  opt_quality: 0.1
  density_convergence: 7
  solvent_model: cosmo

calculation_wcosmo:
  job: energy
  density_convergence: 7
  basisset: def2-TZVPD
  solvent_model: cosmo
""" 
        return content

    def generateJob(self, fileName, ncpu = 8, ram=32, walltimeHs = 10, queue = "elixir"):
        name = re.sub(r"\.sdf", "", fileName) # Remove file suffix
        if "elixir" in queue.lower():
            queue = QUEUE_ELIXIR
        else:
            queue = QUEUE_METACENTRUM
        
        content = f"""#!/bin/bash
#PBS -q {queue} 
#PBS -l select=1:ncpus={ncpu}:mem={ram}gb:scratch_shm=true
#PBS -l walltime={walltimeHs}:00:00
#PBS -o {LOGPATH_PS}{name}.out
#PBS -e {LOGPATH_PS}{name}.err
#PBS -N {JOB_NAME_PREFIX}OPT_{name}
#PBS -m n
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR || exit 1

module add ruby/3.1.0-gcc-10.2.1-ysif6sl
module add turbomole-7.6-smp
module add turbomole-7.6-mpi
module add turbomole-7.6

SDF_WORKDIR={SDFPATH_PS}
SCRIPT_WORKDIR={SCRIPTPATH_PS}
OUTPUT_WORKDIR=$SCRIPT_WORKDIR/OUTPUT

# Adds cuby as executable command
export PATH=/storage/praha5-elixir/home/xjur2k/SOFTWARE/cuby4_new/cuby4:$PATH

# Copy SDF
cp $SDF_WORKDIR/{fileName} .
# Copy yaml
cp $SCRIPT_WORKDIR/{name}.yaml .

# Run calculation
cuby4 {name}.yaml &> LOG

# Check if exists result OUTPUT folder
mkdir -p $OUTPUT_WORKDIR

# Rm job id
rm -f $SCRIPT_WORKDIR/*.run

# Copy result to OUTPUT directory
cp -r * $OUTPUT_WORKDIR/
        """

        return content