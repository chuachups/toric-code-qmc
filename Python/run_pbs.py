# This file takes all the jobs in jobs.py and creates PBS files for each

from jobs import jobs
from make_jobs import make_job
import subprocess
import shlex

# Loop through jobs
for job in jobs:
    # Convert to command
    command = make_job(job)

    # Construct the PBS file
    pbs_filename = "PBS/" + command[3].replace("RawData/", "").replace("/","__").replace("____","__").replace(".h5", ".pbs")
    with open(pbs_filename, "w") as file:
        file.write("#!/bin/bash" + "\n")
        file.write("#PBS -P CFP05-SF-050" + "\n")
        file.write("#PBS -N " + pbs_filename.replace("PBS/", "").replace(".pbs", "") + "\n")
        file.write("#PBS -q serial" + "\n")
        file.write("#PBS -l select=1:ncpus=1" + "\n")
        file.write("#PBS -j oe" + "\n")
        file.write("cd ${PBS_O_WORKDIR}" + "\n")
        file.write(shlex.join(command) + "\n")

    # Call qsub
    result = subprocess.run(["qsub", pbs_filename], check=True)
    print(f"{pbs_filename}: {result.stdout.strip()}")

### PBS Script Format:
###
### #!/bin/bash
### #PBS -P <project name>
### #PBS -N <job title>
### #PBS -l select=1:ncpus=1
### #PBS -j oe
### 
### cd ${PBS_O_WORKDIR};
### <Julia command>