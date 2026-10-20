# This file takes all the jobs in jobs.py and runs them locally 

from jobs import jobs
from make_jobs import make_job
from datetime import datetime
import subprocess

# Loop through jobs
for job in jobs:
    # Convert to command
    command = make_job(job)

    # Print start time and run
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}]")
    print(f"Starting: L={job["L"]}, beta={job["beta"]} J1={job["J_1"]}, J2={job["J_2"]}, basis={job["basis"]}")
    subprocess.run(command, check=True)
# Print end time
print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}]")
print("End")

