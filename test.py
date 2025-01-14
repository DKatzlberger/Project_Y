import os
import sys
import stat
import shutil
import subprocess
import re
import time

def wait_for_job(job_id, poll_interval=10):
    # Poll qstat until the job ID is no longer listed
    while True:
        try:
            result = subprocess.run(["qstat"], capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"qstat command failed: {result.stderr}")
            # Check if job ID is in qstat output
            if job_id not in result.stdout:
                print(f"Job {job_id} has completed.")
                break
        except Exception as e:
            print(f"Error checking job status: {e}")
        time.sleep(poll_interval)


# TESTING
job_id = "PanCancerAtlas"
result = subprocess.run(["qstat"], capture_output=True, text=True)

result.stdout.splitlines()[2]