import pyvo
import pandas as pd
from pandas import DataFrame
import numpy as np


   

def retrieve(filemy): 

#
# Recreate the job from url and session (token)
#

# read the url
    with open('job_url.txt', 'r') as fd:
       job_url = fd.readline()

    # recreate the job 
    job = pyvo.dal.AsyncTAPJob(job_url)

    #
    # Check the job status
    #
    print('JOB %s: %s' % (job.job.runid, job.phase))

    # if still running --> exit
    if job.phase not in ("COMPLETED", "ERROR", "ABORTED"):
       exit(0)

    #
    # Fetch the results
    #
    job.raise_if_error()
    print('\nfetching the results...')
    tap_results = job.fetch_result()
    print('\n...DONE\n')
    frames1 = (tap_results.to_table().to_pandas())  
    print(frames1)
    df=DataFrame(frames1)
    df.to_csv(filemy)



retrieve('sel_wise_brightTEST.csv')
