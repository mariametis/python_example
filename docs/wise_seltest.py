import pyvo
import pandas as pd
from pandas import DataFrame
import numpy as np
import math

from astroquery.vizier import Vizier
from astropy import units as u
from astroquery.xmatch import XMatch
from astropy.table import Table
from astroquery.utils.tap.core import TapPlus

   
def wise_seltest(): 
   tap_service_url = "https://TAPVizieR.cds.unistra.fr/TAPVizieR/tap"  
   tap_service = pyvo.dal.TAPService(tap_service_url)

   query_name = "pippo"

   frames1 = ()
   
   query="""
         SELECT "II/328/allwise".AllWISE, "II/328/allwise".RAJ2000, "II/328/allwise".DEJ2000, 
         "II/328/allwise".Jmag, "II/328/allwise".e_Jmag, 
         "II/328/allwise".Hmag, "II/328/allwise".e_Hmag, 
         "II/328/allwise".Kmag, "II/328/allwise".e_Kmag,
         "II/328/allwise".W1mag, "II/328/allwise".e_W1mag,
         "II/328/allwise".W2mag, "II/328/allwise".e_W2mag,
         "II/328/allwise".W3mag, "II/328/allwise".e_W3mag,
         "II/328/allwise".W4mag, "II/328/allwise".e_W4mag,
         "II/328/allwise".ccf, "II/328/allwise".ex, "II/328/allwise".var, "II/328/allwise".qph,
         "II/328/allwise".W3mag-"II/328/allwise".W4mag as w3_w4
         FROM "II/328/allwise"
         WHERE ( "II/328/allwise".W3mag-"II/328/allwise".W4mag > -0.2 and "II/328/allwise".W3mag-"II/328/allwise".W4mag < 6.0 and "II/328/allwise".Kmag< 4 )
         """
 

   job = tap_service.submit_job(query,  runid=query_name, queue="2h")
   job.run()

   print('JOB %s: SUBMITTED' % (job.job.runid,))
   print('JOB %s: %s' % (job.job.runid, job.phase))

   #
   # Save the job's url in a file to later retrieve results.
   #
   print('URL: %s' % (job.url,))

   with open('job_url.txt', 'w') as fd:
       fd.write(job.url)



wise_seltest()
