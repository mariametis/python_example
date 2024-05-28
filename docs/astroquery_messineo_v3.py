from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
from astroquery.gaia import GaiaClass
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.xmatch import XMatch
from astropy.table import Table

import pandas as pd
from pandas import DataFrame

import numpy as np
from numpy import array

import pyvo

from gaiaxpy import calibrate
from gaiaxpy import plot_spectra

# Hide warnings
import warnings
warnings.catch_warnings()
warnings.simplefilter("ignore")
#%matplotlib inline
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
from pylab import *
import datetime
import pathlib




def get_twomassid(mycatalog,outputc_csv): 
    suarez = Vizier(catalog=mycatalog,columns=["*", 'RAJ2000', 'DEJ2000'],row_limit=500,vizier_server='vizier.cds.unistra.fr').query_constraints(RAJ2000=">0 ")[0]
    print(suarez)                       

    frames1 =suarez.to_pandas()
    c=SkyCoord(frames1['RAJ2000'],frames1['DEJ2000'], unit=(u.hourangle, u.deg))
    frames1['RA_deg']=c.ra.degree
    frames1['DEC_deg']=c.dec.degree
    frames1.to_csv('frames1.csv')

    #XMatch only works with coordinates in degrees
    input_table = Table.read('frames1.csv')
    tabletwo = XMatch.query(cat1=input_table,
                     cat2='vizier:II/246/out',
                     max_distance=2.5 * u.arcsec, colRA1='RA_deg',
                     colDec1='DEC_deg')
    frames2 =tabletwo.to_pandas()
    frames2.to_csv(outputc_csv)



def get_gaiaid(input_csv,outputc_csv): 
    tap_service_url = "https://gaia.ari.uni-heidelberg.de/tap"
    tap_service = pyvo.dal.TAPService(tap_service_url)

    datacsv = pd.read_csv(input_csv)
    data = []
    frames = ()
    id3=datacsv['2MASS']
    nlen=len(id3)


    for  i in range(1, nlen):
       nametwo = str(id3[i])
# Submit queries
       query=f"SELECT source_id,  clean_tmass_psc_xsc_oid, original_ext_source_id  FROM gaiadr3.tmass_psc_xsc_best_neighbour    WHERE original_ext_source_id LIKE '{nametwo}' "
       print (query)
       if nametwo != 'none': 
          tap_result = tap_service.run_sync(query)                 
       if nametwo != 'none': 
          frames =frames + (tap_result.to_table().to_pandas(),) #turple
#
# Contatenate into a pandas.DataFrame
#
    df_results = pd.concat(frames)
    df_results.head()
    print (type(df_results))  #<class 'pandas.core.frame.DataFrame'>

    df3=df_results.set_index('original_ext_source_id').join(datacsv.set_index('2MASS'),how='inner')
# Pandas join on column
    print(df3)
    
    df3.to_csv(outputc_csv)



def get_gaia_par(input_csv,outputc_csv): 
    #Authenticated access (TAP+ only)

    Gaia.login_gui()
    job = Gaia.upload_table(upload_resource=input_csv,
                            table_name="user_table100", format="csv")

    query=f"SELECT  gaia.source_id,gaia.phot_g_mean_mag,gaia.phot_bp_mean_mag,phot_rp_mean_mag,bp_rp,phot_variable_flag,non_single_star,has_xp_continuous,phot_bp_n_obs,phot_rp_n_obs,ra,dec \
    FROM gaiadr3.gaia_source AS gaia \
    JOIN user_mmessi01.user_table100 AS ul \
    ON gaia.source_id = ul.source_id \
    WHERE ( gaia.source_id = ul.source_id)"
    print (query)


    job = Gaia.launch_job(query=query)
    results = job.get_results()

    frames =(results.to_pandas(),)  

    #
    # Contatenate into a pandas.DataFrame
    #
    df_results = pd.concat(frames)
    df_results.head()
    print (df_results)
    
    datacsv = pd.read_csv(input_csv)
    df3=df_results.set_index('source_id').join(datacsv.set_index('source_id'),how='inner')
    # Pandas join on column
    print(df3)
    df3.to_csv(outputc_csv)

    job = Gaia.delete_user_table("user_table100")
    
def get_bprp(input_csv,dir): 

# Connect to Gaia archive
  gaia = GaiaClass(gaia_tap_server='https://gea.esac.esa.int/', gaia_data_server='https://gea.esac.esa.int/')
  #gaia.login()

  datacsv = pd.read_csv(input_csv)
  data = []
  frames = ()
  id3=datacsv['source_id']
  continuous=datacsv['has_xp_continuous']
  nlen=len(id3)

  for  i in range(0, nlen):
    namedr3 = str(id3[i])
    has_cont= str(continuous[i])
    #print(i, namedr3,has_cont)

    if namedr3 != 'none' and has_cont == '1' : 
     
      # Now retrieve the BP/RP mean spectra in the continuous representation
      result = gaia.load_data(ids=namedr3, format='csv', data_release='Gaia DR3', data_structure='raw', retrieval_type='XP_CONTINUOUS',avoid_datatype_check=True)

      # Result will be a dictionary, so you can check the available keys by running result.keys()
      # In this example we are looking in particular for the XP_CONTINUOUS_RAW key
      #print(result)
      continuous_key = [key for key in result.keys() if 'continuous' in key.lower()][0]
      # The first element is the result we want as an Astropy table
      data = result[continuous_key][0]
      # Astropy has a ’write’ method for tables
      # Write the table to CSV
      #filename=dir+namedr3+'.csv'
      #data.write(filename, format='csv',overwrite=True)

###https://gaia-dpci.github.io/GaiaXPy-website/tutorials/Calibrator%20tutorial.html
      sampling = np.geomspace(330,1049.9999999999, 361)
      data2=data.to_pandas()
      calibrated_spectra, sampling = calibrate(data2, sampling=sampling)

##you can save the plots 
      spectrum_plot(dir,sampling,calibrated_spectra)

#get also light curve
#      data_structure = 'INDIVIDUAL'   # Options are: 'INDIVIDUAL', 'COMBINED', 'RAW'
#      data_release   = 'Gaia DR3'     # Options are: 'Gaia DR3' (default), 'Gaia DR2'
#      #https://www.cosmos.esa.int/web/gaia-users/archive/datalink-products#datalink_jntb_get_all
#      datalink = Gaia.load_data(ids=namedr3, data_release = data_release, retrieval_type=retrieval_type, data_structure = data_structure, verbose = False, output_file = None)
#      dl_keys  = [inp for inp in datalink.keys()]
#      dl_keys.sort()

#      print()
#      print(f'The following Datalink products have been downloaded:')
#      for dl_key in dl_keys:
#         print(f' * {dl_key}')
#         dl_out  = extract_dl_ind(datalink, dl_key, figsize=[20,7])   # Change the figsize to e.g. figsize=[20,7] to increase the size of the displayed image.
#         dl_out[0:5] 
         
         
mycatalog="J/ApJ/822/L5/table1"
print, "*********************************"
print, "get  2MASS IDs from the VIZIER database using coordinates from Vizier cat"
print, "get_twomassid(mycatalog,'tabletwo.csv')"
get_twomassid(mycatalog,'tabletwo.csv') 

print, "*********************************"
print, "get Gaia ID from 2MASS IDs, only stars a 2MASS are used"
print, "get_gaiaid('tabletwo.csv','gaiaID.csv')"
get_gaiaid('tabletwo.csv','gaiaID.csv')


print, "*********IT requires the Gaia accout name and password for authentication, to allow file uploads"
print, "*********enter the module and change 'user_mmessi01.user' "
print, "get Gaia parameters from Gaia IDs"
print, "get_gaia_par('gaiaID.csv','gaia_par.csv')"
get_gaia_par('gaiaID.csv','gaia_par.csv')




