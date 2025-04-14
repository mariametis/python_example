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


############################
def spectrum_plot(dirout,sampling,calibrated_spectra): 

    tit=str(calibrated_spectra['source_id'].iloc[0])
    file=dirout+str(calibrated_spectra['source_id'].iloc[0])+'.eps'
     
    #print(sampling.dtype) 
    #print(type(sampling)) 
    #print(sampling)

#    sys.exit(0)
 
    wave=sampling.flatten().tolist()  #from nparray to list
    flux=calibrated_spectra['flux'].values #now nparray with ndim=1
    flux=list(flatten(flux)) #from nparray to list
    lala=[wave,flux]
    
    lala=np.asarray(lala)
    columns = ['wave','flux']
    dfsel = pd.DataFrame(lala.T,columns=columns)
    print(dfsel)

    aa=dfsel.loc[dfsel['wave'].between(350,1050 ),['wave']].values
    bb=dfsel.loc[dfsel['wave'].between(350,1050 ),['flux']].values
    wave=aa
    flux=bb

    xywidth=1.1
    linesize=1.1

    xytitlesize=8
    xylenght1=3
    xylenght2=2

    rc('axes', linewidth=xywidth)
    mpl.rcParams['axes.titlesize'] = xytitlesize #?? no effect
    mpl.rcParams['axes.labelsize'] = xytitlesize # this are x- and ytitle
    mpl.rcParams['axes.labelpad'] = 2 # distance from the axis
    mpl.rcParams['figure.labelweight'] = 9 # ????

    mpl.rcParams['xtick.major.size'] = xylenght1
    mpl.rcParams['xtick.major.width'] =xywidth
    mpl.rcParams['xtick.minor.size'] = xylenght2
    mpl.rcParams['xtick.minor.width'] = xywidth
    mpl.rcParams['ytick.major.size'] = xylenght1
    mpl.rcParams['ytick.major.width'] =xywidth
    mpl.rcParams['ytick.minor.size'] = xylenght2
    mpl.rcParams['ytick.minor.width'] = xywidth
    mpl.rcParams['xtick.labelsize'] = xytitlesize  #this are the number labels
    mpl.rcParams['ytick.labelsize'] = xytitlesize

    fontsize = 16
    mpl.rcParams['lines.linewidth'] = linesize
    mpl.rcParams['lines.markersize'] = xywidth
    mpl.rcParams['grid.color'] = '#E6E6E6'   #light gray

    

    rect_2 = [ 0.15, .08, 0.82, 0.21]
    rect_1 = [ 0.15, .08, 0.82, .82]
    fig = plt.figure(figsize=(4, 2))
      
    
    
    box1 = fig.add_subplot(rect_1) #plt.axes(rect_1)
    plt.ylabel('Flux')
    plt.xlabel('wavelenght[nm]')
    plt.title(tit)
    box1.set_xlim([390,1050])
    box1.plot(wave, flux,color ="#000000",)
    box1.tick_params(axis='both', which='major',direction='in', left=True, top=False, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=True)
    box1.tick_params(axis='both', which='minor', direction='in',left=True, top=False, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=True)
#    box1.minorticks_on
#    box1.xaxis.set_major_locator(MultipleLocator(5))
#    box1.xaxis.set_major_formatter('{x:.0f}')
#    box1.xaxis.set_minor_locator(MultipleLocator(1))
#    box1.yaxis.set_minor_locator(AutoMinorLocator())
    #%plt.grid()
    plt.savefig(file,dpi=1000)
   

    
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
         
         


get_bprp('myflag3.csv','./')


