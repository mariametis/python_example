import pyvo as vo
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.utils.data import download_file
from astropy.io import fits
import wget

#read table
data=pd.read_csv("../LIST/LBV.v2.csv",dtype={'RAD': 'string', 'DECD': 'string'})
myra=data['RAD']
mydec=data['DECD']
myname=data['Name']
nlen=len(myra)
data=pd.read_csv("../LIST/LBV.v2.csv")
myra_v=data['RAD']
mydec_v=data['DECD']

file=["" for x in range(nlen)]
filefits=["" for x in range(nlen)]
for  k in range(0, nlen):
  file[k]=myname[k]+".eps"
  filefits[k]=myname[k]+".fits"
  if(myname[k] == 'Wra15-751/IRAS11065-6026'):
      file[k]='Wra15-751-IRAS11065-6026'+".eps"
      filefits[k]='Wra15-751-IRAS11065-6026'+".fits"
  if(myname[k] == 'WMD14/IRAS16254-4739'):
      file[k]='WMD14-IRAS16254-4739'+".eps"
      filefits[k]='WMD14-IRAS16254-4739'+".fits"
  if(myname[k] == 'IRAS17050-3921/GKF2010-MN51'):
      file[k]='IRAS17050-3921-GKF2010-MN51'+".eps"
      filefits[k]='IRAS17050-3921-GKF2010-MN51'+".fits"

data['filefits']=filefits
data['file']=file
data.to_csv('LBV.v2.massfits.csv',index=0)


#THIS BLOCK makes the download from skyview
# I can decide where to save the fits
for  k in range(0, nlen):
   ra = myra[k]
   dec = mydec[k]
   ra_v = myra_v[k]
   dec_v = mydec_v[k]
   if ((myname[k] != "none")):
    #url='https://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?mode=getImage&RA='+ra+'&DEC='+dec+'&survey=2mass&band=k&subsetsize=1.0&file_type=fits&reproject=true'
    url="https://skyview.gsfc.nasa.gov/cgi-bin/images?Survey=2mass-k&position='"+ra+","+dec+"'&Coordinates=J2000&Return=FITS&size=0.034"  
    if myname[k]=="GKM2012-WS1":
      url="https://skyview.gsfc.nasa.gov/cgi-bin/images?Survey=2mass-h&position='"+ra+","+dec+"'&Coordinates=J2000&Return=FITS&size=0.034"  
    output_directory = "./myfits/"+filefits[k]
    filefits[k] = wget.download(url, out=output_directory)


#THIS BLOCK makes the encapsulated charts
for  k in range(0, nlen):
   ra = myra[k]
   dec = mydec[k]
   ra_v = myra_v[k]
   dec_v = mydec_v[k]

   pos = SkyCoord(ra=ra_v, dec=dec_v, unit='deg')
   fname= "./myfits/"+filefits[k]
   image1 = fits.open(fname)
   wcs = WCS(image1[0].header)
   cutout = Cutout2D(image1[0].data, pos, (60, 60),  wcs=wcs)
   wcs = cutout.wcs

   fig = plt.figure()
   ax = fig.add_subplot(1, 1, 1, projection=wcs)

   if((myname[k] == 'GCIRS34W') or (myname[k] == 'GCIRC-AF') or (myname[k] == 'GCIRS-34W') or (myname[k] == 'GCIRS-33SE') or (myname[k] == 'GCIRS-16NW') \
    or (myname[k] == 'GCIRC-33E') or (myname[k] == 'GCIRS-16SW') \
    or (myname[k] == 'GCIRS-16C') or (myname[k] =='GCIRS-16NE') ):
      ax.imshow(cutout.data, cmap='gray_r', origin='lower', vmax = 7000)

   if((myname[k] != 'GCIRS34W') and (myname[k] != 'GCIRC-AF') and myname[k] != 'GCIRS-34W' and myname[k] != 'GCIRS-33SE' and myname[k] != 'GCIRS-16NW' \
    and myname[k] != 'GCIRC-33E' and myname[k] != 'GCIRS-16SW' \
    and myname[k] != 'GCIRS-16C' and myname[k] !='GCIRS-16NE'):
      ax.imshow(cutout.data, cmap='gray_r', origin='lower', vmax = 1000)
   
   #ax.imshow(cutout.data, cmap='gray', vmin=0,vmax = 255)
   ax.scatter(ra_v, dec_v, transform=ax.get_transform('fk5'), s=500, edgecolor='red', facecolor='none')
   #plt.show()
   plt.savefig(file[k],dpi=1000)
   plt.close(fig) 

   #sys.exit()
