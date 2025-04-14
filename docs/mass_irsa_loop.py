import pyvo as vo
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.utils.data import download_file
from astropy.io import fits


#the problem when crashing is that I cannot find the downloaded fits, download_file utility is unfrendly

#read table
data=pd.read_csv("../LIST/LBV.v2.csv",dtype={'source_id': 'string'})
myra=data['RAD']
mydec=data['DECD']
myname=data['Name']
nlen=len(myra)

#loop the target and save chart
for  k in range(2, nlen):
  ra = myra[k]
  dec = mydec[k]
  file=myname[k]+".eps"
  if ((myname[k] != "W243") & (myname[k] != "etacar")):
   if(myname[k] == 'Wra15-751/IRAS11065-6026'):
      file='Wra15-751-IRAS11065-6026'+".eps"
   pos = SkyCoord(ra=ra, dec=dec, unit='deg')

   twomass_service = vo.dal.SIAService("https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?type=at&ds=asky&")
   im_table = twomass_service.search(pos=pos, size=1.0*u.arcsec)
   im_table.to_table()

   for i in range(len(im_table)):
    if im_table[i]['band'] == 'K':
        break
   print(im_table[i].getdataurl())

   fname = download_file(im_table[i].getdataurl(), cache=True)
   image1 = fits.open(fname)

   wcs = WCS(image1[0].header)
   cutout = Cutout2D(image1[0].data, pos, (60, 60),  wcs=wcs)
   wcs = cutout.wcs

   fig = plt.figure()

   ax = fig.add_subplot(1, 1, 1, projection=wcs)
   ax.imshow(cutout.data, cmap='gray_r', origin='lower', 
          vmax = 1000)
   ax.scatter(ra, dec, transform=ax.get_transform('fk5'), s=500, edgecolor='red', facecolor='none')
   #plt.show()
   plt.savefig(file,dpi=1000)

