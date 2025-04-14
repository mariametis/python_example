import pyvo as vo
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.utils.data import download_file
from astropy.io import fits

ra = 314.30417
dec = 77.595559
pos = SkyCoord(ra=ra, dec=dec, unit='deg')

twomass_service = vo.dal.SIAService("https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?type=at&ds=asky&")
im_table = twomass_service.search(pos=pos, size=1.0*u.arcsec)
im_table.to_table()

for i in range(len(im_table)):
    if im_table[i]['band'] == 'H':
        break
print(im_table[i].getdataurl())

fname = download_file(im_table[i].getdataurl(), cache=True)
image1 = fits.open(fname)

wcs = WCS(image1[0].header)

cutout = Cutout2D(image1[0].data, pos, (60, 60), wcs=wcs)
wcs = cutout.wcs

fig = plt.figure()

ax = fig.add_subplot(1, 1, 1, projection=wcs)
ax.imshow(cutout.data, cmap='gray_r', origin='lower', 
          vmax = 1000)
ax.scatter(ra, dec, transform=ax.get_transform('fk5'), s=500, edgecolor='red', facecolor='none')
plt.show()


