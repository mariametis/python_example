import numpy as np
import matplotlib.pyplot as plt
import astropy.visualization
import sys
import os
import img_scale
from astropy.visualization import make_lupton_rgb
from astropy.wcs import WCS
from astropy.io import fits
from pylab import *

#rc('axes', linewidth=2)
#fontsize = 20

#Python: Spawning another program
#https://stackoverflow.com/questions/34952291/python-spawning-another-program
A='test.fits htan.txt'
cmd='/Users/mmessine/software/Montage_v3.3/bin/mproject '+'IC2391-0074_2001.fits '+A 
print(cmd)
result = os.system(cmd)

hdu1= fits.open('test.fits') #UVES acquisition
hdu1.info()
im1 = hdu1[0].data*1.0

hdu2 = fits.open('dss2-red_IC2391-0074.fits') #dss2
hdu2.info()
im2 = hdu2[0].data*1.0

#https://learn.astropy.org/tutorials/celestial_coords1.html
header = hdu2[0].header
wcs_helix = WCS(header)

#https://www.astrobetter.com/blog/2010/10/22/making-rgb-images-from-fits-files-with-pythonmatplotlib/
#http://astroweb.case.edu/jakub/TA/img_scale.py   img_scale.py
#img = np.zeros((im1.shape[1], im1.shape[0], 3), dtype=float)
#img = np.zeros((60,59, 3), dtype=float)
#img[:,:,0] = im1
#img[:,:,1] = im1
#img[:,:,2] = im1
a = img_scale.range_from_percentile(im1,low_cut=0.1, high_cut=0.65)
b = img_scale.range_from_percentile(im2,low_cut=0.2, high_cut=0.15)
c = img_scale.range_from_percentile(im2,low_cut=0.2, high_cut=0.15)

a = img_scale.linear(im1,scale_min=a[0], scale_max=a[1])
b = img_scale.linear(im2,scale_min=b[0], scale_max=b[1])
c = img_scale.linear(im2,scale_min=c[0], scale_max=c[1])
#https://eteq.github.io/astropy/visualization/lupton_rgb.html
#image_rgb = astropy.visualization.make_lupton_rgb(a,b,c,stretch=1.,filename="qw.jpeg")

image_rgb = astropy.visualization.make_lupton_rgb(a,b,c,filename="qw.jpeg")

fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(projection=wcs_helix) #;;projection here and get_coord_overlay are different
plt.imshow(image_rgb, origin='lower')
#plt.colorbar()
plt.xlabel(r'RA')
plt.ylabel(r'Dec')
overlay = ax.get_coords_overlay('icrs')
#overlay.grid(color='white', ls='dotted')
#plt.show()
plt.savefig('foo.png')

sys.exit()



#reproject
#https://astronomy.stackexchange.com/questions/51766/making-rgb-image-using-fits-files-with-different-dimensions
#ndata, _ = reproject.reproject_interp(hdu2[0], hdu1[0].header)
#datat, _ = reproject.reproject_interp(hdu3[0], hdu1[0].header)


plt.imshow(im1, cmap='gray')
plt.colorbar()
plt.show()
sys.exit()

