import pyfits
import numpy as np
import pylab as py
import img_scale
from astropy import wcs
from astropy.io import fits

im_dir='/Volumes/Q6/NGFS/DECam/stacks/' 
im_file1=im_dir+'fornax_tile1_i.fits'
im_file2=im_dir+'fornax_tile1_g.fits'
im_file3=im_dir+'fornax_tile1_u.fits'
im_coo_wcs = np.array([[55.31831,-36.37102], [53.80851,-34.69870]], np.float_)

r = pyfits.getdata(im_file1)
g = pyfits.getdata(im_file2)
b = pyfits.getdata(im_file3)

im_h = fits.getheader(im_file1)
im_w = wcs.WCS(im_h)
im_coo = np.round(im_w.wcs_world2pix(im_coo_wcs, 1))
im_size=np.array([im_coo[1,0]-im_coo[0,0],im_coo[1,1]-im_coo[0,1]], np.int_)
r = r[im_coo[0,0]:im_coo[0,0]+im_size[0], im_coo[0,1]:im_coo[0,1]+im_size[1]]

im_h = fits.getheader(im_file2)
im_w = wcs.WCS(im_h)
im_coo = np.round(im_w.wcs_world2pix(im_coo_wcs, 1))
g = g[im_coo[0,0]:im_coo[0,0]+im_size[0], im_coo[0,1]:im_coo[0,1]+im_size[1]]

im_h = fits.getheader(im_file3)
im_w = wcs.WCS(im_h)
im_coo = np.round(im_w.wcs_world2pix(im_coo_wcs, 1))
b = b[im_coo[0,0]:im_coo[0,0]+im_size[0], im_coo[0,1]:im_coo[0,1]+im_size[1]]

img = np.zeros((r.shape[0], r.shape[1], 3), dtype=float)
img[:,:,0] = img_scale.sqrt(r, scale_min=-15, scale_max=60)
img[:,:,1] = img_scale.sqrt(g, scale_min=-3, scale_max=15)
img[:,:,2] = img_scale.sqrt(b, scale_min=-2, scale_max=10)

py.clf()
py.imshow(img, aspect='equal')
py.title('Blue = u, Green = g, Red = i')
py.savefig('fornax_tile1_ugi.png')
