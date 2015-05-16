#! /usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

import sys,os
import os.path
import subprocess
import numpy as np
import pyfits
import multiprocessing, Queue
import ctypes
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from astropy.stats import sigma_clipped_stats
from photutils import detect_sources, segment_properties, properties_table

sigmatofwhm=2*np.sqrt(2*np.log(2))
fwhmtosigma=1./sigmatofwhm

class Worker_convolve(multiprocessing.Process):

	def __init__(self, work_queue, result_queue):

		# base class initialization
		multiprocessing.Process.__init__(self)

		# job management stuff
		self.work_queue = work_queue
		self.result_queue = result_queue
		self.kill_received = False

	def run(self):
		while not self.kill_received:

			# get a task
			try:
				x_range, y_range = self.work_queue.get_nowait()
			except Queue.Empty:
				break

			# the actual processing
#			print 'Worker convolve started ', id
			print "Convolving image section - XRANGE=", y_range, " - YRANGE=", x_range

			im_size=np.asarray(shared_im.shape)
			kernel_size=np.asarray(shared_kernel.shape)

			im_pad=(kernel_size-1)/2+10
			x_pad=[0 if x_range[0]==0 else -im_pad[0], 0 if x_range[1]==(im_size[0]-1) else im_pad[0]]
			y_pad=[0 if y_range[0]==0 else -im_pad[1], 0 if y_range[1]==(im_size[1]-1) else im_pad[1]]


			shared_nim[x_range[0]:x_range[1], y_range[0]:y_range[1]]=convolve(shared_im[x_range[0]+x_pad[0]:x_range[1]+x_pad[1], y_range[0]+y_pad[0]:y_range[1]+y_pad[1]], shared_kernel, normalize_kernel=True)[-x_pad[0]:x_range[1]-x_range[0]-x_pad[0], -y_pad[0]:y_range[1]-y_range[0]-y_pad[0]]

			print 'Worker convolve is done', id
			self.result_queue.put(id)


if __name__ == "__main__":

	im_fix_file='/Volumes/Q6/NGFS/DECam/stacks/ss_fornax_tile1_g_long_ALIGNi_FIX.003.fits'
	weight_fix_file='/Volumes/Q6/NGFS/DECam/stacks/ss_fornax_tile1_g_long_ALIGNi_FIX.003.WEIGHT.fits'
	check_fix_file='/Volumes/Q6/NGFS/DECam/stacks/check/ss_fornax_FIX.003.CHECK_SEGMENTATION.fits'
	cat_fix_file='/Volumes/Q6/NGFS/DECam/stacks/check/ss_fornax_FIX.003.ldac'
	
	im_file='/Volumes/Q6/NGFS/DECam/stacks/ss_fornax_tile1_g_long_ALIGNi.003.fits'
	im_convol_file='/Volumes/Q6/NGFS/DECam/stacks/ss_fornax_tile1_g_long_ALIGNi_CONVOL.003.fits'
	im_convol_seg_file='/Volumes/Q6/NGFS/DECam/stacks/ss_fornax_tile1_g_long_ALIGNi_CONVOL_SEG.003.fits'

	pix_scale=0.2637	# arcsec/pixel
	kernel_sigma=20.
	kernel_size=[101,101]
	kernel_small_sigma= 2./pix_scale * fwhmtosigma
	kernel_small_size= np.full(2, round(4*kernel_small_sigma/2.)*2+1, dtype=np.int)
	kernel_large_sigma= 40./pix_scale * fwhmtosigma
	kernel_large_size= np.full(2, round(4*kernel_large_sigma/2.)*2+1, dtype=np.int)
	grid_n=[6,6]
	n_cpu=2
	n_core=4

	n_processes=n_cpu*n_core*1
	
	hdulist = pyfits.open(im_fix_file)
	im_data=hdulist[0].data
	im_h=hdulist[0].header
	hdulist.close()
	im_size=np.asarray(im_data.shape)
	
	hdulist = pyfits.open(weight_fix_file)
	weight_data=hdulist[0].data
	weight_h=hdulist[0].header
	hdulist.close()

	im_data[weight_data == 0.]=np.nan
	im_mask_nan = np.isnan(im_data)

	shared_im_base = multiprocessing.Array(ctypes.c_float, im_size[0]*im_size[1])
	shared_im = np.ctypeslib.as_array(shared_im_base.get_obj())
	shared_im = shared_im.reshape(im_size[0], im_size[1])
	shared_im[:] = im_data
	im_data=0

	if os.path.isfile(im_convol_seg_file) is False:

		shared_nim_base = multiprocessing.Array(ctypes.c_float, im_size[0]*im_size[1])
		shared_nim = np.ctypeslib.as_array(shared_nim_base.get_obj())
		shared_nim = shared_nim.reshape(im_size[0], im_size[1])

		print "Generating kernel to mask small objects"
		print 'Kernel_size: ', kernel_small_size, ' - Kernel_sigma: ', kernel_small_sigma
		kernel_data=np.asarray(Gaussian2DKernel(kernel_small_sigma, x_size=kernel_small_size[0], y_size=kernel_small_size[1], mode='integrate'))
		kernel_size=kernel_data.shape

		shared_kernel_base = multiprocessing.Array(ctypes.c_float, kernel_size[0]*kernel_size[1])
		shared_kernel = np.ctypeslib.as_array(shared_kernel_base.get_obj())
		shared_kernel = shared_kernel.reshape(kernel_size[0], kernel_size[1])
		shared_kernel[:] = kernel_data

		work_queue = multiprocessing.Queue()
		grid_n=np.asarray(grid_n)
		grid_mesh=np.ceil(im_size*1./grid_n).astype(int)
		grid_x=np.append( np.arange(0,im_size[0], grid_mesh[0]), im_size[0])
		grid_y=np.append( np.arange(0,im_size[1], grid_mesh[1]), im_size[1])
	
		for i in range(0,grid_x.size-1):
			for j in range(0,grid_y.size-1):
				if work_queue.full():
					print "Oh no! Queue is full after only %d iterations" % j
	
				x_range=[grid_x[i],grid_x[i+1]]
				y_range=[grid_y[j],grid_y[j+1]]
				work_queue.put( (x_range, y_range) )
	
		# create a queue to pass to workers to store the results
		result_queue = multiprocessing.Queue()
		procs=[]
	
		# spawn workers
		for i in range(n_processes):
			worker = Worker_convolve(work_queue, result_queue)
			procs.append(worker)
			worker.start()
	
		# collect the results off the queue
		for i in range(n_processes):
			result_queue.get()
	
		for p in procs:
			p.join()
	
		print "Computing median and standard deviation"
	
		x_region=np.round( im_size[0]*1./3+[-500,500] )
		y_region=np.round( im_size[1]*1./3+[-500,500] )
		print 'NaN region1 ', np.where( np.isnan(shared_im[x_region[0]:x_region[1], y_region[0]:y_region[1]]) )
		mean1, median1, std1 = sigma_clipped_stats(shared_im[x_region[0]:x_region[1], y_region[0]:y_region[1]], mask=im_mask_nan[x_region[0]:x_region[1], y_region[0]:y_region[1]], sigma=3.0)	
	
		x_region=np.round( im_size[0]*2./3+[-500,500] )
		y_region=np.round( im_size[1]*2./3+[-500,500] )
		print 'NaN region2 ', np.where( np.isnan(shared_im[x_region[0]:x_region[1], y_region[0]:y_region[1]]) )
		mean2, median2, std2 = sigma_clipped_stats(shared_im[x_region[0]:x_region[1], y_region[0]:y_region[1]], mask=im_mask_nan[x_region[0]:x_region[1], y_region[0]:y_region[1]] , sigma=3.0)	
		
		print "Statistics region 1 ", median1, std1
		print "Statistics region 2 ", median2, std2
	
		im_median=np.mean([median1,median2])
		im_stddev=np.sqrt((std1**2+std2**2)/2)
		im_thresh=im_median + 2*im_stddev
		print "Image median, stddev, threshold ", im_median, im_stddev, im_thresh
	
		seg_data = detect_sources(shared_nim, im_thresh, npixels=5)
	
		print "Writing file ", im_convol_file
		
		if os.path.isfile(im_convol_file): os.remove(im_convol_file)
		pyfits.writeto(im_convol_file, shared_nim, header=im_h)
	
		if os.path.isfile(im_convol_seg_file): os.remove(im_convol_seg_file)
		pyfits.writeto(im_convol_seg_file, seg_data, header=im_h)
	else:
		hdulist = pyfits.open(im_convol_seg_file)
		seg_data=hdulist[0].data
		seg_h=hdulist[0].header
		hdulist.close()


# Here we plot the sources detected using the segmentation map
	seg_nid=np.max(seg_data)+1

	cat_data=np.recarray(seg_nid, dtype={'names':['id','npix','flux','mag','fwhm'], 'formats':['i4','i4','f4','f4','f4']})
	cat_data.fill(0)

	print 'Computing the flux and size for the detected sources'
	seg_npix=np.bincount(np.ravel(seg_data))
	seg_flux=np.bincount(np.ravel(seg_data), weights=np.ravel(shared_im))
	for i in range(seg_nid):
		cat_data[i]=(i, seg_npix[i], seg_flux[i], 0., 0.)

	cat_data=cat_data[1:]

	fig = plt.figure()
	ax = plt.gca()
	ax.plot(np.sqrt(cat_data['npix']), cat_data['flux'], linestyle='', marker='o', markersize=5, c='blue', alpha=0.5, markeredgecolor='none')
	ax.set_yscale('log')
	ax.set_xlabel('Sqrt(Area) (pix)')
	ax.set_ylabel('Flux (ADU)')

	fig.savefig('Mag_size_plot_tile1_g.pdf', format='pdf')

	source_props = segment_properties(shared_im, seg_data)
	source_table = properties_table(source_props)
	print source_table

	sys.exit()
	
	hdulist = pyfits.open(check_fix_file)
	check_data=hdulist[0].data
	check_h=hdulist[0].header
	hdulist.close()
	
	hdulist = pyfits.open(cat_fix_file)
	cat_data=hdulist[2].data
	hdulist.close()
	print "Sextractor catalog field names", cat_data.columns.names
	
	
	gv_cat=~( (cat_data['FLUX_RADIUS']>10.) & (cat_data['MAG_AUTO']>19.) )
	gv_mask=np.in1d( check_data, cat_data['NUMBER'][gv_cat]).reshape( check_data.shape )
	im_data[gv_mask]=np.nan
	
	print "Writing file ", im_convol_file
	
	if os.path.isfile(im_convol_file): os.remove(im_convol_file)
	pyfits.writeto(im_convol_file, shared_nim, header=im_h)

