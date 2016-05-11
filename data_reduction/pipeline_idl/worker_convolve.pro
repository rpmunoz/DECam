pro worker_convolve, tile, im_size, xrange, yrange, out
	; do some calculations here
 	in = (shmvar('sm_im_data_tile'+tile))[xrange[0]:xrange[1],yrange[0]:yrange[1]]

  out = (convolve(in, shmvar('gauss_kernel_small'), FT_IMAGE=in_fft, /no_pad) GT 0.3)
  out = temporary(out) + (convolve(in, shmvar('gauss_kernel_medium'), FT_IMAGE=in_fft, /no_pad) GT 0.3)
 	out = temporary(out) + dilate(in,[[0,1,0],[1,1,1],[0,0,0]])
 	out = temporary(out) + byte(in)
 	in = 0

 	out = out[ (xrange[0] EQ 0L ? 0L : 200L):((xrange[1] EQ im_size[0]-1) ? -1L : -201L), (yrange[0] EQ 0L ? 0L : 200L):((yrange[1] EQ im_size[1]-1) ? -1L : -201L) ]

end
