pro worker_sky_surface, chip, chip_name, im_size, chip_offset, chip_n, im_file, im_name, sky_file, sky_name, mask_file, mask_name, im_filter, program, sky_region, out_im, out_sky, plot_pos_chip

	method='tps'
	im_file=strsplit(im_file,'/', /extract)
	sky_file=strsplit(sky_file,'/', /extract)

 	im_data = (shmvar(im_name))[*,*,chip]
 	mask_data = (shmvar(mask_name))[*,*,chip]
	im_size=size(im_data, /dim)

	out_im_masked = im_data
	bv=where(mask_data EQ 1, n_bv, COMPLEMENT=gv, NCOMPLEMENT=n_gv)
	gv=gv[sort(out_im_masked[gv])]
	bv_low=[] ;gv[0:n_gv*0.02] ; & out_im_masked[bv_low] = !VALUES.F_NAN
	bv_high=gv[n_gv*0.99:-1] ; & out_im_masked[bv_high] = !VALUES.F_NAN
	bv=[bv,bv_low,bv_high]
	out_im_masked[bv] = !VALUES.F_NAN

	out_sky=make_array(im_size, value=0., /FLOAT)
	out_sky_x=findgen(im_size[0])#make_array(im_size[1], value=1., /FLOAT)
	out_sky_y=make_array(im_size[0], value=1., /FLOAT)#findgen(im_size[1])
	
	sky_region_size=size(sky_region, /dim)

	for i=0L, sky_region_size[1]-1 do begin

		nsky_region=sky_region[*,i]
		nsky_data=out_im_masked[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]]
		nsky_size=size(nsky_data, /dim)

		if method EQ 'tps' then begin

			median_width=128L

    	nim_data=fltarr( ceil(nsky_size[0]*1./median_width)*median_width, ceil(nsky_size[1]*1./median_width)*median_width)
 	  	nim_data[*]=!VALUES.F_NAN
  	  nim_data[0:nsky_size[0]-1,0:nsky_size[1]-1]=nsky_data
    	nim_size=size(nim_data,/dim)

   		a=reform(rebin(reform(indgen(nim_size[0]),[median_width,1,nim_size[0]/median_width]),[median_width,median_width,nim_size[0]/median_width]), [median_width*median_width,nim_size[0]/median_width,1])
   		a=reform(rebin(a,[median_width*median_width,nim_size[0]/median_width,nim_size[1]/median_width]),[median_width*median_width,nim_size[0]*nim_size[1]/(median_width*median_width)])

			b=reform(rebin(transpose(indgen(nim_size[1])), [median_width,nim_size[1]]),[median_width*median_width,1,nim_size[1]/median_width])
  	  b=reform(rebin(b,[median_width*median_width,nim_size[0]/median_width,nim_size[1]/median_width]),[median_width*median_width,nim_size[0]*nim_size[1]/(median_width*median_width)])

			nim_fraction=reform( total( finite(nim_data[nim_size[0]*b+a]), 1)/median_width^2, nim_size/median_width)
			nim_data=reform( median(nim_data[nim_size[0]*b+a],dim=1), nim_size/median_width)
  	  nim_size=size(nim_data,/dim)

			nim_x=( median_width*findgen(nim_size[0]) + (median_width-1.)/2)#make_array(nim_size[1], value=1., /float)
			nim_y=make_array(nim_size[0], value=1., /float)#( median_width*findgen(nim_size[1]) + (median_width-1.)/2)

;			gv=where(finite(nim_data), n_gv)
			gv=where(nim_fraction GT 0.5, n_gv)
			out_sky[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]]=grid_tps(nim_x[gv], nim_y[gv], nim_data[gv], NGRID=nsky_size, START=[0,0], DELTA=1.)

		endif else $
		if method EQ 'sfit' then begin
	
			temp_x=out_sky_x[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]]
			temp_y=out_sky_y[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]]
			gv=where(finite(nsky_data), n_gv)
			gv_ind=array_indices(nsky_size, gv, /dim)
		
			sfit_result = sfit( [gv_ind,transpose(nsky_data[gv])], 2, /IRREGULAR, /MAX_DEGREE, kx=sfit_coeff)
			out_sky[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]] = sfit_coeff[0] + sfit_coeff[1]*temp_y + sfit_coeff[2]*temp_y^2 + sfit_coeff[3]*temp_x + sfit_coeff[4]*temp_x*temp_y + sfit_coeff[5]*temp_x^2
		endif

	endfor
	
	out_im = im_data - out_sky

 	im_data=0
	nsky_data=0	

end
