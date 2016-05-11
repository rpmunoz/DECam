pro worker_sky_median_global, chip, chip_name, im_size, chip_offset, chip_n, im_file, im_name, sky_file, sky_name, mask_file, mask_name, im_filter, program, out_im, out_sky
	; do some calculations here
	im_file=strsplit(im_file,'/', /extract)
	sky_file=strsplit(sky_file,'/', /extract)

 	im_data = (shmvar(im_name))[*,*,chip]
 	mask_data = (shmvar(mask_name))[*,*,chip]
	gv1=where(mask_data[0:1022,*] EQ 0)
	gv2=where(mask_data[1023:-1,*] EQ 0)
	im_median=fltarr(2)
	im_median[0]=median((im_data[0:1022,*])[gv1])
	im_median[1]=median((im_data[1023:-1,*])[gv2])

	sky_median=fltarr([2,n_elements(sky_name)])
	sky_data=fltarr([im_size,n_elements(sky_name)])

	if chip_name EQ 'N15' OR chip_name EQ 'S18' then begin
		cgloadct, 0
		if n_elements(sky_name) LE 5 then begin
			cgwindow, wxsize=1200, wysize=1200
			plot_pos = cgLayout([2,n_elements(sky_name)], OXMargin=[6,4], OYMargin=[2,2], XGap=4, YGap=4)
		endif else begin
			cgwindow, wxsize=1200, wysize=1600
			plot_pos = cgLayout([2,n_elements(sky_name)], OXMargin=[6,4], OYMargin=[2,2], XGap=4, YGap=2)
		endelse
	endif

	plot_yrange_delta1=list()
	plot_yrange_factor1=list()
	plot_yrange_delta2=list()
	plot_yrange_factor2=list()
	gv=(sort(randomu(seed, n_elements(sky_name))))[0:1]
	foreach i, gv do begin
		temp_data=(shmvar(sky_name[i]))[*,*,chip]

		plot_x=indgen(im_size[0])
		plot_y=median(temp_data, dim=2)
		coeff=robust_poly_fit(plot_x, plot_y, 3, temp_yfit, temp_sigma)
		plot_yrange_delta1.add, (max(temp_yfit)-min(temp_yfit))
		plot_yrange_factor1.add, (median(plot_y) + (plot_yrange_delta1[-1])*[-1,1])/median(plot_y)

		plot_x=indgen(im_size[1])
		plot_y=median(temp_data, dim=1)
		coeff=robust_poly_fit(plot_x, plot_y, 3, temp_yfit, temp_sigma)
		plot_yrange_delta2.add, (max(temp_yfit)-min(temp_yfit))
		plot_yrange_factor2.add, (median(plot_y) + (plot_yrange_delta2[-1])*[-1,1])/median(plot_y)
	endforeach

	temp_max=max(plot_yrange_delta1.toarray(type='float'), gv_max1)
	temp_max=max(plot_yrange_delta2.toarray(type='float'), gv_max2)
	plot_yrange_factor1=plot_yrange_factor1[gv_max1]
	plot_yrange_factor2=plot_yrange_factor2[gv_max2]

	for i=0L, n_elements(sky_name)-1 do begin
	 	sky_data[*,*,i] = (shmvar(sky_name[i]))[*,*,chip]
		sky_median[0,i]=median(sky_data[0:1022,*,i])
		sky_median[1,i]=median(sky_data[1023:-1,*,i])
		sky_size=size(sky_data[*,*,i], /dim)

		if chip_name EQ 'N15' OR chip_name EQ 'S18' then begin
			plot_x=indgen(sky_size[0])
			plot_y=median(sky_data[*,*,i], dim=2)
			plot_xrange=[0,sky_size[0]]
			cgplot, [0], [0], position=plot_pos[*,2*i], xrange=plot_xrange, yrange=median(plot_y)*plot_yrange_factor1, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, xtickformat='(I0)', xtitle=(i EQ n_elements(sky_name)-1 ? 'X (pixels)':''), ytitle='Counts (ADU)', /nodata, /addcmd, /window, /noerase, title='File '+(sky_file[i])[-1], charsize=0.6
			cgplot, plot_x, plot_y, position=plot_pos[*,2*i], thick=2, color='black', /over, /addcmd

			plot_x=indgen(sky_size[1])
			plot_y=median(sky_data[*,*,i], dim=1)
			plot_xrange=[0,sky_size[1]]
			cgplot, [0], [0], position=plot_pos[*,2*i+1], xrange=plot_xrange, yrange=median(plot_y)*plot_yrange_factor2, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, xtickformat='(I0)', xtitle=(i EQ n_elements(sky_name)-1 ? 'Y (pixels)':''), /nodata, /addcmd, /window, /noerase, title='File '+(sky_file[i])[-1], charsize=0.6
			cgplot, plot_x, plot_y, position=plot_pos[*,2*i+1], thick=2, color='black', /over, /addcmd
		endif

		sky_data[0:1022,*,i] /= sky_median[0,i]
		sky_data[1023:-1,*,i] /= sky_median[1,i]
	endfor

	if chip_name EQ 'N15' OR chip_name EQ 'S18' then begin
		cgcontrol, output='results/sky_subtraction_002_'+program+'_'+repstr(im_file[-1],'.fits','')+'_'+chip_name+'.pdf'
		cgdelete, /all
	endif
	
	sky_median_data=median(sky_data,dim=3)
	for i=0L, n_elements(sky_name)-1 do begin
	 	sky_data[*,*,i] = (shmvar(sky_name[i]))[*,*,chip]
		sky_data[0:1022,*,i] /= median(sky_data[0:1022,*,i]/sky_median_data[0:1022,*])
		sky_data[1023:-1,*,i] /= median(sky_data[1023:-1,*,i]/sky_median_data[1023:-1,*])
	endfor
	sky_median_data=0

	if im_filter EQ 'u' then begin
		width=15

    out_sky=fltarr( ceil(im_size[0]*1./width)*width, ceil(im_size[1]*1./width)*width)
    out_sky[*]=!VALUES.F_NAN
    out_sky[0:im_size[0]-1,0:im_size[1]-1]=median(sky_data,dim=3)
		sky_data=0
    nim_size=size(out_sky,/dim)

    a=reform(rebin(reform(indgen(nim_size[0]),[width,1,nim_size[0]/width]),[width,width,nim_size[0]/width]), [width*width,nim_size[0]/width,1])
    a=reform(rebin(a,[width*width,nim_size[0]/width,nim_size[1]/width]),[width*width,nim_size[0]*nim_size[1]/(width*width)])

    b=reform(rebin(transpose(indgen(nim_size[1])), [width,nim_size[1]]),[width*width,1,nim_size[1]/width])
    b=reform(rebin(b,[width*width,nim_size[0]/width,nim_size[1]/width]),[width*width,nim_size[0]*nim_size[1]/(width*width)])

    out_sky=median(out_sky[nim_size[0]*b+a],dim=1)
		gv_ind = total( reform( array_indices( nim_size, nim_size[0]*b+a, /dim), [2,size(a,/dim)] ), 2)/width^2
 		nim_x=transpose(gv_ind[0,*])
		nim_y=transpose(gv_ind[1,*])
		a=0 & b=0 & gv_ind=0

		gv=where(finite(out_sky), n_gv, COMPLEMENT=bv, NCOMPLEMENT=n_bv)
		triangulate, nim_x[gv], nim_y[gv], im_triangles, im_bounds
;		nim_data = griddata( nim_x[gv], nim_y[gv], nim_data[gv], /NATURAL, TRIANGLES=im_triangles, /GRID, XOUT=indgen(im_size[0]), YOUT=indgen(im_size[1]) )
;		out_sky = griddata( nim_x[gv], nim_y[gv], nim_data[gv], /KRIGING, TRIANGLES=im_triangles, /GRID, XOUT=indgen(im_size[0]), YOUT=indgen(im_size[1]), ANISOTROPY=[1,1,0], SEARCH=[512], MAX=12, SECTORS=2, VARIOGRAM=[2,100.,1.,2.], MIN_POINTS=2 )
		param=[500.,8e-5,4e-5]
		out_sky = cgKrig2D(out_sky[gv], nim_x[gv], nim_y[gv], SPHERICAL=param, XOUT=indgen(im_size[0]), YOUT=indgen(im_size[1]))
		
;		bv=where(finite(out_sky) EQ 0, n_bv)
;		if n_bv GT 0 then out_sky[bv]=nim_data[bv]
	endif $
	else begin
		out_sky=median(sky_data,dim=3)
		sky_data=0
	endelse

	gv=where(mask_data NE 0, n_gv)
	out_im=im_data
	out_im[gv]=!VALUES.F_NAN
	out_im /= out_sky
	out_median=[ median(out_im[0:1022,*]),median(out_im[1023:-1,*]) ]

	out_im = im_data - rebin([replicate(out_median[0],im_size[0]/2),replicate(out_median[1],im_size[0]/2)], im_size) * out_sky

 	im_data=0

end
