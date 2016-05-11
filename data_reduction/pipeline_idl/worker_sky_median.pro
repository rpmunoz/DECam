pro worker_sky_median, chip, chip_name, im_size, chip_offset, chip_n, im_name, im_file, im_filter, program, out
	; do some calculations here
	im_file=strsplit(im_file,'/', /extract)

	sky_name=im_name[1:-1]
	sky_file=im_file[1:-1]

 	im_data = (shmvar(im_name[0]))[*,*,chip]
 	mask_data = (shmvar('sm_mask_data'))[*,*,chip]
	gv1=where(mask_data[0:1022,*] EQ 0)
	gv2=where(mask_data[1023:-1,*] EQ 0)
	im_median=fltarr(2)
	im_median[0]=median((im_data[0:1022,*])[gv1])
	im_median[1]=median((im_data[1023:-1,*])[gv2])

	sky_median=fltarr([2,n_elements(sky_name)])
	sky_data=fltarr([im_size,n_elements(sky_name)])

	if chip_name EQ 'N15' OR chip_name EQ 'S18' then begin
		cgloadct, 0
		cgwindow, wxsize=1200, wysize=1200
		plot_pos = cgLayout([2,n_elements(sky_name)], OXMargin=[6,4], OYMargin=[2,2], XGap=4, YGap=4)
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
			cgplot, [0], [0], position=plot_pos[*,2*i], xrange=plot_xrange, yrange=median(plot_y)*plot_yrange_factor1, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, xtickformat='(I0)', xtitle=(i EQ n_elements(sky_name)-1 ? 'X (pixels)':''), ytitle='Counts (ADU)', /nodata, /addcmd, /window, /noerase, title='File '+(sky_file[i])[-1], charsize=1.
			cgplot, plot_x, plot_y, position=plot_pos[*,2*i], thick=2, color='black', /over, /addcmd

			plot_x=indgen(sky_size[1])
			plot_y=median(sky_data[*,*,i], dim=1)
			plot_xrange=[0,sky_size[1]]
			cgplot, [0], [0], position=plot_pos[*,2*i+1], xrange=plot_xrange, yrange=median(plot_y)*plot_yrange_factor2, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, xtickformat='(I0)', xtitle=(i EQ n_elements(sky_name)-1 ? 'Y (pixels)':''), /nodata, /addcmd, /window, /noerase, title='File '+(sky_file[i])[-1], charsize=1.
			cgplot, plot_x, plot_y, position=plot_pos[*,2*i+1], thick=2, color='black', /over, /addcmd
		endif

		sky_data[0:1022,*,i] /= sky_median[0,i]
		sky_data[1023:-1,*,i] /= sky_median[1,i]
	endfor
	
	if chip_name EQ 'N15' OR chip_name EQ 'S18' then begin
		cgcontrol, output='results/sky_subtraction_001_'+program+'_'+repstr((im_file[0])[-1],'.fits','')+'_'+chip_name+'.pdf'
		cgdelete, /all
	endif

	out = im_data - rebin([replicate(im_median[0],im_size[0]/2),replicate(im_median[1],im_size[0]/2)], im_size) * median(sky_data,dim=3)
 	im_data=0
 	sky_data=0

end
