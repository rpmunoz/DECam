pro worker_sky_surface, chip, chip_name, im_size, chip_offset, chip_n, im_file, im_name, sky_file, sky_name, mask_file, mask_name, im_filter, program, sky_region, out_im, out_im_masked, out_im_masked_median, out_im_masked_res, out_sky, plot_pos_chip, plot_x, plot_y, plot_xrange, plot_yrange

	method='tps'
	im_file=strsplit(im_file,'/', /extract)
	sky_file=strsplit(sky_file,'/', /extract)

 	im_data = (shmvar(im_name))[*,*,chip]
 	mask_data = (shmvar(mask_name))[*,*,chip]
	im_size=size(im_data, /dim)

	out_im_masked = im_data
	bv=where(mask_data EQ 1, n_bv, COMPLEMENT=gv, NCOMPLEMENT=n_gv)
;	if n_bv GT 0 then out_im_masked[bv] = !VALUES.F_NAN
	gv=gv[sort(out_im_masked[gv])]
	bv_low=gv[0:n_gv*0.05] ; & out_im_masked[bv_low] = !VALUES.F_NAN
	bv_high=gv[n_gv*0.95:-1] ; & out_im_masked[bv_high] = !VALUES.F_NAN
	bv=[bv,bv_low,bv_high]
	out_im_masked[bv] = !VALUES.F_NAN

;	out_im_masked_median=median(out_im_masked, 5)

	out_sky=make_array(im_size, value=0., /FLOAT)
	out_sky_x=findgen(im_size[0])#make_array(im_size[1], value=1., /FLOAT)
	out_sky_y=make_array(im_size[0], value=1., /FLOAT)#findgen(im_size[1])
	
	sky_region_size=size(sky_region, /dim)
;	plot_x={p1:[1,2,3,4]}
;	plot_y={p1:[5,6,7,8]}
;	plot_xrange={p1:[1,10]}
;	plot_yrange={p1:[1,100]}
;	plot_pos={p1:[0,0,1,1]}

	for i=0L, sky_region_size[1]-1 do begin

		nsky_region=sky_region[*,i]
		nsky_data=out_im_masked[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]]
		nsky_size=size(nsky_data, /dim)

		if method EQ 'tps' then begin

			median_width=50

    	nim_data=fltarr( ceil(nsky_size[0]*1./median_width)*median_width, ceil(nsky_size[1]*1./median_width)*median_width)
 	  	nim_data[*]=!VALUES.F_NAN
  	  nim_data[0:nsky_size[0]-1,0:nsky_size[1]-1]=nsky_data
    	nim_size=size(nim_data,/dim)

   		a=reform(rebin(reform(indgen(nim_size[0]),[median_width,1,nim_size[0]/median_width]),[median_width,median_width,nim_size[0]/median_width]), [median_width*median_width,nim_size[0]/median_width,1])
   		a=reform(rebin(a,[median_width*median_width,nim_size[0]/median_width,nim_size[1]/median_width]),[median_width*median_width,nim_size[0]*nim_size[1]/(median_width*median_width)])

			b=reform(rebin(transpose(indgen(nim_size[1])), [median_width,nim_size[1]]),[median_width*median_width,1,nim_size[1]/median_width])
  	  b=reform(rebin(b,[median_width*median_width,nim_size[0]/median_width,nim_size[1]/median_width]),[median_width*median_width,nim_size[0]*nim_size[1]/(median_width*median_width)])

			nim_data=reform( median(nim_data[nim_size[0]*b+a],dim=1), nim_size/median_width)
  	  nim_size=size(nim_data,/dim)

			nim_x=( median_width*findgen(nim_size[0]) + (median_width-1.)/2)#make_array(nim_size[1], value=1., /float)
			nim_y=make_array(nim_size[0], value=1., /float)#( median_width*findgen(nim_size[1]) + (median_width-1.)/2)

			gv=where(finite(nim_data), n_gv)
			out_sky[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]]=grid_tps(nim_x[gv], nim_y[gv], nim_data[gv], NGRID=nsky_size, START=[0,0], DELTA=1.)

		endif else $
		if method EQ 'sfit' then begin
	
			temp_x=out_sky_x[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]]
			temp_y=out_sky_y[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]]
			gv=where(finite(nsky_data), n_gv)
			gv_ind=array_indices(nsky_size, gv, /dim)
		
			sfit_result = sfit( [gv_ind,transpose(nsky_data[gv])], 2, /IRREGULAR, /MAX_DEGREE, kx=sfit_coeff)
			out_sky[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]] = sfit_coeff[0] + sfit_coeff[1]*temp_y + sfit_coeff[2]*temp_y^2 + sfit_coeff[3]*temp_x + sfit_coeff[4]*temp_x*temp_y + sfit_coeff[5]*temp_x^2
;;		sfit_result = sfit( [gv_ind,transpose(nsky_data[gv])], 3, /IRREGULAR, /MAX_DEGREE, kx=sfit_coeff)
;;		out_sky[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]] = sfit_coeff[0] + sfit_coeff[1]*temp_y + sfit_coeff[2]*temp_y^2 + sfit_coeff[3]*temp_y^3 + sfit_coeff[4]*temp_x + sfit_coeff[5]*temp_x*temp_y + sfit_coeff[6]*temp_x*temp_y^2 + sfit_coeff[7]*temp_x^2 + sfit_coeff[8]*temp_x^2*temp_y + sfit_coeff[9]*temp_x^3
;;		sfit_result = sfit( [gv_ind,transpose(nsky_data[gv])], 2, /IRREGULAR, kx=sfit_coeff)
;;		out_sky[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]] = sfit_coeff[0] + sfit_coeff[1]*temp_y + sfit_coeff[2]*temp_y^2 + sfit_coeff[3]*temp_x + sfit_coeff[4]*temp_x*temp_y + sfit_coeff[5]*temp_x*temp_y^2 + sfit_coeff[6]*temp_x^2 + sfit_coeff[7]*temp_x^2*temp_y + sfit_coeff[8]*temp_x^2*temp_y^2
		endif

		; First, plot against X axis	
;		plot_x.add, nsky_region[0]+indgen(nsky_size[0])
;		plot_y.add, median(nsky_data, dim=2)
;		plot_xrange.add, [min(plot_x[-1]), max(plot_x[-1])] ;  (sky_region[i])[0]+[0,nsky_size[0]]
;		plot_yrange.add, median(plot_y[-1])*[0.98,1.02]
;		plot_pos.add, plot_pos_chip[*,2*i+0]

		; Second, plot against Y axis	
;    plot_x.add, nsky_region[1]+indgen(nsky_size[1])
;    plot_y.add, median(nsky_data, dim=1)
;    plot_xrange.add, [min(plot_x[-1]), max(plot_x[-1])]
;    plot_yrange.add, median(plot_y[-1])*[0.98,1.02]
;		plot_pos.add, plot_pos_chip[*,2*i+1]
	
;		if chip_name EQ 'S4' OR chip_name EQ 'N4' OR chip_name EQ 'N15' OR chip_name EQ 'S18' then begin

			; First, plot against X axis	
;			case chip_name of
;				'S4': nplot_pos=plot_pos[*,0]
;				'N4': nplot_pos=plot_pos[*,2]
;				'S18': nplot_pos=plot_pos[*,4]
;				'N15': nplot_pos=plot_pos[*,6]
;			endcase
		
;			delta=(nplot_pos[3]-nplot_pos[1])/sky_region_size[1]
;			nplot_pos[1] = nplot_pos[1] + i*delta
;			nplot_pos[3] = nplot_pos[1] + delta

;			plot_x=nsky_region[0]+indgen(nsky_size[0])
;			plot_y=median(nsky_data, dim=2)

;			plot_xrange=[min(plot_x), max(plot_x)] ;  (sky_region[i])[0]+[0,nsky_size[0]]
;			plot_yrange=median(plot_y)*[0.98,1.02]
;;;			cgplot, [0], [0], position=nplot_pos, xrange=plot_xrange, yrange=plot_yrange, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, xtickformat='(I0)', xtitle=(i EQ 0 ? 'X (pixels)':''), ytitle='Counts (ADU)', /nodata, /addcmd, /window, /noerase, title='CHIP '+chip_name, charsize=0.6
;;;			cgplot, plot_x, plot_y, position=nplot_pos, thick=2, color='black', /over, /addcmd
;;			plot, [0], [0], position=nplot_pos, xrange=plot_xrange, yrange=plot_yrange, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, xtickformat='(I0)', xtitle=(i EQ 0 ? 'X (pixels)':''), ytitle='Counts (ADU)', /nodata, /noerase, title='CHIP '+chip_name, charsize=0.6, xtickname=replicate(' ', 30)
;;			oplot, plot_x, plot_y, thick=2
;			p=plot(plot_x, plot_y, position=nplot_pos, xrange=plot_xrange, yrange=plot_yrange, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, /current)

			; Second, plot against Y axis
;			case chip_name of
;				'S4': nplot_pos=plot_pos[*,1]
;				'N4': nplot_pos=plot_pos[*,3]
;				'S18': nplot_pos=plot_pos[*,5]
;				'N15': nplot_pos=plot_pos[*,7]
;			endcase
;      delta=(nplot_pos[3]-nplot_pos[1])/sky_region_size[1]
;      nplot_pos[1] = nplot_pos[1] + i*delta
;      nplot_pos[3] = nplot_pos[1] + delta
;      nsky_data=out_im_masked[nsky_region[0]:nsky_region[2] , nsky_region[1]:nsky_region[3]]
;      nsky_size=size(nsky_data, /dim)

;      plot_x=nsky_region[1]+indgen(nsky_size[1])
;      plot_y=median(nsky_data, dim=1)
      
;      plot_xrange=[min(plot_x), max(plot_x)]
;      plot_yrange=median(plot_y)*[0.98,1.02]
;;;      cgplot, [0], [0], position=nplot_pos, xrange=plot_xrange, yrange=plot_yrange, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, xtickformat='(I0)', xtitle=(i EQ 0 ? 'Y (pixels)':''), ytitle='Counts (ADU)', /nodata, /addcmd, /window, /noerase, title='CHIP '+chip_name, charsize=0.6
;;;      cgplot, plot_x, plot_y, position=nplot_pos, thick=2, color='black', /over, /addcmd
;;      plot, [0], [0], position=nplot_pos, xrange=plot_xrange, yrange=plot_yrange, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, xtickformat='(I0)', xtitle=(i EQ 0 ? 'Y (pixels)':''), ytitle='Counts (ADU)', /nodata, /noerase, title='CHIP '+chip_name, charsize=0.6, xtickname=replicate(' ', 30)
;;			oplot, plot_x, plot_y, thick=2
;			p=plot(plot_x, plot_y, position=nplot_pos, xrange=plot_xrange, yrange=plot_yrange, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, /current)

;		endif

	endfor
	
	out_im = im_data - out_sky
;	out_im_masked_res = median(out_im_masked - out_sky, 5)

 	im_data=0
	nsky_data=0	

end
