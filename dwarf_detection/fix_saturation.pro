
if recipe EQ 'dwarf galaxy detection' then begin

  tile_uniq=input_target[uniq(input_target.tile, sort(input_target.tile))].tile
  filter_uniq=input_target[uniq(input_target.filter, sort(input_target.filter))].filter
  type_uniq=input_target[uniq(input_target.type, sort(input_target.type))].type

	locus_mag_range=[16,24]
	locus_flux_radius_min=1.
	locus_ellipticity_min=0.6
	locus_flag_max=3
	locus_radius_bin=0.5
	plot_mag_range=[28,15]

	vig_diam=101

  for i=0L, n_elements(tile_uniq)-1 do begin
    for j=0L, n_elements(filter_uniq)-1 do begin
      for k=0L, n_elements(type_uniq)-1 do begin

  			if do_align_filter NE '' then begin
	 			 if filter_uniq[j] NE do_align_filter then do_align='_ALIGN'+do_align_filter $
   			 else do_align=''
  			endif

        if do_sky_method EQ 'median global' then begin
          swarp_list_file=output_stack_swarp_dir+'/swarp_ss_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.dat'
          swarp_xml_file=output_stack_swarp_dir+'/swarp_ss_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.xml'
          swarp_im_out=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.fits'
          swarp_weight_out=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.WEIGHT.fits'
          sex_stack_cat_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.ldac'
          sex_stack_xml_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.xml'
          sex_stack_checkimage_file=output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.CHECK_SEGMENTATION.fits'
        endif else $
        if do_sky_method EQ 'gradient tps' then begin
          swarp_list_file=output_stack_swarp_dir+'/swarp_ss_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.dat'
          swarp_xml_file=output_stack_swarp_dir+'/swarp_ss_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.xml'

          stack_im_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.fits'
          stack_weight_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.WEIGHT.fits'
					stack_im_fix_file=output_stack_check_dir+'/ss_fornax_FIX.003.fits'
					stack_weight_fix_file=output_stack_check_dir+'/ss_fornax_FIX.003.WEIGHT.fits'
					stack_im_convol_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_CONVOL.003.fits'

          sex_stack_cat_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.ldac'
          sex_stack_xml_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.xml'
          sex_stack_checkimage_file=output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.CHECK_SEGMENTATION.fits'

					stack_im_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.fits'
					stack_weight_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.WEIGHT.fits'
          stack_cat_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.ldac'
					stack_xml_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.xml'
          stack_checkimage_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.CHECK_SEGMENTATION.fits'
					stack_checkimage_sex_type='SEGMENTATION'

          sex_stack_cat_fix_file=output_stack_check_dir+'/ss_fornax_FIX.003.ldac'
					sex_stack_xml_fix_file=output_stack_check_dir+'/ss_fornax_FIX.003.xml'
          sex_stack_checkimage_fix_file=output_stack_check_dir+'/ss_fornax_FIX.003.CHECK_SEGMENTATION.fits'
					sex_stack_checkimage_type='SEGMENTATION'
        endif else begin
          swarp_list_file=output_stack_swarp_dir+'/swarp_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.dat'
          swarp_xml_file=output_stack_swarp_dir+'/swarp_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.xml'
          swarp_im_out=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.fits'
          swarp_weight_out=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.WEIGHT.fits'
          sex_stack_cat_file=output_stack_sex_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.ldac'
          sex_stack_xml_file=output_stack_sex_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.xml'
          sex_stack_checkimage_file=output_stack_check_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.CHECK_SEGMENTATION.fits'
        endelse

				print
				print, 'DWARF GALAXY DETECTION - Processing image ', stack_im_file

				if file_test(sex_stack_cat_fix_file) EQ 0 OR do_overwrite then begin

					im_data=readfits_big(stack_im_file, im_h)
					wim_data=readfits_big(stack_weight_file, wim_h)
					im_size=size(im_data, /dim)
	
					nim_region=[ im_size[0]/2+[-1000,1000] > 0 < (im_size[0]-1), im_size[1]/2+[-1000,1000] > 0 < (im_size[1]-1) ]
					nim_data=im_data[nim_region[0]:nim_region[1], nim_region[2]:nim_region[3]]
					nwim_data=wim_data[nim_region[0]:nim_region[1], nim_region[2]:nim_region[3]]
					nim_size=size(nim_data, /dim)
					writefits, stack_im_sex_file, nim_data, im_h
					writefits, stack_weight_sex_file, nwim_data, im_h
	
					command='sex '+stack_im_sex_file+' -c sex_config/ctio_decam_stack.sex -CATALOG_NAME '+stack_cat_sex_file+' -WEIGHT_IMAGE '+stack_weight_sex_file+' -MAG_ZEROPOINT 30. -XML_NAME '+stack_xml_sex_file+' -CHECKIMAGE_TYPE '+stack_checkimage_sex_type+' -CHECKIMAGE_NAME '+stack_checkimage_sex_file+' -BACK_SIZE 128 -FILTER_NAME sex_config/gauss_2.5_5x5.conv -DEBLEND_NTHRESH 16 -DEBLEND_MINCONT 0.01'
					print, command
					spawn, command
	
		      cat_sex=mrdfits(stack_cat_sex_file, 2, cat_sex_h, COLUMNS=['NUMBER','X_IMAGE','Y_IMAGE','FLUX_RADIUS','MAG_AUTO','MAGERR_AUTO','FLAGS','A_IMAGE','B_IMAGE'], /silent)
	
		      wset, 0
	  	    plot, cat_sex.flux_radius, cat_sex.mag_auto, psym=1, xrange=[1,6], yrange=plot_mag_range
	    	  oplot, [0,100], locus_mag_range[0]*[1,1], line=2, color=100
	      	oplot, [0,100], locus_mag_range[1]*[1,1], line=2, color=100
	     		gv_stars=where(cat_sex.mag_auto GT locus_mag_range[0] and cat_sex.mag_auto LT locus_mag_range[1] AND cat_sex.flux_radius GT locus_flux_radius_min AND cat_sex.flags LE locus_flag_max AND cat_sex.b_image/cat_sex.a_image GT locus_ellipticity_min, n_gv_stars)
	     	 	if n_gv_stars LT 10 then begin
	      	  print, 'IQ - Error, there is not enough number of stars'
	     	  	continue
	    	  endif
	      	oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=100
	
		      plothist, cat_sex[gv_stars].flux_radius, temp_xhist, temp_yhist, bin=locus_radius_bin, /noplot
	  	    temp=max(temp_yhist, gv) & locus_radius=temp_xhist[gv]
	    	  locus_radius=median((cat_sex[where(cat_sex.mag_auto GT locus_mag_range[0] AND cat_sex.mag_auto LT locus_mag_range[1] AND cat_sex.flux_radius GT locus_radius*0.9 AND cat_sex.flux_radius LT locus_radius*1.1 , n_gv)]).flux_radius)
	   	  	print, locus_radius, FORMAT='("Median flux_radius ", F0.1)'
	
	      	oplot, locus_radius*[1,1], [0,100], color=200
	      	oplot, locus_radius*[0.9,0.9], [0,100], line=2, color=200
	      	oplot, locus_radius*[1.1,1.1], [0,100], line=2, color=200
	
					gv_stars=where(cat_sex.mag_auto GT locus_mag_range[0] AND cat_sex.mag_auto LT locus_mag_range[1]-1 AND cat_sex.flux_radius GT locus_radius*0.9 AND cat_sex.flux_radius LT locus_radius*1.1 AND cat_sex.flags LE locus_flag_max AND cat_sex.x_image GT vig_diam/2. AND cat_sex.x_image LT (nim_size[0]-vig_diam/2.) AND cat_sex.y_image GT vig_diam/2. AND cat_sex.y_image LT (nim_size[1]-vig_diam/2.), n_gv_stars)
					oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=200
	
					gv_sort=reverse(sort(cat_sex[gv_stars].mag_auto))
					gv_stars=gv_stars[gv_sort]
	
	        temp_n=n_gv_stars<5
	        temp_x=cat_sex[gv_stars].x_image
	        temp_y=cat_sex[gv_stars].y_image
					psf_fwhm=list()
	
	        for k=0L, temp_n-1 do begin
	          x_range=[ floor(temp_x[k]-(vig_diam-1.)/2), ceil(temp_x[k]+(vig_diam-1.)/2) ]
	          y_range=[ floor(temp_y[k]-(vig_diam-1.)/2), ceil(temp_y[k]+(vig_diam-1.)/2) ]
	          vig_data = nim_data[x_range[0]:x_range[1],y_range[0]:y_range[1]]; - im_sky
	          vig_size=size(vig_data, /dim)
	
	          x_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
	          y_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
	          vig_center_data=vig_data[x_center_range[0]:x_center_range[1],y_center_range[0]:y_center_range[1]]
	          vig_center_size=size(vig_center_data, /dim)
	
						gcntrd, vig_center_data, vig_center_size[0]/2., vig_center_size[1]/2., vig_cx, vig_cy, 3.
						im_c = [x_center_range[0]+vig_cx, y_center_range[0]+vig_cy]
	;          im_max=max(vig_center_data, gv_max)
	;          im_c= [gv_max mod vig_center_size[0], gv_max/vig_center_size[0]]
	;          im_c += [x_center_range[0], y_center_range[0]]
	
						dist_circle, vig_mask, vig_size, im_c[0], im_c[1]
						vig_sky=median(vig_data[where(vig_mask GT 20., n_sky)])
						vig_mag=25.-2.5*alog10(total(vig_data[where(vig_mask LE 20., n_star)]) - n_star*vig_sky )
	
						vig_res=abs( vig_data-vig_sky - max(vig_center_data-vig_sky)/2. )
						gv=sort(vig_res*vig_mask^2)
						vig_fwhm=2*median(vig_mask[gv[1:5]])
						psf_fwhm.add, vig_fwhm
						;print, 'Radius for computing FWHM ', vig_mask[gv[1:5]]
						plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
						oplot, [0,100], max(vig_data-vig_sky)/2.*[1,1], line=2, color=200
						oplot, vig_fwhm/2.*[1,1], [-1e5,1e5], line=2, color=200
	
					endfor
	
					vig_fwhm=median(psf_fwhm.toarray(type='float'))
	   	  	print, vig_fwhm, FORMAT='("Median FWHM ", F0.2)'
	
					vig_skyrad=4*vig_fwhm < 40.
					vig_psfrad=3*vig_fwhm < 50.
					vig_fitrad=vig_fwhm < 50.
	
					getpsf, nim_data, cat_sex[gv_stars[0:9]].x_image-1, cat_sex[gv_stars[0:9]].y_image-1, cat_sex[gv_stars[0:9]].mag_auto, indgen(10)*0., decam_ron, decam_gain, psf_param, psf_residuals, indgen(10), vig_psfrad, vig_fitrad, output_stack_check_dir+'/im_sextractor_psf.fits' 
					psf_sigma=sqrt( (psf_param[3]^2 + psf_param[4]^2)/2 )
	   	  	print, psf_sigma, FORMAT='("Median PSF SIGMA ", F0.2)'
	
					; Now we use the weightmap to look for blobs				
					gv_blank=where(wim_data EQ 0., n_gv_blank)
					im_size=size(wim_data, /DIM)
					im_blob=Obj_New('Blob_Analyzer', (wim_data EQ 0.))
					im_blob_stats=replicate({id:0L, npix:0L, indices:list()}, im_blob->numberofblobs())
					for ii=0L, n_elements(im_blob_stats)-1 do begin
						im_blob_stats[ii].id=ii
						im_blob_stats[ii].npix=n_elements( im_blob->getindices(ii) )
					endfor
					gv_sort=reverse(sort(im_blob_stats.npix))
					im_blob_stats=im_blob_stats[gv_sort]
					gv_blob=where(im_blob_stats.npix GT 1 AND im_blob_stats.npix LT 1e5, n_gv_blob)
					gv_border=im_blob->getindices( im_blob_stats[0].id )
	
					wim_data[gv_blank]=1.
					wim_data[gv_border]=0.
					writefits_big, stack_weight_fix_file, wim_data, im_h
	
					im_data[gv_blank]=!values.f_nan
					im_crop_border=20
					im_star_radius=50 ; Box size to fit the saturated stars
	
					for ii=0L, n_gv_blob-1 do begin
						if ii mod 10 EQ 0 then print, 'Processing blob ', strn(ii), ' of ', strn(n_gv_blob)
						gv_ind=array_indices(im_size, im_blob->getindices( im_blob_stats[gv_blob[ii]].id ), /dim)
						im_blob_region=[min(gv_ind[0,*]), max(gv_ind[0,*]), min(gv_ind[1,*]), max(gv_ind[1,*])]
						im_crop_region=[min(gv_ind[0,*])-im_crop_border, max(gv_ind[0,*])+im_crop_border, min(gv_ind[1,*])-im_crop_border, max(gv_ind[1,*])+im_crop_border]
	
						im_crop_data=im_data[im_crop_region[0]:im_crop_region[1], im_crop_region[2]:im_crop_region[3]]
						im_crop_size=size(im_crop_data, /dim)
	
						im_crop_datax=total(im_crop_data, 2, /nan)/total(finite(im_crop_data), 2, /nan)
						im_crop_x=im_crop_region[0]+indgen(n_elements(im_crop_datax))
						im_crop_datay=total(im_crop_data, 1, /nan)/total(finite(im_crop_data), 1, /nan)
						im_crop_y=im_crop_region[2]+indgen(n_elements(im_crop_datay))
	
	;					temp=mpfitpeak(im_crop_x, im_crop_datax, gauss_coeff, /gaussian)
	;					if im_crop_size[0] GE 5 then begin
							gauss_kernel_size = 31 < floor(im_crop_size[0]/2.-1)*2+1
							gauss_kernel=gauss1(findgen(gauss_kernel_size),[(gauss_kernel_size-1)/2.,5*psf_sigma,1e2])
	
							im_crop_datax_convol=convol(im_crop_datax, gauss_kernel, /edge_truncate)/1e2
							temp=max(im_crop_datax_convol, max_id)
							im_crop_x_max=im_crop_region[0]+max_id
	
							im_crop_center_region=[im_crop_x_max-im_crop_border, im_crop_x_max+im_crop_border, min(gv_ind[1,*])-im_crop_border, max(gv_ind[1,*])+im_crop_border]
							im_crop_center_data=im_data[im_crop_center_region[0]:im_crop_center_region[1], im_crop_center_region[2]:im_crop_center_region[3]]
							im_crop_datay=total(im_crop_center_data, 1, /nan)/total(finite(im_crop_center_data), 1, /nan)
	
							im_crop_datay_smooth=smooth(im_crop_datay,3, /edge_truncate, /nan)
							gv_peaks=lclxtrem(im_crop_datay_smooth, 10, /maxima, count=n_gv_peaks)
							if n_gv_peaks GT 1 then begin
								if im_crop_datay_smooth[gv_peaks[1]]/im_crop_datay_smooth[gv_peaks[0]] GT 0.2 then begin
									gv_peaks_sort=(gv_peaks[0:1])[sort(gv_peaks[0:1])]
									im_crop_datay_smooth[ gv_peaks_sort[0]+1:gv_peaks_sort[1]-1]=!values.f_nan
									im_crop_y_max=im_crop_center_region[2]+total(indgen(im_crop_size[0])*im_crop_datay_smooth, /nan)/total(im_crop_datay_smooth, /nan) 
								endif $
								else begin
									im_crop_y_max=im_crop_center_region[2]+gv_peaks[0]
								endelse
							endif else $
							if n_gv_peaks EQ 1 then begin
								im_crop_y_max=im_crop_center_region[2]+gv_peaks
							endif $
							else begin
								print, 'WARNING - No local maxima was found'
								im_crop_y_max=0
							endelse	
	
	;						im_crop_datay_convol=convol(im_crop_datay, gauss_kernel, /edge_truncate)/1e2
	;						temp=max(im_crop_datay_convol, max_id)
	;						im_crop_y_max=im_crop_region[2]+max_id
	;;					endif else begin
	;;						im_crop_datax_convol=0.
	;;					endelse
	
						im_star_center=round([ im_crop_x_max, im_crop_y_max ]) ;mean(im_blob_region[2:3]) ])
	
						cgDisplay, 600, 600, wid=0, location=[0,(get_screen_size())[1]-600.]
						pos = cgLayout([1,2], OXMargin=[5, 3], OYMargin=[5, 8], XGap=5, YGap=5)
						cgplot, im_crop_datax, pos=pos[*,0], xtitle='X (pixels)', charsize=0.8
						cgplot, (im_star_center[0]-im_crop_region[0])*[1,1], [0,1e6], line=2, color='red', /over
						cgplot, im_crop_datay, pos=pos[*,1], xtitle='Y (pixels)', charsize=0.8, /noerase
						cgplot, im_crop_datay_smooth, color='blue', /over
						cgplot, (im_star_center[1]-im_crop_region[2])*[1,1], [0,1e6], line=2, color='red', /over
	
						if max(im_crop_datax_convol) GT 10. AND n_gv_peaks GT 0 then begin
							print, ii, im_star_center[0], im_star_center[1], max(im_crop_datax_convol), FORMAT='("YES saturated star - ID: ", I0, " - Location: ", F0.1, " ", F0.1, " - Max im_convol: ", F0.1)'
	
							im_star_region=[im_star_center[0]+[-1,1]*im_star_radius, im_star_center[1]+[-1,1]*im_star_radius]
							im_star_data=im_data[im_star_region[0]:im_star_region[1], im_star_region[2]:im_star_region[3]]
							im_star_size=size(im_star_data, /dim)
							dist_circle, im_star_mask, im_star_size
							gv=where(finite(im_star_data) AND im_star_mask LT im_star_radius, n_gv)
							bv=where(finite(im_star_data) EQ 0 AND im_star_mask LT im_star_radius, n_bv)
							im_star_ind=array_indices(im_star_data, gv)
							im_star_data_error=1.+im_star_mask
						
	;					gauss2d_coeff=[(im_star_size-1)/2., 2*psf_sigma, 100*max(im_star_data)]
	;					gauss2d_parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 4)
							model2d_coeff=[(im_star_size-1)/2., 2*psf_sigma, 100*max(im_star_data), 2.]
							model2d_parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 5)
							model2d_parinfo[0:2].limited=[1,1]
							model2d_parinfo[0].limits=model2d_coeff[0]+[-5,5]
							model2d_parinfo[1].limits=model2d_coeff[1]+[-5,5]
							model2d_parinfo[2].limits=psf_sigma*[1,5]
							im_star_coeff=mpfit2dfun('MOFFAT2', im_star_ind[0,*], im_star_ind[1,*], im_star_data[gv], im_star_data_error[gv], model2d_coeff, parinfo=model2d_parinfo, /quiet)
	
	;					im_star_model=gauss2( indgen(im_star_size[0])#make_array(im_star_size[1], value=1., /float), make_array(im_star_size[0], value=1., /float)#indgen(im_star_size[1]), im_star_coeff )
							im_star_model=moffat2( indgen(im_star_size[0])#make_array(im_star_size[1], value=1., /float), make_array(im_star_size[0], value=1., /float)#indgen(im_star_size[1]), im_star_coeff )
							temp_data=im_star_data*!VALUES.F_NAN
							temp_data[bv]=im_star_model[bv]
	
							plot_yrange=[0,2*max(im_star_data)]
							plot_zrange=[0,2*max(im_star_data)]
							cgDisplay, 600, 600, wid=1, location=get_screen_size()-600.*[3,1]>0
							cgsurf, im_star_data, xtitle='X', ytitle='Y', ztitle='Flux', rotx=30., rotz=30., color='black', zrange=plot_zrange, /save, /zstyle
							cgsurf, temp_data, color='red', /noerase, /horizontal, zrange=plot_zrange, /t3d, /zstyle
	;					cgsurf, im_star_fit_data, color='blue', /noerase, zrange=[0, 3*max(im_star_data)], /t3d, /zstyle
	;						cgsurf, im_star_model, color='red', /noerase, /horizontal, zrange=[0, 3*max(im_star_data)], /t3d, /zstyle
	
							cgDisplay, 600, 600, wid=2, location=get_screen_size()-600.*[2,1.5]>0
							cgplot, im_star_data[im_star_coeff[0],*], color='black', xtitle='Y', ytitle='Flux', yrange=plot_yrange, charsize=0.8
							cgplot, im_star_model[im_star_coeff[0],*], color='red', /over
	
							temp=max(total(im_star_data, 1, /nan), max_id) 
							cgDisplay, 600, 600, wid=3, location=get_screen_size()-600.*[1,2]>0
							cgplot, im_star_data[*,max_id], color='black', xtitle='X', ytitle='Flux', yrange=plot_yrange, charsize=0.8
							cgplot, im_star_model[*,max_id], color='red', /over
	
							; HERE we replace the missing values with the model values
							im_star_data[bv]=im_star_model[bv]
							im_data[im_star_region[0]:im_star_region[1], im_star_region[2]:im_star_region[3]]=im_star_data
	
							; HERE we interpolate the pixels beyond the star region
	;						temp_ind1=(im_star_region[0]-im_blob_region[0]) GT 0 ? im_blob_region[0]+indgen(im_star_region[0]-im_blob_region[0]) : []
	;						temp_ind2=(im_blob_region[1]-im_star_region[1]) GT 0 ? im_star_region[1]+indgen(im_blob_region[1]-im_star_region[1]) : []
	;						im_blank_ind=[temp_ind1,temp_ind2] ;[im_crop_region[0]+indgen(im_star_region[0]-im_crop_region[0]),im_star_region[1]+indgen(im_crop_region[1]-im_star_region[1])]
	
							im_blank_ind=(im_blob_region[1]-im_blob_region[0]) GT 0 ? (im_blob_region[0]+indgen(im_blob_region[1]-im_blob_region[0])) : [im_blob_region[0]]
							foreach jj, im_blank_ind do begin ;im_crop_region[0]+indgen(im_crop_region[1]-im_crop_region[0]+1) do begin					
								xdata=im_blob_region[2]+indgen(im_blob_region[3]-im_blob_region[2]+1)
								ydata=reform(im_data[jj,im_blob_region[2]:im_blob_region[3]], n_elements(xdata))
								xdata_fit=im_crop_region[2]+indgen(im_crop_region[3]-im_crop_region[2]+1)
								ydata_fit=reform(im_data[jj,im_crop_region[2]:im_crop_region[3]], n_elements(xdata_fit))
								gv=where(finite(ydata_fit), n_gv)
								bv=where(finite(ydata) EQ 0, n_bv)
								if n_bv GT 0 then begin
									if n_gv GE 10 then begin
										im_crop_coeff=robust_poly_fit(xdata_fit[gv], ydata_fit[gv], 2)
										im_data[jj,xdata[bv]] = im_crop_coeff[0] + im_crop_coeff[1]*xdata[bv] + im_crop_coeff[2]*xdata[bv]^2
									endif else $
									if n_gv GE 3 then begin
										im_data[jj,xdata[bv]]=biweight_mean(ydata_fit[gv])
									endif else begin
										im_data[jj,xdata[bv]]=median(ydata_fit[gv])
									endelse
								endif
							endforeach
	
						endif $
						else begin
							print, ii, max(im_crop_datax_convol), FORMAT='("NO saturated star - ID: ", I0, " - Max im_crop_datax_convol: ", F0.1)'
	
							im_blank_ind=(im_blob_region[1]-im_blob_region[0]) GT 0 ? (im_blob_region[0]+indgen(im_blob_region[1]-im_blob_region[0])) : [im_blob_region[0]]
							foreach jj, im_blank_ind do begin ;im_crop_region[0]+indgen(im_crop_region[1]-im_crop_region[0]+1) do begin					
								xdata=im_blob_region[2]+indgen(im_blob_region[3]-im_blob_region[2]+1)
								ydata=reform(im_data[jj,im_blob_region[2]:im_blob_region[3]], n_elements(xdata))
								xdata_fit=im_crop_region[2]+indgen(im_crop_region[3]-im_crop_region[2]+1)
								ydata_fit=reform(im_data[jj,im_crop_region[2]:im_crop_region[3]], n_elements(xdata_fit))
								gv=where(finite(ydata_fit), n_gv)
								bv=where(finite(ydata) EQ 0, n_bv)
								if n_bv GT 0 then begin
									if n_gv GE 10 then begin
										im_crop_coeff=robust_poly_fit(xdata_fit[gv], ydata_fit[gv], 2)
										im_data[jj,xdata[bv]] = im_crop_coeff[0] + im_crop_coeff[1]*xdata[bv] + im_crop_coeff[2]*xdata[bv]^2
									endif else $
									if n_gv GE 3 then begin
										im_data[jj,xdata[bv]]=biweight_mean(ydata_fit[gv])
									endif else begin
										im_data[jj,xdata[bv]]=median(ydata_fit[gv])
									endelse
								endif
							endforeach
	
						endelse
	
					endfor
	
					im_data[gv_border]=0.
					writefits_big, stack_im_fix_file, im_data, im_h
	
					command='sex '+stack_im_fix_file+' -c sex_config/ctio_decam_stack.sex -CATALOG_NAME '+sex_stack_cat_fix_file+' -WEIGHT_IMAGE '+stack_weight_fix_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_fix_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_fix_file+' -BACK_SIZE 128 -FILTER_NAME sex_config/gauss_2.5_5x5.conv -DEBLEND_NTHRESH 16 -DEBLEND_MINCONT 0.01'
					print, command
					spawn, command

				endif
	
				cat_sex=mrdfits(sex_stack_cat_fix_file, 2, cat_sex_h, COLUMNS=['NUMBER','X_IMAGE','Y_IMAGE','ALPHA_J2000','DELTA_J2000','FLUX_RADIUS','MAG_AUTO','ELONGATION'])
				gv_dwarf=where( cat_sex.flux_radius GT 10. AND cat_sex.mag_auto GT 19., n_gv_dwarf, COMPLEMENT=bv_dwarf, NCOMPLEMENT=n_bv_dwarf)

				im_data=readfits_big(stack_im_file, im_h)
				im_size=size(im_data, /dim)
				check_data=long(readfits_big(sex_stack_checkimage_fix_file, im_h))
				temp_h = Histogram(check_data, MIN=1, REVERSE_INDICES=check_ri, BINSIZE=1)

				temp_time=systime(1)
				for ii=0L, n_bv_dwarf-1 do begin
					gv_mask=check_ri[check_ri[cat_sex[bv_dwarf[ii]].number-1]:check_ri[cat_sex[bv_dwarf[ii]].number]-1]
					im_data[gv_mask]=!VALUES.F_NAN
				endfor
    		print, 'Partial time to mask the image ', strn((systime(1)-temp_time)/60.), ' min'

;				dwarf_kernel=make_gauss_kernel(FWHM=[40,40], NPIX=[121,121])
				temp_x=indgen(101)#make_array(101, value=1., /float)
				temp_y=make_array(101, value=1., /float)#indgen(101)
				temp_sigma=40./(2.0d* sqrt( 2.0d* aLog(2.0d) ))
				dwarf_kernel=float(gauss2(temp_x, temp_y, [50., 50., temp_sigma, 1e2]))
				dwarf_kernel_size=size(dwarf_kernel, /dim)

				im_ratio=ceil(im_size*1./dwarf_kernel_size)
				nim_data=make_array(im_ratio*dwarf_kernel_size, value=0., /float)
				nim_data[0,0]=im_data
				nim_data=convolve(nim_data, dwarf_kernel, /no_pad)/total(dwarf_kernel)
    		print, 'Partial time to convolve the image ', strn((systime(1)-temp_time)/60.), ' min'

				print, 'DWARF GALAXY DETECTION - Writing convolved image ', stack_im_convol_file
				writefits_big, stack_im_convol_file, nim_data, im_h
				stop

			endfor
		endfor
	endfor


