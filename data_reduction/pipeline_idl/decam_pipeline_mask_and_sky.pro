if recipe EQ 'create mask' then begin

	tile_uniq=input_target[uniq(input_target.tile, sort(input_target.tile))].tile
	filter_uniq=input_target[uniq(input_target.filter, sort(input_target.filter))].filter
	type_uniq=input_target[uniq(input_target.type, sort(input_target.type))].type
	swarp_resample_dir='temp'
	swarp_vmem_dir='/Volumes/Q6/rmunoz/temp'

	if not file_test(swarp_resample_dir, /directory) then file_mkdir, swarp_resample_dir
	if not file_test(swarp_vmem_dir, /directory) then file_mkdir, swarp_vmem_dir

	decam_data
	for i=0L, n_elements(tile_uniq)-1 do begin
		for k=0L, n_elements(type_uniq)-1 do begin

			for j=0L, n_elements(filter_uniq)-1 do begin

				stack_im_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+'.fits'
				stack_weight_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+'.WEIGHT.fits'
				stack_check_file=output_stack_check_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+'.CHECK_SEGMENTATION.fits'
				stack_mask_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+'.MASK.fits'

				if do_overwrite then file_delete, stack_mask_file, /allow_non
				print, 'CREATE MASK - Processing object mask file '+stack_mask_file

				if file_test(stack_mask_file, /regular) EQ 0 AND file_test(stack_check_file, /regular) EQ 1 then begin

					print, 'CREATE MASK - Creating object mask file '+stack_mask_file

					fits_open, stack_check_file, fcb;, /fpack
					im_size=fcb.axis[0:1]
					im_h=headfits(stack_check_file)

      		shmmap, 'sm_im_data_tile'+tile_uniq[i], /float, dim=im_size, /destroy

					im_data=make_array(im_size, value=0., /long)	
					temp_lines=ceil(im_size[1]/5.);  4000L ;2000L
					temp_max=long(product(im_size)-1.)
					ii=0L
					repeat begin
						print, 'Reading SEGMENTATION file - Iteration ', strn(ii+1), ' of ', strn(ceil(im_size[1]*1./temp_lines))
						ii1=long(ii*temp_lines*im_size[0])
						ii2=long((ii+1)*temp_lines*im_size[0] - 1) < temp_max
						fits_read, fcb, temp_data, temp_h, first=ii1, last=ii2, exten_no=0
				    im_data[*,ii1/im_size[0]:(ii2+1)/im_size[0]-1]=reform(temp_data, [im_size[0],(ii2-ii1+1)/im_size[0]])
				    ii++
				  endrep until long(ii*temp_lines*im_size[0]) GT temp_max
		  		fits_close, fcb

    			temp_data=shmvar('sm_im_data_tile'+tile_uniq[i])
    			temp_data[0,0]=(im_data GT 0)*1.
					pout=ptr_new( (im_data GT 0) )
					temp_data=0
					im_data=0

					temp_time=systime(1)
					if n_elements(bridges) GT 0 then burn_bridges, bridges
					bridges = build_bridges(do_n_cpu)

					for ii=0L, n_elements(bridges)-1 do begin
						(bridges[ii])->execute, '.r worker_convolve'
				    (bridges[ii])->execute, '.r callback_convolve'
  				  (bridges[ii])->setproperty, callback='callback_convolve'
				    (bridges[ii])->execute, "SHMMap, 'sm_im_data_tile"+tile_uniq[i]+"', /float,  dim=["+strjoin(strtrim(im_size,2)+'L',',')+"]"
				    (bridges[ii])->execute, "SHMMap, 'gauss_kernel_small', /float,  dim=[51L,51L]"
  				  (bridges[ii])->execute, "SHMMap, 'gauss_kernel_medium', /float,  dim=[301L,301L]"
					endfor

					power_x=alog(im_size[0])/alog(2.)
					power_y=alog(im_size[1])/alog(2.)
					print, 'Power of two for image size along X and Y: ', power_x, power_y
					power_x=floor(power_x)-1
					power_y=floor(power_y)-1

					temp_xrange = [ 2L^power_x*indgen(ceil(im_size[0]*1./2L^power_x))-[0,200+0,400+200,400+600], im_size[0] ]
					temp_yrange = [ 2L^power_y*indgen(ceil(im_size[1]*1./2L^power_y))-[0,200+0,400+200,400+600], im_size[1] ]

					print, 'Image size ', im_size
					print, 'Convolve xrange ', temp_xrange
					print, 'Convolve yrange ', temp_yrange
;			    temp_xrange=round( im_size[0]*[0.,0.25,0.5,0.75,1.] )
;  			  temp_yrange=round( im_size[1]*[0.,0.25,0.5,0.75,1.] )

    			for ii=0L, n_elements(temp_xrange)-2 do begin
			      for jj=0L, n_elements(temp_yrange)-2 do begin
      			  print, 'Iteration '+strn(ii)+','+strn(jj)

    			    xrange_section=[(temp_xrange[ii]-200)>0L,(temp_xrange[ii+1]-1+200)<(im_size[0]-1)]
      			  yrange_section=[(temp_yrange[jj]-200)>0L,(temp_yrange[jj+1]-1+200)<(im_size[1]-1)]
      			  ud = {xrange:temp_xrange[ii:ii+1], yrange:temp_yrange[jj:jj+1], pout:pout, im_size:im_size}
      			  bridge = get_idle_bridge(bridges)
     			  	bridge->setproperty, userdata=ud
      			  bridge->setvar, 'tile', tile_uniq[i]
      			  bridge->setvar, 'im_size', im_size
    			    bridge->setvar, 'xrange', xrange_section
       				bridge->setvar, 'yrange', yrange_section
        			bridge->execute, /nowait, 'worker_convolve, tile, im_size, xrange, yrange, out'

			      endfor
   				endfor
    			barrier_bridges, bridges

			    for ii=0L,n_elements(bridges)-1 do begin
			      (bridges[ii])->execute, "SHMunmap, 'sm_im_data_tile"+tile_uniq[i]+"'"
   					(bridges[ii])->execute, "SHMunmap, 'gauss_kernel_small'"
   					(bridges[ii])->execute, "SHMunmap, 'gauss_kernel_medium'"
    			endfor
					shmunmap, 'sm_im_data_tile'+tile_uniq[i]
 
					print
    			print, 'CREATE MASK - Partial time up to Convolution ', strn((systime(1)-temp_time)/60.), ' min'

					sxaddpar, im_h, 'BITPIX', 8, /pdu
					sxaddpar, im_h, 'BZERO', 0, /pdu
    			writefits, stack_mask_file, (*pout) GT 0, im_h
    			ptr_free, pout

					burn_bridges, bridges
					temp=temporary(bridges)
					temp=0

				endif $
				else begin

					print, 'CREATE MASK - Object mask already exists'
				endelse

			endfor

			stack_mask_file=file_search(output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_?_'+type_uniq[k]+'.MASK.fits', COUNT=n_stack_mask)
			stack_master_mask_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+type_uniq[k]+'.MASTER_MASK.fits'
			swarp_list_file=output_stack_swarp_dir+'/swarp_fornax_tile'+tile_uniq[i]+'_'+type_uniq[k]+'.MASTER_MASK.lst'
		
			if do_overwrite then file_delete, stack_master_mask_file, /allow_non

			if file_test(stack_master_mask_file, /regular) EQ 0 then begin

				if n_stack_mask EQ 1 then begin
					command='ln -s '+stack_mask_file+' '+stack_master_mask_file
					print, command
					spawn, command
				endif else $
				if n_stack_mask GT 1 then begin

					forprint, stack_mask_file, textout=swarp_list_file, FORMAT='(A)', /NOCOMMENT
					command='swarp @'+swarp_list_file+' -c swarp_config/ctio_decam.swarp'+' -IMAGEOUT_NAME '+stack_master_mask_file+' -COMBINE_TYPE SUM -WEIGHTOUT_NAME temp.fits -DELETE_TMPFILES Y -RESAMPLING_TYPE NEAREST -OVERSAMPLING 0 -RESAMPLE_DIR '+swarp_resample_dir+' -RESAMPLE_SUFFIX .resample.fits -SATLEV_DEFAULT 35000 -SUBTRACT_BACK N -WEIGHT_TYPE NONE -WEIGHT_SUFFIX .WEIGHT.fits -WEIGHT_THRESH 0. -RESCALE_WEIGHTS N -BLANK_BADPIXELS N -WRITE_XML N -FSCALASTRO_TYPE NONE -FSCALE_KEYWORD NONE -FSCALE_DEFAULT 1. -GAIN_KEYWORD NONE -GAIN_DEFAULT 0. -VERBOSE_TYPE QUIET -VMEM_DIR '+swarp_vmem_dir
					print, command
					spawn, command
  				file_delete, 'temp.fits', /noexpand, /allow_non, /quiet

				endif else begin
					print, 'There are no images to compute the master mask image'
				endelse

			endif

			for j=0L, n_elements(filter_uniq)-1 do begin
				; Now, project the MASK file into the individual images

				gv_im=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.type EQ type_uniq[k], n_gv_im)
;				output_mask_file=repstr(input_target[gv_im].im_file, '.fits', '.MASK.fits')
;				output_mask_head_file=repstr(input_target[gv_im].im_file, '.fits', '.MASK.head')

				temp_time=systime(1)
				if n_elements(bridges) GT 0 then burn_bridges, bridges
				bridges = build_bridges(15)

				for ii=0L, n_elements(bridges)-1 do begin
					(bridges[ii])->execute, '.r worker'
					(bridges[ii])->execute, '.r callback'
					(bridges[ii])->setproperty, callback='callback'
				endfor

				for ii=0L, n_gv_im-1 do begin
					print, 'CREATE MASK - Processing image ', input_target[gv_im[ii]].im_file

					if file_test(input_target[gv_im[ii]].mask_file, /noexpand, /regular) then begin
						;fits_info, input_target[gv_im[ii]].mask_file, n_ext=mask_n_ext, /silent $
						fits_open, input_target[gv_im[ii]].im_file, im_fcb
						fits_close, im_fcb
						fits_open, input_target[gv_im[ii]].mask_file, mask_fcb
						fits_close, mask_fcb

 						do_exist=0
						if mask_fcb.nextend EQ input_target[gv_im[ii]].n_chip then begin
							if total( mask_fcb.axis[0:1,1:-1] EQ im_fcb.axis[0:1,1:-1] ) EQ input_target[gv_im[ii]].n_chip*2 then do_exist=1
						endif
					endif $
					else do_exist=0

					if do_exist EQ 1 AND do_overwrite EQ 0 then begin
						print, 'CREATE MASK - Individual mask already exists'
						continue
					endif

					print, 'Building object mask for file '+input_target[gv_im[ii]].im_file+'  '+strn(ii+1)+'/'+strn(n_gv_im)
					print, 'Tile '+input_target[gv_im[ii]].tile
					print, 'Filter '+input_target[gv_im[ii]].filter

					fits_open, input_target[gv_im[ii]].im_file, fcb
					im_h=headfits(input_target[gv_im[ii]].im_file, exten=0)

					openr, head_lun, input_target[gv_im[ii]].scamp_head_file, /get_lun
					head_scamp=list()
					temp_data= ''
					while not EOF(head_lun) do begin
						readf, head_lun, temp_data
						head_scamp.add, temp_data
					endwhile
					head_scamp=head_scamp.toarray(type='string')
					free_lun, head_lun
					gv_end=where(head_scamp EQ 'END     ')
					gv_end=[-1,gv_end]

					sxaddpar, im_h, 'BITPIX', 8, /pdu
					sxaddpar, im_h, 'BZERO', 0, /pdu
					mwrfits, 3, input_target[gv_im[ii]].mask_file, im_h, /create

					swarp_head_file='temp/swarp_im_tile'+tile_uniq[i]+'_chip'+string(indgen(input_target[gv_im[ii]].n_chip)+1, FORMAT='(I0)')+'.head'
					swarp_im_file='temp/swarp_im_tile'+tile_uniq[i]+'_chip'+string(indgen(input_target[gv_im[ii]].n_chip)+1, FORMAT='(I0)')+'.fits'
					swarp_wim_file='temp/swarp_wim_tile'+tile_uniq[i]+'_chip'+string(indgen(input_target[gv_im[ii]].n_chip)+1, FORMAT='(I0)')+'.fits'
					swarp_resample_suffix='.tile'+tile_uniq[i]+'_chip'+string(indgen(input_target[gv_im[ii]].n_chip)+1, FORMAT='(I0)')+'.fits'

					for jj=0L, input_target[gv_im[ii]].n_chip-1 do begin
						im_h=headfits(input_target[gv_im[ii]].im_file, exten=jj+1)
						im_size=fcb.axis[0:1,jj+1]

						mkhdr, temp_h, 1, im_size, /IMAGE
						gv=where(stregex(temp_h,'^END.*') EQ 0, n_gv)
						head_mask=list()
						head_mask.add, temp_h[0:gv[0]-1], /extract
						head_mask.add, head_scamp[gv_end[jj]+1:gv_end[jj+1]], /extract
						head_mask=head_mask.toarray(type='string')

						openw, head_lun, swarp_head_file[jj], /GET_LUN
						printf, head_lun, head_mask, FORMAT='(A)'
						close, head_lun, /all

						command='swarp '+stack_master_mask_file+' -c swarp_config/ctio_decam.swarp'+' -IMAGEOUT_NAME '+swarp_im_file[jj]+' -WEIGHTOUT_NAME '+swarp_wim_file[jj]+' -DELETE_TMPFILES Y -RESAMPLING_TYPE NEAREST -OVERSAMPLING 0 -RESAMPLE_DIR '+swarp_resample_dir+' -RESAMPLE_SUFFIX '+swarp_resample_suffix[jj]+' -SATLEV_DEFAULT 35000 -SUBTRACT_BACK N -WEIGHT_TYPE NONE -WEIGHT_SUFFIX .WEIGHT.fits -WEIGHT_THRESH 0. -RESCALE_WEIGHTS N -BLANK_BADPIXELS N -WRITE_XML N -FSCALASTRO_TYPE NONE -FSCALE_KEYWORD NONE -FSCALE_DEFAULT 1. -GAIN_KEYWORD NONE -GAIN_DEFAULT 0. -VERBOSE_TYPE QUIET -MEM_MAX 3072 -COMBINE_BUFSIZE 3072 -VMEM_DIR '+swarp_vmem_dir

						bridge = get_idle_bridge(bridges)
						bridge->setvar, 'command', command

						print, 'CREATE MASK - Running command ', command
						bridge->execute, /nowait, 'worker, command, out'

					endfor
					fits_close, fcb
	
					barrier_bridges, bridges

					for jj=0L, input_target[gv_im[ii]].n_chip-1 do begin

						im_data=(readfits(swarp_im_file[jj], /silent) GT 0)
						im_h=headfits(input_target[gv_im[ii]].im_file, exten=jj+1)
						sxaddpar, im_h, 'BITPIX', 8, /pdu
						sxaddpar, im_h, 'BZERO', 0, /pdu
						mwrfits, im_data, input_target[gv_im[ii]].mask_file, im_h, /silent

					endfor
					print, 'CREATE MASK - Individual mask image was built ', input_target[gv_im[ii]].mask_file + string(10b)

					file_delete, swarp_head_file
					file_delete, swarp_im_file
					file_delete, swarp_wim_file
	

				endfor

				burn_bridges, bridges
				; If IDL hangs on burn_bridges, then execute the following command in the terminal "pkill -f idl_opserver"
				temp=temporary(bridges)
				temp=0
    		
				print, 'CREATE MASK - Partial time for creating individual masks ', strn((systime(1)-temp_time)/60.), ' min'

			endfor
		endfor
	endfor

	shmunmap, 'gauss_kernel_small'
	shmunmap, 'gauss_kernel_medium'
	shmunmap, 'gauss_kernel_large'


endif else $
if recipe EQ 'sky subtraction' then begin


	for i=0L, n_elements(input_target)-1 do begin
			
		print, 'SKY SUBTRACTION - Iteration '+strn(i+1)+' of '+strn(n_elements(input_target))
		print, 'Processing image ', input_target[i].im_file
		temp_time1=systime(1)

		case do_sky_method of
			'median running': begin
				do_worker='worker_sky_median_running'
				do_callback='callback_sky_median_running'
				ss_im_file=repstr(input_target[i].ss_im_file,'.fits','.001.fits')
				sky_im_file=repstr(input_target[i].sky_im_file,'.fits','.001.fits')
				end
			'median global': begin
				do_worker='worker_sky_median_global'
				do_callback='callback_sky_median_global'
				ss_im_file=repstr(input_target[i].ss_im_file,'.fits','.002.fits')
				sky_im_file=repstr(input_target[i].sky_im_file,'.fits','.002.fits')
				end
			'surface': begin
				do_sky_nim=0
				do_worker='worker_sky_surface'
				do_callback='callback_sky_surface'
				masked_im_file=repstr(input_target[i].masked_im_file,'.fits','.003.fits')
				masked_median_im_file=repstr(input_target[i].masked_median_im_file,'.fits','.003.fits')
				masked_res_im_file=repstr(input_target[i].masked_res_im_file,'.fits','.003.fits')
				ss_im_file=repstr(input_target[i].ss_im_file,'.fits','.003.fits')
				sky_im_file=repstr(input_target[i].sky_im_file,'.fits','.003.fits')
				end
			else: stop
		endcase

		if file_test(ss_im_file, /noexpand, /regular) then begin
			fits_open, input_target[i].im_file, im_fcb
			fits_close, im_fcb
			fits_open, ss_im_file, ss_fcb
			fits_close, ss_fcb

 			do_exist=0
			if ss_fcb.nextend EQ input_target[i].n_chip then begin
				if total( ss_fcb.axis[0:1,1:-1] EQ im_fcb.axis[0:1,1:-1] ) EQ input_target[i].n_chip*2 then do_exist=1
			endif
		endif $
		else do_exist=0

		if do_exist EQ 1 AND do_overwrite EQ 0 then begin
			print, 'SKY SUBTRACTION - Individual sky subtracted image already exists'
			continue
		endif

		fits_open, input_target[i].im_file, fcb
		fits_close, fcb
		im_size=fcb.axis[0:1,1]
		im_chip_n=fcb.nextend

		im_chip_range=indgen(ceil(im_chip_n*1./do_sky_nchip))*do_sky_nchip
		pout_im=ptr_new( fltarr([im_size,im_chip_n]) )
;		pout_im_masked=ptr_new( fltarr([im_size,im_chip_n]) )
;		pout_im_masked_median=ptr_new( fltarr([im_size,im_chip_n]) )
;		pout_im_masked_res=ptr_new( fltarr([im_size,im_chip_n]) )
		pout_sky=ptr_new( fltarr([im_size,im_chip_n]) )

		if do_sky_nim GT 0 then begin

			if do_sky_use_target then begin
				gv_sky=where( input_target.filter EQ input_target[i].filter AND input_target.program EQ input_target[i].program, n_gv_sky )
			endif else begin
				if input_target[i].filter EQ 'u' then begin
					gv_sky=where( input_target.tile	NE do_sky_tile_exclude AND input_target.filter EQ input_target[i].filter AND input_target.program EQ input_target[i].program, n_gv_sky )
				endif else begin
					gv_sky=where( input_target.tile	NE do_sky_tile_exclude AND input_target.tile NE input_target[i].tile AND input_target.filter EQ input_target[i].filter AND input_target.program EQ input_target[i].program, n_gv_sky )
				endelse
			endelse

			n_sky=do_sky_nim < n_gv_sky
			gv_sky=(gv_sky[sort( abs(input_target[gv_sky].mjd-input_target[i].mjd) )])[0:n_sky-1]
			sm_sky_name='sm_im_'+['sky1','sky2','sky3','sky4','sky5','sky6','sky7','sky8','sky9','sky10']
			sm_sky_name=sm_sky_name[0:n_sky-1]

			print, 'SKY MEDIAN - Using the following images to compute the sky ', strn(n_sky)
			print, input_target[gv_sky].im_file, FORMAT='(A)'
			print
		endif else begin

			n_sky=0
			temp=temporary(sm_sky_name)
		endelse

		for j=0L, n_elements(im_chip_range)-1 do begin

			chip_offset=im_chip_range[j]
			chip_n=j LT n_elements(im_chip_range)-1 ? (im_chip_range[j+1]-chip_offset) : im_chip_n-chip_offset

			print, 'SKY SUBTRACTION - Using the following chip range ', strn(chip_offset),+' - '+strn(chip_offset+chip_n)
			print

			im_data=fltarr([im_size,chip_n])
			mask_data=bytarr([im_size,chip_n])
			sm_im_name='sm_im_target'
			sm_mask_name='sm_mask_target'
	
		  shmmap, sm_mask_name, /byte, dim=[im_size,chip_n], /destroy
			for k=0L, (1+n_sky)-1 do begin

				if k EQ 0 then begin
					gv=i
					sm_name=sm_im_name
				endif $
				else begin
					gv=gv_sky[k-1]
					sm_name=sm_sky_name[k-1]
				endelse	

				print, 'SKY MEDIAN - Reading image into memory ', input_target[gv].im_file
				temp_time2=systime(1)
		
		  	shmmap, sm_name, /float, dim=[im_size,chip_n], /destroy
	
				for l=0L, chip_n-1 do begin ;input_target[gv_input[k]].n_chip-1 do begin
					im_data[*,*,l]=readfits(input_target[gv].im_file, im_h, ext=chip_offset+l+1, /silent)
					mask_data[*,*,l]=readfits(input_target[gv].mask_file, im_h, ext=chip_offset+l+1, /silent)
				endfor
	
				if k EQ 0 then begin
		    	temp_data=shmvar(sm_mask_name)
	 	 	  	temp_data[0,0,0]=mask_data
					temp_data=0
				endif else begin
					gv_object=where(mask_data NE 0, n_gv_object)
					if n_gv_object GT 0 then im_data[gv_object]=!VALUES.F_NAN
				endelse
	
		    temp_data=shmvar(sm_name)
	 	 	  temp_data[0,0,0]=im_data
				temp_data=0
	
			endfor

			if n_elements(bridges) GT 0 then burn_bridges, bridges
			bridges = build_bridges(do_n_cpu)
			im_data=0
			mask_data=0
	
			for ii=0L, n_elements(bridges)-1 do begin
				(bridges[ii])->execute, '.r '+do_worker ;'worker_sky_median'
				(bridges[ii])->execute, '.r '+do_callback  ;'callback_sky_median'
	  		(bridges[ii])->setproperty, callback=do_callback  ;'callback_sky_median'
				
				(bridges[ii])->execute, "SHMMap, 'sm_im_target', /float,  dim=["+strjoin(strtrim([im_size,chip_n],2)+'L',',')+"]"
				(bridges[ii])->execute, "SHMMap, 'sm_mask_target', /byte,  dim=["+strjoin(strtrim([im_size,chip_n],2)+'L',',')+"]"
				for k=0L, n_elements(sm_sky_name)-1 do begin
					(bridges[ii])->execute, "SHMMap, '"+sm_sky_name[k]+"', /float,  dim=["+strjoin(strtrim([im_size,chip_n],2)+'L',',')+"]"
				endfor
			endfor
	
	
			print
	    print, 'SKY SUBTRACTION - Partial time for reading the target and sky images ', strn((systime(1)-temp_time2)/60.), ' min'
	
;			cgloadct, 0
;			cgwindow, wxsize=1200, wysize=1200, winid=1
			loadct, 0

			wBase = WIDGET_BASE(/COLUMN)
			wDraw = WIDGET_WINDOW(wbase, UVALUE='draw', UNAME='DRAW')
			WIDGET_CONTROL, wBase, /REALIZE
			WIDGET_CONTROL, wDraw, GET_VALUE=plot_window

;			plot_window.Select

;			window, 1, xsize=1200, ysize=1600
			plot_pos_all = cgLayout([2,8], OXMargin=[6,4], OYMargin=[2,2], XGap=4, YGap=2)
		
			for l=0L, chip_n-1 do begin ;input_target[i].n_chip-1 do begin
				chip_name=fcb.extname[chip_offset+l+1]
				plot_do = (chip_name EQ 'S4' OR chip_name EQ 'N4' OR chip_name EQ 'N15' OR chip_name EQ 'S18') ? 1 : 0
				case chip_name of
      		'S4': plot_pos_chip=plot_pos_all[*,0:3]
      		'N4': plot_pos_chip=plot_pos_all[*,4:7]
      		'S18': plot_pos_chip=plot_pos_all[*,8:11]
      		'N15': plot_pos_chip=plot_pos_all[*,12:15]
					else: plot_pos_chip=[0,0,0,0]
    		endcase

	    	ud = ptr_new({chip:l, chip_name:chip_name, im_size:im_size, chip_offset:chip_offset, chip_n:chip_n, pout_im:pout_im, pout_sky:pout_sky, plot_do:plot_do, plot_window:plot_window})

	    	bridge = get_idle_bridge(bridges)
	     	bridge->setproperty, userdata=ud
	      bridge->setvar, 'chip', l
	      bridge->setvar, 'chip_name', chip_name
	      bridge->setvar, 'im_size', im_size
	      bridge->setvar, 'chip_offset', chip_offset
	      bridge->setvar, 'chip_n', chip_n
	      bridge->setvar, 'im_file', input_target[i].im_file
	      bridge->setvar, 'im_name', sm_im_name
	      bridge->setvar, 'sky_file', n_sky GT 0 ? input_target[gv_sky].im_file : ''
	      bridge->setvar, 'sky_name', n_sky GT 0 ? sm_sky_name : ''
	      bridge->setvar, 'mask_file', input_target[i].mask_file
	      bridge->setvar, 'mask_name', sm_mask_name
	      bridge->setvar, 'im_filter', input_target[i].filter
	      bridge->setvar, 'program', input_target[i].program
	      bridge->setvar, 'plot_pos_chip', plot_pos_chip
	      bridge->setvar, 'sky_region', do_sky_region
	      bridge->execute, nowait=1, do_worker+', chip, chip_name, im_size, chip_offset, chip_n, im_file, im_name, sky_file, sky_name, mask_file, mask_name, im_filter, program, sky_region, out_im, out_sky, plot_pos_chip'
	
			endfor
	    barrier_bridges, bridges

;			cgcontrol, output='results/sky_subtraction_003_'+program+'_'+repstr(im_file[-1],'.fits','')+'.pdf'
;  		cgdelete, /all

			for ii=0L,n_elements(bridges)-1 do begin
				for jj=0L, n_elements(sm_sky_name)-1 do begin
					(bridges[ii])->execute, "SHMunmap, '"+sm_sky_name[jj]+"'"
				endfor
				(bridges[ii])->execute, "SHMunmap, '"+sm_im_name+"'"
				(bridges[ii])->execute, "SHMunmap, '"+sm_mask_name+"'"
	    endfor
	
			for jj=0L, n_elements(sm_sky_name)-1 do begin
				shmunmap, sm_sky_name[jj]
			endfor
			shmunmap, sm_im_name ;'sm_mask_data'
			shmunmap, sm_mask_name ;'sm_mask_data'
	
			burn_bridges, bridges
			temp=temporary(bridges)
			temp=0

		endfor
						
		im_h=headfits(input_target[i].im_file, exten=0)
;		mwrfits, 3, masked_im_file, im_h, /create
;		mwrfits, 3, masked_res_im_file, im_h, /create
;		mwrfits, 3, masked_median_im_file, im_h, /create
		mwrfits, 3, ss_im_file, im_h, /create
;		mwrfits, 3, sky_im_file, im_h, /create
	
		for l=0L, input_target[i].n_chip-1 do begin
			im_h=headfits(input_target[i].im_file, exten=l+1)
;			mwrfits, (*pout_im_masked)[*,*,l], masked_im_file, im_h, /silent
;			mwrfits, (*pout_im_masked_res)[*,*,l], masked_res_im_file, im_h, /silent
;			mwrfits, (*pout_im_masked_median)[*,*,l], masked_median_im_file, im_h, /silent
			mwrfits, (*pout_im)[*,*,l], ss_im_file, im_h, /silent
;			mwrfits, (*pout_sky)[*,*,l], sky_im_file, im_h, /silent
		endfor
		print, 'SKY MEDIAN - Image was saved on disk ', ss_im_file

	  ptr_free, pout_im
;	  ptr_free, pout_im_masked
;	  ptr_free, pout_im_masked_median
;	  ptr_free, pout_im_masked_res
	  ptr_free, pout_sky

		print
	  print, 'SKY MEDIAN - Total execution time for reading, computing the median and subtracting ', strn((systime(1)-temp_time1)/60.), ' min'
		WIDGET_CONTROL, wBase, /DESTROY
		if do_debug then stop
	
	endfor

endif else $

