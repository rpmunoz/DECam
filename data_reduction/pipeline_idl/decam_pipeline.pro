; Licensed under a 3-clause BSD style license - see LICENSE.rst
; DECam pipeline is being developed by Roberto Pablo Munoz, PhD
; Munoz et al. 2015, ApJ, 813, L15

forward_function decam_data

pro decam_pipeline, recipe, PROGRAM=program, FILTER=filter, OVERWRITE=overwrite, TILE=tile, TYPE=type, DEBUG=debug, STANDARD=standard, AHEAD=ahead, SKY_METHOD=sky_method, SKY_NIM=sky_nim, SKY_NCHIP=sky_nchip, SKY_USE_TARGET=sky_use_target, ALIGN_FILTER=align_filter, NDITHER=ndither

	common decam_share, gauss_kernel_large, gauss_kernel_small

do_program= (n_elements(program) GT 0) ? program : '2013B-0613'
do_filter= (n_elements(filter) GT 0) ? filter : 'g'
do_tile= (n_elements(tile) GT 0) ? tile : '1'
do_type= (n_elements(type) GT 0) ? type : 'long'
do_standard= (n_elements(standard) GT 0) ? standard : 'SDSS'
do_n_chip= (do_program EQ '2013B-0613') ? '61' : '60'
do_scamp_ast='Y'
do_scamp_phot='Y'
do_overwrite= (n_elements(overwrite) GT 0) ? keyword_set(overwrite) : 0
do_ahead= (n_elements(ahead) GT 0) ? keyword_set(ahead) : 0
do_debug= (n_elements(debug) GT 0) ? keyword_set(debug) : 0
do_filter_orig=do_filter
do_program_orig=do_program
do_n_cpu=6
do_ndither= (n_elements(ndither) GT 0) ? ndither : 0
do_sky_nchip= (n_elements(sky_nchip) GT 0) ? sky_nchip : 31
do_sky_nim= (n_elements(sky_nim) GT 0) ? sky_nim : 5
do_sky_method = (n_elements(sky_method) GT 0) ? sky_method : ''
do_sky_tile_exclude =  (n_elements(sky_tile_exclude) GT 0) ? sky_tile_exclude : '1'
do_sky_use_target = (n_elements(sky_use_target) GT 0) ? keyword_set(sky_use_target) : 0
do_sky_region = (n_elements(sky_region) GT 0) ? sky_region : [[0,0,1022,4093],[1023,0,2045,4093]]
do_align_filter = (n_elements(align_filter) GT 0) ? align_filter : ''

input_dir='/Volumes/Q6/NGFS/DECam'
output_dir='/Volumes/Q6/NGFS/DECam'

decam_gain=4.3
decam_ron=6.2

if do_program EQ 'all' then begin
	temp_program=strsplit(file_search(input_dir+'/*', /TEST_DIRECTORY),'/',/extract)
	do_program=strarr(n_elements(temp_program))
	for i=0L, n_elements(temp_program)-1 do do_program[i]=(temp_program[i])[-1]
	gv=where(do_program NE 'stacks', n_gv)
	if n_gv GT 0 then do_program=do_program[gv]
endif

input_im_dir=input_dir+'/'+do_program+'/processed'
input_calib_dir=input_dir+'/'+do_program+'/calib'
input_chip_fwhm=[28,35] ; Where to measure FWHM
input_chip_flux_radius=[27,28,29,34,35,36] ; Where to measure FWHM
input_chip_iq=[27,28,35,36] ; Where to measure FWHM

output_im_dir=output_dir+'/'+do_program+'/pipeline/images'
output_calib_dir=output_dir+'/'+do_program+'/pipeline/calib'
output_sex_dir=output_dir+'/'+do_program+'/pipeline/sextractor'
output_sex_check_dir=output_dir+'/'+do_program+'/pipeline/sextractor/check'
output_scamp_dir=output_dir+'/'+do_program+'/pipeline/scamp'
output_stack_swarp_dir=output_dir+'/stacks';output_dir+'/'+do_program+'/pipeline/swarp'
output_stack_sex_dir=output_dir+'/stacks';output_dir+'/'+do_program+'/pipeline/swarp'
output_stack_psfex_dir=output_dir+'/stacks/psfex' ;output_dir+'/'+do_program+'/pipeline/psfex'
output_stack_check_dir=output_dir+'/stacks/check';output_dir+'/'+do_program+'/pipeline/swarp'

survey_info = list( {tile:'1', coo:['03:38:13.20965','-35:31:53.1400'], filter:['u','g','i']}, $
	{tile:'2', coo:['03:36:41.11946','-33:43:14.5466'], filter:['g','i']}, $
	{tile:'3', coo:['03:29:59.24609','-34:52:23.7734'], filter:['g','i']}, $
	{tile:'4', coo:['03:31:30.90216','-36:40:56.1030'], filter:['u','g','i']}, $
	{tile:'5', coo:['03:39:45.68829','-37:20:26.5726'], filter:['g','i']}, $
	{tile:'6', coo:['03:46:31.32753','-36:11:29.8345'], filter:['g','i']}, $
	{tile:'7', coo:['03:45:00.28988','-34:23:02.7505'], filter:['g','i']}, $
	{tile:'10', coo:['03:28:30.08975','-33:03:57.4157'], filter:['g','i']}, $
	{tile:'13', coo:['03:24:30.95950','-37:49:37.4418'], filter:['g','i']} )

survey_sequence = { seq1:['1','2','3','4','5','6','7','10','13'], seq2:['1','3','4','13','5','6'], seq3:['1','4'] }

survey_cutout=list( {name:'FCC47', tile:'12', ra:51.63327409100075, dec:-35.712338493680505, filter:'g'}, {name:'FCC148', tile:'1', ra:53.8200186062367, dec:-35.265494056646766, filter:'g'}, {name:'FCC170', tile:'1', ra:54.131709345353414, dec:-35.295328080905946, filter:'g'}, {name:'FCC177', tile:'1', ra:54.197351080007834, dec:-34.73808026260103, filter:'g'}, {name:'FCC310', tile:'6', ra:56.55701775714639, dec:-36.69529847345809, filter:'g'} )


for i=0L, n_elements(do_program)-1 do begin

	if file_test('survey_target_'+do_program[i]+'.dat') AND recipe NE 'database' then begin
		readcol, 'survey_target_'+do_program[i]+'.dat', temp_im_orig_file, temp_filter, temp_tile, temp_dither, temp_weight_orig_file, temp_zp, temp_fwhm, temp_mjd, temp_exptime, temp_n_chip, temp_type, FORMAT='A,A,A,A,A,F,F,D,F,I,A', COMMENT='#'
		n_survey_info=n_elements(temp_im_orig_file)
		create_struct, temp_input_target, '', ['im_orig_file','program','filter','tile','dither','weight_orig_file','mjd','mjd_floor','zp','im_file','weight_file','fwhm','exptime','n_chip','sex_cat_file','sex_xml_file','sex_check_file','scamp_cat_file','scamp_head_file','scamp_ahead_file','swarp_im_file','swarp_head_file','sex_zp_file','type','mask_file','masked_im_file','masked_median_im_file','masked_res_im_file','sky_im_file','ss_im_file','ss_weight_file','ss_sex_cat_file','ss_sex_xml_file','ss_scamp_cat_file','ss_scamp_head_file','ss_scamp_ahead_file','ss_swarp_head_file','ss_swarp_im_file'], 'A,A,A,A,A,A,D,D,F,A,A,F,F,I,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A', dim=n_survey_info
		temp_input_target.im_orig_file=input_im_dir[i]+'/'+repstr(temp_im_orig_file,'.fits.fz','.fits')
		temp_input_target.weight_orig_file=input_im_dir[i]+'/'+repstr(temp_weight_orig_file,'.fits.fz','.fits')
		temp_input_target.program=do_program[i]
		temp_input_target.filter=temp_filter
		temp_input_target.tile=strtrim(temp_tile,2)
		temp_input_target.dither=strtrim(temp_dither,2)
		temp_input_target.mjd=temp_mjd
		temp_input_target.mjd_floor=floor(temp_mjd+0.2)
		temp_input_target.zp=temp_zp
		temp_input_target.fwhm=temp_fwhm
		temp_input_target.exptime=temp_exptime
		temp_input_target.n_chip=temp_n_chip
		temp_input_target.type=temp_type
		temp_type=( string(rebin(transpose(temp_type EQ 'short'),[(size(byte('_'+temp_type),/dim))[0],n_elements(temp_type)])*byte('_'+temp_type)) )
		temp_input_target.im_file=output_im_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.fits'
		temp_input_target.weight_file=output_im_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.WEIGHT.fits'
		temp_input_target.mask_file=output_im_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.MASK.fits'
		temp_input_target.sex_cat_file=output_sex_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.ldac'
		temp_input_target.sex_xml_file=output_sex_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.xml'
		temp_input_target.sex_check_file=output_sex_check_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.BACKGROUND.fits'
		temp_input_target.sex_zp_file=output_sex_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'_zp.dat'
		temp_input_target.scamp_cat_file=output_scamp_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'_stars.ldac'
		temp_input_target.scamp_head_file=output_scamp_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'_stars.head'
		temp_input_target.scamp_ahead_file=output_scamp_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'_stars.ahead'
		temp_input_target.swarp_im_file=output_im_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.fits'
		temp_input_target.swarp_head_file=output_im_dir[i]+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.head'

		temp_input_target.masked_im_file=output_im_dir[i]+'/'+'masked_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.fits'
		temp_input_target.masked_median_im_file=output_im_dir[i]+'/'+'masked_median_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.fits'
		temp_input_target.masked_res_im_file=output_im_dir[i]+'/'+'masked_res_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.fits'
		temp_input_target.sky_im_file=output_im_dir[i]+'/'+'sky_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.fits'
		temp_input_target.ss_im_file=output_im_dir[i]+'/'+'ss_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.fits'
		temp_input_target.ss_weight_file=output_im_dir[i]+'/'+'ss_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.WEIGHT.fits'
		temp_input_target.ss_sex_cat_file=output_sex_dir[i]+'/'+'ss_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.ldac'
		temp_input_target.ss_sex_xml_file=output_sex_dir[i]+'/'+'ss_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.xml'
		temp_input_target.ss_scamp_head_file=output_scamp_dir[i]+'/'+'ss_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'_stars.head'
		temp_input_target.ss_scamp_ahead_file=output_scamp_dir[i]+'/'+'ss_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'_stars.ahead'
		temp_input_target.ss_scamp_cat_file=output_scamp_dir[i]+'/'+'ss_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'_stars.ldac'
		temp_input_target.ss_swarp_head_file=output_im_dir[i]+'/'+'ss_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.head'
		temp_input_target.ss_swarp_im_file=output_im_dir[i]+'/'+'ss_fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+temp_type+'.fits'

		if i EQ 0 then input_target=temp_input_target $
		else input_target=[input_target,temp_input_target]

	endif
endfor
	
if n_elements(input_target) GT 0 then begin
	
	do_tile_full=strsplit(do_tile,',', /extract)
	do_filter_full=strsplit(do_filter,',', /extract)
	do_type_full=strsplit(do_type,',', /extract)

	if n_elements(do_tile_full) GT 1 then begin
		if n_elements(do_filter_full) GT 1 then begin
			temp_tile1=byte(input_target.tile)
			temp_tile2=byte(do_tile_full)
			temp_filter1=byte(input_target.filter)
			temp_filter2=byte(do_filter_full)

			gv=where( total( strcmp( string(rebin(temp_tile1,[(size(temp_tile1,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])), string(rebin(reform(temp_tile2,[(size(temp_tile2,/dim))[0],1,(size(temp_tile2,/dim))[1]]),[(size(temp_tile2,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])) ), 2) AND total( strcmp( string(rebin(temp_filter1,[(size(temp_filter1,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])), string(rebin(reform(temp_filter2,[(size(temp_filter2,/dim))[0],1,(size(temp_filter2,/dim))[1]]),[(size(temp_filter2,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])) ), 2) GT 0, n_gv)
			
		endif $
		else begin
 			temp_tile1=byte(input_target.tile)
      temp_tile2=byte(do_tile_full)

			if do_filter EQ 'all' then begin
				gv=where( total( strcmp( string(rebin(temp_tile1,[(size(temp_tile1,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])), string(rebin(reform(temp_tile2,[(size(temp_tile2,/dim))[0],1,(size(temp_tile2,/dim))[1]]),[(size(temp_tile2,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])) ), 2) GT 0, n_gv)
			endif $
			else begin
				gv=where( (total( strcmp( string(rebin(temp_tile1,[(size(temp_tile1,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])), string(rebin(reform(temp_tile2,[(size(temp_tile2,/dim))[0],1,(size(temp_tile2,/dim))[1]]),[(size(temp_tile2,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])) ), 2) GT 0) AND input_target.filter EQ do_filter, n_gv)
			endelse
	
		endelse

	endif $
	else begin

		if n_elements(do_filter_full) GT 1 then begin
			temp_filter1=byte(input_target.filter)
			temp_filter2=byte(do_filter_full)

			if do_tile EQ 'all' then begin
				gv=where( total( strcmp( string(rebin(temp_filter1,[(size(temp_filter1,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])), string(rebin(reform(temp_filter2,[(size(temp_filter2,/dim))[0],1,(size(temp_filter2,/dim))[1]]),[(size(temp_filter2,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])) ), 2) GT 0, n_gv)
			endif $
			else begin
				gv=where( input_target.tile EQ do_tile AND (total( strcmp( string(rebin(temp_filter1,[(size(temp_filter1,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])), string(rebin(reform(temp_filter2,[(size(temp_filter2,/dim))[0],1,(size(temp_filter2,/dim))[1]]),[(size(temp_filter2,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])) ), 2) GT 0), n_gv)
			endelse
		
		endif $
		else begin

			if do_tile EQ 'all' AND do_filter EQ 'all' then begin
				n_gv=n_elements(input_target)
				gv=indgen(n_elements(input_target))
			endif else $
			if do_tile EQ 'all' then begin
				gv=where(input_target.filter EQ do_filter, n_gv)
			endif else $
			if do_filter EQ 'all' then begin
				gv=where(input_target.tile EQ do_tile, n_gv)
			endif $
			else begin
				gv=where(input_target.tile EQ do_tile AND input_target.filter EQ do_filter, n_gv)
			endelse

		endelse
	endelse

	if n_gv GT 0 then input_target=input_target[gv] $
	else stop

	if do_type EQ 'all' then begin
		n_gv=n_elements(input_target)
		gv=indgen(n_elements(input_target))
	endif $
	else begin
		gv=where(input_target.type EQ do_type, n_gv)
	endelse

	if n_gv GT 0 then input_target=input_target[gv] $
	else stop

	forprint, input_target.im_file, input_target.tile, input_target.filter, FORMAT='A,2X,A,2X,A', textout=2 
endif

for i=0L, n_elements(do_program)-1 do begin

	if file_test('survey_iq_'+do_program[i]+'.dat') AND recipe NE 'iq' then begin
		readcol, 'survey_iq_'+do_program[i]+'.dat', temp_im_file, temp_filter, temp_tile, temp_dither, temp_sky, temp_fwhm, temp_type, FORMAT='A,A,A,A,F,F,A', COMMENT='#'
		n_survey=n_elements(temp_im_file)
		
		if n_survey GT 0 then begin
			create_struct, temp_input_iq, '', ['im_file','filter','tile','dither','sky','fwhm','type'], 'A,A,A,A,F,F,A', dim=n_survey
			temp_input_iq.im_file=temp_im_file
			temp_input_iq.filter=temp_filter
			temp_input_iq.tile=temp_tile
			temp_input_iq.dither=temp_dither
			temp_input_iq.sky=temp_sky
			temp_input_iq.fwhm=temp_fwhm
			temp_input_iq.type=temp_type

			if n_elements(input_iq) EQ 0 then input_iq=temp_input_iq $
			else input_iq=[input_iq,temp_input_iq]
		endif
	endif
endfor

if n_elements(input_iq) GT 0 AND recipe NE 'database' AND recipe NE 'setup' AND recipe NE 'sextractor' AND recipe NE 'iq' then begin
	print, 'Removing images with LOW IQ'

	gv_iq=where( (input_iq.type EQ 'long' AND input_iq.sky LT 12000.) OR (input_iq.type EQ 'short' AND input_iq.sky LT 400.), n_gv_iq, COMPLEMENT=bv_iq)
	print, 'Removing the following files'
	forprint, input_iq[bv_iq].im_file, input_iq[bv_iq].sky, FORMAT='A,4X,F0.1', text=2
	
	temp_file1=byte(input_target.im_file)
	temp_file2=byte(input_iq[gv_iq].im_file)
	
	gv=where( total( strcmp( string(rebin(temp_file1,[(size(temp_file1,/dim))[0],n_elements(input_target.im_file),n_gv_iq])), string(rebin(reform(temp_file2,[(size(temp_file2,/dim))[0],1,(size(temp_file2,/dim))[1]]),[(size(temp_file2,/dim))[0],n_elements(input_target.im_file),n_gv_iq])) ), 2) GT 0, n_gv)
	input_target=input_target[gv]

endif

for i=0L, n_elements(do_program)-1 do begin

	if file_test('survey_calib_'+do_program[i]+'.dat') AND recipe NE 'database' then begin
		readcol, 'survey_calib_'+do_program[i]+'.dat', temp_im_orig_file, temp_object, temp_filter, temp_weight_orig_file, temp_zp, temp_fwhm, temp_mjd, temp_exptime, temp_n_chip, temp_photflag, FORMAT='A,A,A,A,F,F,D,F,I,A', COMMENT='#'
		n_survey=n_elements(temp_im_orig_file)
		create_struct, temp_input_calib, '', ['im_orig_file','object','filter','tile','weight_orig_file','mjd','mjd_floor','zp','airmass','im_file','weight_file','fwhm','exptime','n_chip','sex_cat_file','sex_xml_file','sex_check_file','sex_zp_file','sex_zp_full_file','photflag','scamp_cat_file','scamp_head_file','scamp_ahead_file','swarp_im_file','swarp_head_file'], 'A,A,A,A,A,D,D,F,F,A,A,F,F,I,A,A,A,A,A,A,A,A,A,A,A', dim=n_survey
		temp_input_calib.im_orig_file=input_calib_dir[i]+'/'+repstr(temp_im_orig_file,'.fits.fz','.fits')
		temp_input_calib.weight_orig_file=input_calib_dir[i]+'/'+repstr(temp_weight_orig_file,'.fits.fz','.fits')
		temp_input_calib.object=temp_object
		temp_input_calib.filter=temp_filter
		temp_input_calib.mjd=temp_mjd
		temp_input_calib.mjd_floor=floor(temp_mjd+0.2)
		temp_input_calib.zp=temp_zp
		temp_input_calib.airmass=0.
		temp_input_calib.fwhm=temp_fwhm
		temp_input_calib.exptime=temp_exptime
		temp_input_calib.n_chip=temp_n_chip
		temp_input_calib.photflag=temp_photflag
		temp_mjd=string(temp_mjd, FORMAT='(F0.8)')
		temp_input_calib.im_file=output_im_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.fits'
		temp_input_calib.weight_file=output_im_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.WEIGHT.fits'
		temp_input_calib.sex_cat_file=output_sex_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.ldac'
		temp_input_calib.sex_xml_file=output_sex_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.xml'
		temp_input_calib.sex_check_file=output_sex_check_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.BACKGROUND.fits'
		temp_input_calib.sex_zp_file=output_calib_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_zp.dat'
		temp_input_calib.sex_zp_full_file=output_calib_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_zp_full.dat'
		temp_input_calib.scamp_cat_file=output_scamp_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_stars.ldac'
		temp_input_calib.scamp_head_file=output_scamp_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_stars.head'
		temp_input_calib.scamp_ahead_file=output_scamp_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_stars.ahead'
		temp_input_calib.swarp_im_file=output_im_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.fits'
		temp_input_calib.swarp_head_file=output_im_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.head'

		if i EQ 0 then input_calib=temp_input_calib $
		else input_calib=[input_calib,temp_input_calib]

	endif
endfor

if n_elements(input_calib) GT 0 then begin

	do_filter_split=strsplit(do_filter,',', /extract)
	if n_elements(do_filter_split) GT 1 then begin
		gv=where( total(rebin(byte(input_calib.filter),[n_elements(do_filter_full),n_elements(input_calib.filter)]) EQ rebin(transpose(byte(do_filter_split)),[n_elements(do_filter_full),n_elements(input_calib.filter)]),1) GT 0, n_gv)
		if n_gv GT 0 then input_calib=input_calib[gv] $
		else stop
	endif $
	else begin
		case do_filter of
			'all': begin
				n_gv=n_elements(input_calib)
				gv=indgen(n_elements(input_calib))
			end
			else: gv=where(input_calib.filter EQ do_filter, n_gv)
		endcase
		if n_gv GT 0 then input_calib=input_calib[gv] $
		else stop
	endelse

endif

for i=0L, n_elements(do_program)-1 do begin

	if file_test('survey_zp_'+do_program[i]+'.dat') then begin
		readcol, 'survey_zp_'+do_program[i]+'.dat', temp_mjd, temp_filter, temp_zp, temp_k, temp_zp_error, temp_k_error, temp_photflag, FORMAT='F,A,F,F,F,F,A', COMMENT='#'
		n_survey=n_elements(temp_mjd)
		create_struct, temp_input_zp, '', ['mjd','filter','zp','k','zp_error','k_error','photflag'], 'F,A,F,F,F,F,A', dim=n_survey
		temp_input_zp.mjd=temp_mjd
		temp_input_zp.filter=temp_filter
		temp_input_zp.zp=temp_zp
		temp_input_zp.k=temp_k
		temp_input_zp.zp_error=temp_zp_error
		temp_input_zp.k_error=temp_k_error
		temp_input_zp.photflag=temp_photflag

		if i EQ 0 then input_zp=temp_input_zp $
		else input_zp=[input_zp,temp_input_zp]

	endif
endfor

wim_file=file_search(output_stack_swarp_dir+'/*.WEIGHT.fits')
im_file=repstr(wim_file, '.WEIGHT.fits', '.fits')
n_gv_im=n_elements(im_file)
;gv_wim=where(stregex(im_file, '.*WEIGHT.fits', /boolean) EQ 1, n_gv_wim);, COMPLEMENT=gv_im, NCOMPLEMENT=n_gv_im)
;gv_im=where(stregex(im_file, '.*(WEIGHT.fits|MASK.fits)', /boolean) EQ 0, n_gv_im);, COMPLEMENT=gv_im, NCOMPLEMENT=n_gv_im)

create_struct, input_stack, '', ['im_file','filter','tile','type','weight_file'], 'A,A,A,A,A', dim=n_gv_im

input_stack.im_file=im_file ;[gv_im]
input_stack.weight_file=wim_file ;[gv_wim]
temp_im_file=strsplit(im_file,'/', /extract)
for i=0L, n_gv_im-1 do begin
	temp=strsplit( (temp_im_file[i])[-1], '_', /extract)
	input_stack[i].tile=repstr(temp[1],'tile','')
	input_stack[i].filter=temp[2]
	input_stack[i].type=repstr(temp[3],'.fits','')
endfor

do_program_plain=repstr(do_program_orig,',','_')
do_filter_plain=repstr(do_filter_orig,',','_')
loadct, 12


if recipe EQ 'database' then begin

	for ii=0L, n_elements(do_program)-1 do begin

		im_file=file_search(input_im_dir[ii], '*.fz')
		n_im_file=n_elements(im_file)
	
		create_struct, temp_target, '', ['dir','im_file','weight_file','expnum','type','date','ra','dec','mjd','mjd_floor','exptime','tile','dither','filter','zp','fwhm','n_chip','type_exptime'], 'A,A,A,I,A,A,D,D,D,D,F,A,A,A,F,F,I,A', dim=1
	
		for i=0L, n_im_file-1 do begin
			print, 'DATABASE - Processing file ', (strsplit(im_file[i],'/',/extract))[-1]
	
			im_h=headfits(im_file[i])
			fits_info, im_file[i], n_ext=n_ext, /silent
			temp_target.dir=input_im_dir[ii]
			temp_target.im_file=(strsplit(im_file[i],'/',/extract))[-1]
			temp_target.expnum=fxpar(im_h, 'EXPNUM')
			temp_target.type=strtrim(fxpar(im_h, 'PRODTYPE'),2)
			temp_target.date=fxpar(im_h, 'DATE-OBS')
			temp_target.mjd=fxpar(im_h, 'MJD-OBS')
			temp_target.mjd_floor=floor(fxpar(im_h, 'MJD-OBS')+0.2)
			temp_target.ra=ten(fxpar(im_h, 'RA'))*360./24.
			temp_target.dec=ten(fxpar(im_h, 'DEC'))
			temp_target.exptime=float(fxpar(im_h, 'EXPTIME'))
			temp_target.filter=strmid(fxpar(im_h, 'FILTER'),0,1)
			temp=size(fxpar(im_h, 'MAGZERO'),/type)
			temp_target.zp= (temp EQ 4 OR temp EQ 5) ? fxpar(im_h, 'MAGZERO') : 0.
			temp_target.n_chip=n_ext
	
			if i EQ 0 then input_target=temp_target else input_target=[input_target,temp_target]
		endfor
		gv_sort=sort(input_target.date+' '+input_target.type)
		input_target=input_target[gv_sort]
	
		if file_test(input_im_dir[ii]+'/pipeline_database.dat') EQ 0 then begin
			openw, lun, input_im_dir[ii]+'/pipeline_database.dat', /get_lun
			printf, lun, '# FILE	 DATE-OBS	 EXPNUM	 OBJECT	 FILTER	 EXPTIME	 SEQID	 OBSTYPE	 PROCTYPE		PRODTYPE	 MAGZERO'
			for i=0L, n_im_file-1 do begin
				command='dfits '+im_file[gv_sort[i]]+' | fitsort -d DATE-OBS EXPNUM OBJECT FILTER EXPTIME SEQID OBSTYPE PROCTYPE PRODTYPE MAGZERO'
				spawn, command, result
				result=repstr(result, input_im_dir[ii]+'/', '')
				printf, lun, result
			endfor
			free_lun, lun
		endif
	
		gv_im=where(input_target.type EQ 'image' AND input_target.exptime GT 30., n_gv_im)
		gv_weight=where(input_target.type EQ 'wtmap' AND input_target.exptime GT 30., n_gv_weight)
		gv_sort=sort(input_target[gv_im].date+' '+input_target[gv_im].type)
		mjd_uniq=input_target[gv_im[uniq(input_target[gv_im].mjd_floor, sort(input_target[gv_im].mjd_floor))]].mjd_floor
		filter_uniq=input_target[gv_im[uniq(input_target[gv_im].filter, sort(input_target[gv_im].filter))]].filter
	
		input_target[gv_im].type_exptime='long'
		for i=0L, n_gv_im-1 do begin
			gv=where(input_target[gv_weight].expnum EQ input_target[gv_im[i]].expnum, n_gv)
			if n_gv EQ 1 then input_target[gv_im[i]].weight_file=input_target[gv_weight[gv]].im_file $
			else stop
		endfor
	
;		input_target=input_target[gv_im[gv_sort]]
	
		for i=0L, n_elements(survey_info)-1 do begin
			for j=0L, n_elements((survey_info[i]).filter)-1 do begin
				gv=gv_im[where( sqrt( (input_target[gv_im].ra-ten((survey_info[i]).coo[0])*360./24)^2 + (input_target[gv_im].dec-ten((survey_info[i]).coo[1]))^2 ) LE 0.1 AND input_target[gv_im].filter EQ (survey_info[i]).filter[j], n_gv)]
				if n_gv GT 0 then begin
					input_target[gv].tile=(survey_info[i]).tile
					print, 'Found images for filter and tile: ', (survey_info[i]).filter[j], ' ', (survey_info[i]).tile
				endif
	;			input_target[gv_im[gv_sort[gv]]].dither=string(indgen(n_gv)+1,FORMAT='(I0)')
			endfor
		endfor
    gv_im=gv_im[where(input_target[gv_im].tile NE ' ', n_gv_im)]
    tile_uniq=input_target[gv_im[uniq(input_target[gv_im].tile, sort(input_target[gv_im].tile))]].tile
	
		for i=0L, n_elements(filter_uniq)-1 do begin
			i_dither=1L
			for j=0L, n_elements(mjd_uniq)-1 do begin
				gv=gv_im[where( input_target[gv_im].filter EQ filter_uniq[i] AND input_target[gv_im].mjd_floor EQ mjd_uniq[j], n_gv)]
				if n_gv GT 0 then begin
					for k=0L, n_tags(survey_sequence)-1 do begin
						gv_dither=gv[where(input_target[gv].dither EQ ' ', n_gv_dither)]
						if n_gv_dither GT 0 then begin
							temp_n_seq=n_elements(survey_sequence.(k))
							temp_tile=[input_target[gv_dither].tile,strarr(temp_n_seq-1)]
							temp_match=intarr(n_gv_dither)
							for l=0L, n_gv_dither-1 do begin
								temp_data=(shift(temp_tile,-l))[0:temp_n_seq-1]
								temp_match[l]=total(temp_data EQ survey_sequence.(k))
							endfor
							gv_match=where(temp_match GE n_elements(survey_sequence.(k))/2, n_gv_match)
							gv_match=[gv_dither[gv_match[where( gv_match - [-n_elements(survey_sequence.(k))/2,gv_match] GE n_elements(survey_sequence.(k))/2, n_gv_match)]],gv_dither[-1]+1]
							for l=0L, n_gv_match-1 do begin
								m=gv_match[l]
								temp_tile=''
								while m LT gv_match[l+1] do begin
									if input_target[m].filter EQ filter_uniq[i] AND input_target[m].mjd_floor EQ mjd_uniq[j] then begin
										gv_tile=where(input_target[m].tile EQ temp_tile, n_gv_tile)
										if n_gv_tile EQ 0 then input_target[m].dither=strn(i_dither) $
										else input_target[m].dither=strn(i_dither)+'b'
										temp_tile=[temp_tile,input_target[m].tile]
									endif
									m++
								endwhile
								i_dither++
							endfor
						endif
					endfor
				endif
				
			endfor
		endfor
		forprint, input_target[gv_im].im_file, input_target[gv_im].mjd_floor, input_target[gv_im].filter, input_target[gv_im].tile, input_target[gv_im].dither, format='A,4X,F0,4X,A,4X,A,4X,A', textout=2
	
;		gv=where(input_target.tile NE ' ' AND input_target.dither NE ' ', n_gv)
		if n_gv_im GT 0 then begin
			openw, lun, 'survey_target_'+do_program[ii]+'.dat', /get_lun
			printf, lun, '#  im_file        filter  tile  dither  weight_file     zp      FWHM    MJD   EXPTIME   NCHIP   TYPE_EXPTIME'
			for i=0L, n_gv_im-1 do begin
				printf, lun, input_target[gv_im[i]].im_file, input_target[gv_im[i]].filter, input_target[gv_im[i]].tile, input_target[gv_im[i]].dither, input_target[gv_im[i]].weight_file, input_target[gv_im[i]].zp, input_target[gv_im[i]].fwhm, input_target[gv_im[i]].mjd, input_target[gv_im[i]].exptime, input_target[gv_im[i]].n_chip, input_target[gv_im[i]].type_exptime, FORMAT='(A,4X,A,4X,A,4X,A,4X,A,4X,F5.2,4X,F4.1,4X,F0.6,4X,F0.1,4X,I0,4X,A)'
			endfor
			free_lun, lun
		endif

; Short exposures
		gv_im=where(input_target.type EQ 'image' AND input_target.exptime LE 30., n_gv_im)
		gv_weight=where(input_target.type EQ 'wtmap' AND input_target.exptime LE 30., n_gv_weight)
		gv_sort=sort(input_target[gv_im].date+' '+input_target[gv_im].type)
		mjd_uniq=input_target[gv_im[uniq(input_target[gv_im].mjd_floor, sort(input_target[gv_im].mjd_floor))]].mjd_floor
		filter_uniq=input_target[gv_im[uniq(input_target[gv_im].filter, sort(input_target[gv_im].filter))]].filter
	
		input_target[gv_im].type_exptime='short'
		for i=0L, n_gv_im-1 do begin
			gv=where(input_target[gv_weight].expnum EQ input_target[gv_im[i]].expnum, n_gv)
			if n_gv EQ 1 then input_target[gv_im[i]].weight_file=input_target[gv_weight[gv]].im_file $
			else stop
		endfor
	
		for i=0L, n_elements(survey_info)-1 do begin
			for j=0L, n_elements((survey_info[i]).filter)-1 do begin
				gv=gv_im[where( sqrt( (input_target[gv_im].ra-ten((survey_info[i]).coo[0])*360./24)^2 + (input_target[gv_im].dec-ten((survey_info[i]).coo[1]))^2 ) LE 0.1 AND input_target[gv_im].filter EQ (survey_info[i]).filter[j], n_gv)]
				if n_gv GT 0 then begin
					input_target[gv].tile=(survey_info[i]).tile
					print, 'Found images for filter and tile: ', (survey_info[i]).filter[j], ' ', (survey_info[i]).tile
				endif
			endfor
		endfor
    gv_im=gv_im[where(input_target[gv_im].tile NE ' ', n_gv_im)]
    tile_uniq=input_target[gv_im[uniq(input_target[gv_im].tile, sort(input_target[gv_im].tile))]].tile

    for i=0L, n_elements(filter_uniq)-1 do begin
      i_dither=1L
      for j=0L, n_elements(tile_uniq)-1 do begin
        gv=where( input_target[gv_im].filter EQ filter_uniq[i] AND input_target[gv_im].tile EQ tile_uniq[j], n_gv)
        if n_gv GT 0 then begin
          input_target[gv_im[gv]].dither=string(indgen(n_gv)+1, FORMAT='(I0)')
        endif
      endfor
    endfor
    forprint, input_target[gv_im].im_file, input_target[gv_im].mjd_floor, input_target[gv_im].filter, input_target[gv_im].tile, input_target[gv_im].dither, format='A,4X,F0,4X,A,4X,A,4X,A', textout=2
	
		if n_gv_im GT 0 then begin
			openw, lun, 'survey_target_'+do_program[ii]+'.dat', /get_lun, /append
			for i=0L, n_gv_im-1 do begin
				printf, lun, input_target[gv_im[i]].im_file, input_target[gv_im[i]].filter, input_target[gv_im[i]].tile, input_target[gv_im[i]].dither, input_target[gv_im[i]].weight_file, input_target[gv_im[i]].zp, input_target[gv_im[i]].fwhm, input_target[gv_im[i]].mjd, input_target[gv_im[i]].exptime, input_target[gv_im[i]].n_chip, input_target[gv_im[i]].type_exptime, FORMAT='(A,4X,A,4X,A,4X,A,4X,A,4X,F5.2,4X,F4.1,4X,F0.6,4X,F0.1,4X,I0,4X,A)'
			endfor
			free_lun, lun
		endif


		im_file=file_search(input_calib_dir[ii], '*.fz')
		n_im_file=n_elements(im_file)
		im_date=list()
		im_type=list()
	
		for i=0L, n_im_file-1 do begin
			im_h=headfits(im_file[i])
			im_date.add, fxpar(im_h, 'DATE-OBS')
			im_type.add, fxpar(im_h, 'PRODTYPE')
		endfor
		im_date=im_date.toarray(type='string')
		im_type=im_type.toarray(type='string')
		gv=sort(im_date+' '+im_type)
	
		openw, lun, input_calib_dir[ii]+'/pipeline_database.dat', /get_lun
		printf, lun, '# FILE	 DATE-OBS   MJD-OBS	 EXPNUM	 OBJECT	 FILTER	 EXPTIME	 SEQID	 OBSTYPE	 PROCTYPE		PRODTYPE	 MAGZERO'
		for i=0L, n_im_file-1 do begin
			command='dfits '+im_file[gv[i]]+' | fitsort -d DATE-OBS MJD-OBS EXPNUM OBJECT FILTER EXPTIME SEQID OBSTYPE PROCTYPE PRODTYPE MAGZERO'
			spawn, command, result
			result=repstr(result, input_calib_dir[ii]+'/', '')
			printf, lun, result
		endfor
		free_lun, lun

	endfor

endif else $
if recipe EQ 'setup' then begin

	for i=0L, n_elements(do_program)-1 do begin
		if file_test(output_im_dir[i], /directory) EQ 0 then file_mkdir, output_im_dir[i]
		if file_test(output_calib_dir[i], /directory) EQ 0 then file_mkdir, output_calib_dir[i]
		if file_test(output_sex_dir[i], /directory) EQ 0 then file_mkdir, output_sex_dir[i]
		if file_test(output_sex_check_dir[i], /directory) EQ 0 then file_mkdir, output_sex_check_dir[i]
		if file_test(output_scamp_dir[i], /directory) EQ 0 then file_mkdir, output_scamp_dir[i]
	endfor
	if file_test(output_stack_swarp_dir, /directory) EQ 0 then file_mkdir, output_stack_swarp_dir
	if file_test(output_stack_sex_dir, /directory) EQ 0 then file_mkdir, output_stack_sex_dir
	if file_test(output_stack_check_dir, /directory) EQ 0 then file_mkdir, output_stack_check_dir
	if file_test(output_stack_psfex_dir, /directory) EQ 0 then file_mkdir, output_stack_psfex_dir

	if do_overwrite then begin
		file_delete, input_target.im_file, /noexpand, /allow_non, /quiet
		file_delete, input_target.weight_file, /noexpand, /allow_non, /quiet
	endif

	for i=0, n_elements(input_target)-1 do begin

		print, 'SETUP - Processing file ', input_target[i].im_orig_file

		if file_test(input_target[i].im_orig_file+'.fz', /regular) EQ 0 OR file_test(input_target[i].weight_orig_file+'.fz', /regular) EQ 0 then begin
			print, 'Some files are missing'
			stop
		endif

		if file_test(input_target[i].im_orig_file, /regular) EQ 0 then begin
			command='funpack '+input_target[i].im_orig_file+'.fz'
			print, command
			spawn, command
		endif

		if file_test(input_target[i].weight_orig_file, /regular) EQ 0 then begin
			command='funpack '+input_target[i].weight_orig_file+'.fz'
			print, command
			spawn, command
		endif

		if file_test(input_target[i].im_file, /regular) EQ 0 then begin
			command='ln -s '+input_target[i].im_orig_file+' '+input_target[i].im_file
			print, command
			spawn, command
		endif

		if file_test(input_target[i].weight_file, /regular) EQ 0 then begin
			command='ln -s '+input_target[i].weight_orig_file+' '+input_target[i].weight_file
			print, command
			spawn, command
		endif
	endfor

	for i=0, n_elements(input_calib)-1 do begin

		if file_test(input_calib[i].im_orig_file+'.fz', /regular) EQ 0 OR file_test(input_calib[i].weight_orig_file+'.fz', /regular) EQ 0 then begin
			print, 'Some files are missing'
			stop
		endif

		if file_test(input_calib[i].im_orig_file, /regular) EQ 0 then begin
			command='funpack '+input_calib[i].im_orig_file+'.fz'
			print, command
			spawn, command
		endif

		if file_test(input_calib[i].weight_orig_file, /regular) EQ 0 then begin
			command='funpack '+input_calib[i].weight_orig_file+'.fz'
			print, command
			spawn, command
		endif

		if file_test(input_calib[i].im_file, /regular) EQ 0 then begin
			command='ln -s '+input_calib[i].im_orig_file+' '+input_calib[i].im_file
			print, command
			spawn, command
		endif

		if file_test(input_calib[i].weight_file, /regular) EQ 0 then begin
			command='ln -s '+input_calib[i].weight_orig_file+' '+input_calib[i].weight_file
			print, command
			spawn, command
		endif
	endfor

endif else $
if recipe EQ 'compute zp' then begin

	vig_diam=101
	sex_mag_bin=0.5
	sex_radius_bin=0.2

	if do_overwrite then file_delete, input_calib.sex_cat_file, /noexpand, /allow_non, /quiet

	for i=0L, n_elements(input_calib)-1 do begin

		im_h=headfits(input_calib[i].im_file)
		im_filter=fxpar(im_h, 'FILTER')
		im_exptime=fxpar(im_h, 'EXPTIME')
		im_zp=fxpar(im_h, 'MAGZERO')
		im_ra=ten(fxpar(im_h, 'RA'))*360./24.
		im_dec=ten(fxpar(im_h, 'DEC'))
		im_mjd=fxpar(im_h, 'MJD-OBS')
		im_airmass=tai2airmass(im_ra,im_dec,2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)
		input_calib[i].airmass=im_airmass

		print, input_calib[i].im_file, date_conv(im_mjd+2400000.5D,'S'), im_zp, im_exptime, im_airmass, im_filter, FORMAT='("Image filename ",A,2X,A,2X,F0.2,2X,F0.2,2X,F0.2,2X,A)'

		case input_calib[i].filter of
			'u': begin
				sex_satur_level='30000.'
				sex_mag_range=[7,15.]
				sex_flux_radius_min=1.4
				plot_mag_range=[18,5]
				plot_flux_radius_range=[2,8]
				end
			'g': begin
				sex_satur_level='30000.'
				sex_mag_range=[11,19.]
				sex_flux_radius_min=1.4
				plot_mag_range=[20,8]
				plot_flux_radius_range=[2,8]
				end
			'i': begin
				sex_satur_level='30000.'
				sex_mag_range=[12,18.]
				sex_flux_radius_min=1.4
				plot_mag_range=[20,8]
				plot_flux_radius_range=[2,8]
				end
			'z': begin
				sex_satur_level='30000.'
				sex_mag_range=[12,20.]
				sex_flux_radius_min=1.4
				plot_mag_range=[20,8]
				plot_flux_radius_range=[2,8]
				end
			else: stop
		endcase

		im_h=headfits(input_calib[i].im_file)
		im_ra=ten(fxpar(im_h, 'RA'))*360./24.
		im_dec=ten(fxpar(im_h, 'DEC'))
		im_dec_sign=im_dec GT 0. ? '+' : '-'
		im_exptime=fxpar(im_h, 'EXPTIME')
		im_zp=fxpar(im_h, 'MAGZERO')
		im_mjd=fxpar(im_h, 'MJD-OBS')
		im_airmass=tai2airmass(im_ra,im_dec,2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)
		input_calib[i].airmass=im_airmass

		file_sex_cat=file_info(input_calib[i].sex_cat_file, /noexpand)
		file_sex_zp=file_info(input_calib[i].sex_zp_file, /noexpand)
		file_sex_zp_full=file_info(input_calib[i].sex_zp_full_file, /noexpand)

		if file_sex_cat.exists AND file_sex_cat.size GT 0 then fits_info, input_calib[i].sex_cat_file, n_ext=sex_n_ext, /silent $
		else sex_n_ext=0

		if sex_n_ext EQ input_calib[i].n_chip*2 AND file_sex_zp.exists AND file_sex_zp.size GT 100 AND file_sex_zp_full.exists AND file_sex_zp_full.size GT 100 then continue

		im_fwhm=list()
		openw, lun_zp, input_calib[i].sex_zp_file, /get_lun
		openw, lun_zp_full, input_calib[i].sex_zp_full_file, /get_lun
		printf, lun_zp, '# Ext  DECam_zp   airmass   zp   zp_err  zp_nstars   FWHM'
		printf, lun_zp_full, '# Ext  NUMBER   RA   DEC   FLUX_APER   FLUX_AUTO   mag_ref   mag_diff   airmass'

		if file_sex_cat.exists AND file_sex_cat.size GT 0 then fits_info, input_calib[i].sex_cat_file, n_ext=sex_n_ext, /silent $
		else sex_n_ext=0

		if sex_n_ext NE input_calib[i].n_chip*2 OR do_overwrite EQ 1 then begin

			for j=0L, n_elements(input_chip_fwhm)-1 do begin
				im_data=readfits(input_calib[i].im_file, im_h, ext=input_chip_fwhm[j])
				writefits, output_sex_dir+'/im_sextractor.fits', im_data, im_h
				wim_data=readfits(input_calib[i].weight_file, wim_h, ext=input_chip_fwhm[j])
				writefits, output_sex_dir+'/wim_sextractor.fits', wim_data, wim_h
		
				im_size=size(im_data, /dim)
				im_gain=fxpar(im_h, 'GAINA')
				im_ron=fxpar(im_h, 'RDNOISEA')
	
				command = 'sex '+output_sex_dir+'/im_sextractor.fits' +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+output_sex_dir+'/im_sextractor.ldac'+' -WEIGHT_IMAGE '+output_sex_dir+'/wim_sextractor.fits'+' -XML_NAME '+input_calib[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(input_calib[i].zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_calib[i].sex_check_file
				print, command
				spawn, command
		
				cat_sex=mrdfits(output_sex_dir+'/im_sextractor.ldac', 2, cat_sex_h, COLUMNS=['NUMBER','X_IMAGE','Y_IMAGE','FLUX_RADIUS','MAG_AUTO','FLUX_AUTO','FLAGS'], /silent)
				plothist, (cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] and cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_flux_radius_min AND cat_sex.flags LT 8, n_gv)]).flux_radius, temp_xhist, temp_yhist, bin=sex_radius_bin, /noplot
				temp=max(temp_yhist, gv) & sex_radius=temp_xhist[gv]
				sex_radius=median((cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 , n_gv)]).flux_radius)

				gv_plot=where(cat_sex.flags LT 8, n_gv_plot)
				plot, cat_sex[gv_plot].flux_radius, cat_sex[gv_plot].mag_auto, psym=1, xrange=plot_flux_radius_range, yrange=plot_mag_range, /ystyle, /xstyle
				oplot, sex_radius*[1,1], [0,100], color=200
				oplot, sex_radius*[0.9,0.9], [0,100], line=2, color=200
				oplot, sex_radius*[1.1,1.1], [0,100], line=2, color=200
				oplot, [0,100], sex_mag_range[0]*[1,1], line=2, color=100
				oplot, [0,100], sex_mag_range[1]*[1,1], line=2, color=100
				wait, 0.2
	
				gv_stars=where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1]-1 AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 AND cat_sex.flags LE 1 AND cat_sex.x_image GT vig_diam/2. AND cat_sex.x_image LT (im_size[0]-vig_diam/2.) AND cat_sex.y_image GT vig_diam/2. AND cat_sex.y_image LT (im_size[1]-vig_diam/2.), n_gv_stars)
	
		    if n_gv_stars GT 0 then begin
					gv_stars=gv_stars[sort(cat_sex[gv_stars].mag_auto)]
	
					temp_n=n_gv_stars<5
	  	    temp_x=cat_sex[gv_stars].x_image
	   	  	temp_y=cat_sex[gv_stars].y_image
	
	  	    for k=0L, temp_n-1 do begin
		        x_range=[ floor(temp_x[k]-(vig_diam-1.)/2), ceil(temp_x[k]+(vig_diam-1.)/2) ]
	    	    y_range=[ floor(temp_y[k]-(vig_diam-1.)/2), ceil(temp_y[k]+(vig_diam-1.)/2) ]
	      	  vig_data = im_data[x_range[0]:x_range[1],y_range[0]:y_range[1]]; - im_sky
						vig_size=size(vig_data, /dim)
	
		        x_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
	    	    y_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
						vig_center_data=vig_data[x_center_range[0]:x_center_range[1],y_center_range[0]:y_center_range[1]]
						vig_center_size=size(vig_center_data, /dim)
	
		        im_max=max(vig_center_data, gv_max)
	        	im_c= [gv_max mod vig_center_size[0], gv_max/vig_center_size[0]]
						im_c += [x_center_range[0], y_center_range[0]]

						dist_circle, vig_mask, vig_size, im_c[0], im_c[1]
						vig_sky=median(vig_data[where(vig_mask GT 20., n_sky)])
						vig_mag=25.-2.5*alog10(total(vig_data[where(vig_mask LE 20., n_star)]) - n_star*vig_sky )

						vig_res=abs( vig_data-vig_sky - max(vig_center_data-vig_sky)/2. )
						gv=sort(vig_res*vig_mask^2)
						vig_fwhm=2*median(vig_mask[gv[1:5]])
						print, 'Radius for computing FWHM ', vig_mask[gv[1:5]]
						plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
						oplot, [0,100], max(vig_data-vig_sky)/2.*[1,1], line=2, color=200
						oplot, vig_fwhm/2.*[1,1], [-1e5,1e5], line=2, color=200

						vig_skyrad=4*vig_fwhm < 40.
						vig_psfrad=3*vig_fwhm < 50.
						vig_fitrad=vig_fwhm < 50.
	
						gcntrd, vig_data, im_c[0], im_c[1], im_cx, im_cy, vig_fwhm
						dist_circle, vig_mask, vig_size, im_cx, im_cy
						vig_sky=median(vig_data[where(vig_mask GT vig_skyrad, n_sky)])
						vig_mag=25.-2.5*alog10(total(vig_data[where(vig_mask LE vig_skyrad, n_star)]) - n_star*vig_sky )
						plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
						oplot, [0,100], max(vig_data-vig_sky)/2.*[1,1], line=2, color=200
						oplot, vig_fwhm/2.*[1,1], [-1e5,1e5], line=2, color=200

						getpsf, vig_data, im_cx, im_cy, vig_mag, vig_sky, im_ron, im_gain, psf_param, psf_residuals, [0], vig_psfrad, vig_fitrad, output_sex_dir+'/im_sextractor_psf.fits' 
		        im_fwhm.add, 2*sqrt(2*alog(2))*sqrt( (psf_param[3]^2 + psf_param[4]^2)/2. )
						print, 'FWHM ', im_fwhm[-1],' pixels'
	
						dist_circle, vig_mask, vig_size, im_cx, im_cy
						plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
						x=findgen(100)/10.
						y=gaussian(x,[psf_param[0], 0., mean(psf_param[3:4])])
						oplot, x, y, line=2, color=200
						oplot, im_fwhm[-1]/2.*[1,1], [-1e5,1e5], line=2, color=200
						
	  	  	endfor
				endif	
	
			endfor
	
			im_fwhm=median(im_fwhm.toarray(type='float'))
			print, 'Median FWHM ', im_fwhm,' pixels'

;		endif


;		im_h=headfits(input_calib[i].im_file)
;		im_ra=ten(fxpar(im_h, 'RA'))*360./24.
;		im_dec=ten(fxpar(im_h, 'DEC'))
;		im_dec_sign=im_dec GT 0. ? '+' : '-'
;		im_exptime=fxpar(im_h, 'EXPTIME')
;		im_zp=fxpar(im_h, 'MAGZERO')
;		im_mjd=fxpar(im_h, 'MJD-OBS')
;		im_airmass=tai2airmass(im_ra,im_dec,2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)
;		input_calib[i].airmass=im_airmass


;		if file_test(input_calib[i].sex_cat_file, /regular) then fits_info, input_calib[i].sex_cat_file, n_ext=sex_n_ext, /silent $
;		else sex_n_ext=0
;
;		if sex_n_ext NE input_n_chip*2 OR do_overwrite EQ 1 then begin

			command = 'sex '+input_calib[i].im_file +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+input_calib[i].sex_cat_file+' -WEIGHT_IMAGE '+input_calib[i].weight_file+' -XML_NAME '+input_calib[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(input_calib[i].zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_calib[i].sex_check_file+' -PHOT_APERTURES '+strjoin(string(2*[1,2,3,4,5]*im_fwhm, FORMAT='(F0.2)'),',')
			print, command
			spawn, command

		endif

;		pipe_airmass=fltarr(input_n_chip)
;		for j=0L, input_n_chip-1 do begin
;			im_h=headfits(input_calib[i].im_file, ext=j+1)
;			extast, im_h, im_ast
;			xy2ad, im_ast.naxis[0]/2., im_ast.naxis[1]/2., im_ast, pipe_ra, pipe_dec
;			pipe_airmass[j]=tai2airmass(pipe_ra,pipe_dec,2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)
;		endfor

		im_h=headfits(input_calib[i].im_file)
		im_ra=ten(fxpar(im_h, 'RA'))*360./24.
		im_dec=ten(fxpar(im_h, 'DEC'))
		im_dec_sign=im_dec GT 0. ? '+' : '-'

		if do_standard EQ 'southern' then begin

			cat_file=file_search('standard/southern_stars', '*.dat.trimmed')
			for j=0L, n_elements(cat_file)-1 do begin
				readcol, cat_file[j], id, ra, dec, mag_u, mag_u_err, mag_u_n, mag_g, mag_g_err, mag_g_n, mag_r, mag_r_err, mag_r_n, mag_i, mag_i_err, mag_i_n, mag_z, mag_z_err, mag_z_n, FORMAT='A,A,A,F,F,I,F,F,I,F,F,I,F,F,I,F,F,I,X,X', comment='#'
				create_struct, cat_ref_partial, '', ['ra','dec','mag','magerr'], 'F,F,5F,5F', dim=n_elements(id)

				cat_ref_partial.ra=tenv(ra)*360./24
				cat_ref_partial.dec=tenv(dec)
				cat_ref_partial.mag[0]=mag_u
				cat_ref_partial.mag[1]=mag_g
				cat_ref_partial.mag[2]=mag_r
				cat_ref_partial.mag[3]=mag_i
				cat_ref_partial.mag[4]=mag_z
				cat_ref_partial.magerr[0]=mag_u_err
				cat_ref_partial.magerr[1]=mag_g_err
				cat_ref_partial.magerr[2]=mag_r_err
				cat_ref_partial.magerr[3]=mag_i_err
				cat_ref_partial.magerr[4]=mag_z_err
			
				if j EQ 0 then cat_ref=cat_ref_partial $
				else cat_ref=[cat_ref,cat_ref_partial]
				
			endfor
			n_gv_ref=n_elements(cat_ref)
	
		endif else $
		if do_standard EQ 'sdss' then begin
			command='aclient_cgi cocat1.u-strasbg.fr sdss9 -c '+string(im_ra, im_dec_sign, abs(im_dec), FORMAT='(F0.6,A,F0.6)')+' -r 80. -lmg 10.,20. -m 10000000'
			spawn, command, cat_ref_data

			cat_ref_data=repstr(cat_ref_data,'---','NaN')
			gv_ref=where( strmid(cat_ref_data,0,1) NE '#', n_gv_ref)
			cat_ref_data=cat_ref_data[gv_ref]
			create_struct, cat_ref, '', ['ra','dec','mag','magerr'], 'F,F,5F,5F', dim=n_gv_ref
			cat_ref.ra=float(strmid(cat_ref_data, 26, 10))
			cat_ref.dec=float(strmid(cat_ref_data, 36, 10))
			cat_ref.mag[0]=float(strmid(cat_ref_data, 71, 6))
			cat_ref.mag[1]=float(strmid(cat_ref_data, 84, 6))
			cat_ref.mag[2]=float(strmid(cat_ref_data, 97, 6))
			cat_ref.mag[3]=float(strmid(cat_ref_data, 110, 6))
			cat_ref.mag[4]=float(strmid(cat_ref_data, 123, 6))
			cat_ref.magerr[0]=float(strmid(cat_ref_data, 78, 5))
			cat_ref.magerr[1]=float(strmid(cat_ref_data, 91, 5))
			cat_ref.magerr[2]=float(strmid(cat_ref_data, 104, 5))
			cat_ref.magerr[3]=float(strmid(cat_ref_data, 117, 5))
			cat_ref.magerr[4]=float(strmid(cat_ref_data, 130, 5))
		endif else $
			stop

		flux_radius_median=list()
		foreach j, input_chip_flux_radius do begin
			cat_sex=mrdfits(input_calib[i].sex_cat_file, 2*(j+1), cat_sex_h, COLUMNS=['NUMBER','ALPHA_J2000','DELTA_J2000','FLUX_RADIUS','MAG_AUTO','FLUX_AUTO','MAG_APER','FLUX_APER','FLAGS'], /silent)
			plothist, (cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] and cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_flux_radius_min AND cat_sex.flags LT 8, n_gv)]).flux_radius, temp_xhist, temp_yhist, bin=sex_radius_bin, /noplot
			temp=max(temp_yhist, gv) & sex_radius=temp_xhist[gv]
			sex_radius=median((cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 , n_gv)]).flux_radius)
			print, 'Flux radius for chip ', j, sex_radius
			flux_radius_median.add, sex_radius
		endforeach
		flux_radius_median=median(flux_radius_median.toarray(type='FLOAT'))

		
		for j=0L, input_calib[i].n_chip-1 do begin
			im_h=headfits(input_calib[i].im_file, ext=j+1)
			extast, im_h, im_ast
			xy2ad, im_ast.naxis[0]/2., im_ast.naxis[1]/2., im_ast, pipe_ra, pipe_dec
			pipe_airmass=tai2airmass(pipe_ra,pipe_dec,2000., mjd=input_calib[i].mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)

			cat_sex=mrdfits(input_calib[i].sex_cat_file, 2*(j+1), cat_sex_h, COLUMNS=['NUMBER','ALPHA_J2000','DELTA_J2000','FLUX_RADIUS','MAG_AUTO','FLUX_AUTO','MAG_APER','FLUX_APER','FLAGS'], /silent)
;			plothist, (cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] and cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_flux_radius_min AND cat_sex.flags LT 8, n_gv)]).flux_radius, temp_xhist, temp_yhist, bin=sex_radius_bin, /noplot
;			temp=max(temp_yhist, gv) & sex_radius=temp_xhist[gv]
			sex_radius=flux_radius_median
			sex_radius=median((cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 , n_gv)]).flux_radius)
			gv_sex=where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 AND cat_sex.flags LE 1, n_gv_sex)

			gv_plot=where(cat_sex.flags LT 8, n_gv_plot)
			plot, [cat_sex[gv_plot].flux_radius], [cat_sex[gv_plot].mag_auto], psym=1, xrange=plot_flux_radius_range, yrange=plot_mag_range, /ystyle, /xstyle
			oplot, [cat_sex[gv_sex].flux_radius], [cat_sex[gv_sex].mag_auto], psym=1, color=200
			oplot, sex_radius*[1,1], [0,100], color=200
			oplot, sex_radius*[0.9,0.9], [0,100], line=2, color=200
			oplot, sex_radius*[1.1,1.1], [0,100], line=2, color=200
			oplot, [0,100], sex_mag_range[0]*[1,1], line=2, color=100
			oplot, [0,100], sex_mag_range[1]*[1,1], line=2, color=100
			wait, 0.2
			if do_debug then wait, 2

		 	gv_match= where( min( (cat_sex[gv_sex].alpha_j2000#make_array(n_gv_ref,value=1,/double)-make_array(n_gv_sex,value=1.,/double)#cat_ref[*].ra)^2 + (cat_sex[gv_sex].delta_j2000#make_array(n_gv_ref,value=1,/double)-make_array(n_gv_sex,value=1.,/double)#cat_ref[*].dec)^2, id_match, DIM=1) LT (0.6/3600.)^2, n_gv_match)

		  if n_gv_match GT 2 then begin
				if do_debug then print, 'Running MATCH'
        gv_sex_match=gv_sex[(id_match[gv_match] mod n_gv_sex)]
        gv_ref_match=(id_match[gv_match]/n_gv_sex)

				case input_calib[i].filter of
					'u': begin
						gv=where(cat_ref[gv_ref_match].mag[0] GT 0 and cat_sex[gv_sex_match].flux_aper[4] GT 0., n_gv)

						mag_diff=cat_ref[gv_ref_match[gv]].mag[0] + 2.5*alog10( cat_sex[gv_sex_match[gv]].flux_aper[4]/im_exptime)
						pipe_zp=biweight_mean(mag_diff, pipe_zperr)
						plot, cat_ref[gv_ref_match].mag[0], mag_diff, psym=1, xrange=sex_mag_range+[0,2], yrange=median(mag_diff)+[-1.,1.]
						oplot, [0,100], pipe_zp*[1,1], color=200, line=2
						print, 'Pipeline zero point ', string(pipe_zp, pipe_zperr, FORMAT='(F0.2,4X,F0.2)')

						for k=0L, n_gv-1 do begin
							printf, lun_zp_full, j+1, cat_sex[gv_sex_match[gv[k]]].number, cat_sex[gv_sex_match[gv[k]]].alpha_j2000, cat_sex[gv_sex_match[gv[k]]].delta_j2000, cat_sex[gv_sex_match[gv[k]]].flux_aper[4], cat_sex[gv_sex_match[gv[k]]].flux_auto, cat_ref[gv_ref_match[gv[k]]].mag[0], mag_diff[k], pipe_airmass, FORMAT='(I4,4X,I4,4X,F0.6,4X,F0.6,4X,F15.3,4X,F15.3,4X,F0.3,4X,F0.3,4X,F0.3)'
						endfor
						end
					'g': begin
						gv=where(cat_ref[gv_ref_match].mag[1] GT 0 and cat_sex[gv_sex_match].flux_aper[4] GT 0., n_gv)

						mag_diff=cat_ref[gv_ref_match[gv]].mag[1] + 2.5*alog10( cat_sex[gv_sex_match[gv]].flux_aper[4]/im_exptime)
						pipe_zp=biweight_mean(mag_diff, pipe_zperr)
						plot, cat_ref[gv_ref_match].mag[1], mag_diff, psym=1, xrange=sex_mag_range+[0,2], yrange=median(mag_diff)+[-1.,1.]
						oplot, [0,100], pipe_zp*[1,1], color=200, line=2
						print, 'Pipeline zero point ', string(pipe_zp, pipe_zperr, FORMAT='(F0.2,4X,F0.2)')

						for k=0L, n_gv-1 do begin
							printf, lun_zp_full, j+1, cat_sex[gv_sex_match[gv[k]]].number, cat_sex[gv_sex_match[gv[k]]].alpha_j2000, cat_sex[gv_sex_match[gv[k]]].delta_j2000, cat_sex[gv_sex_match[gv[k]]].flux_aper[4], cat_sex[gv_sex_match[gv[k]]].flux_auto, cat_ref[gv_ref_match[gv[k]]].mag[1], mag_diff[k], pipe_airmass, FORMAT='(I4,4X,I4,4X,F0.6,4X,F0.6,4X,F15.3,4X,F15.3,4X,F0.3,4X,F0.3,4X,F0.3)'
						endfor
						end
					'i': begin
						gv=where(cat_ref[gv_ref_match].mag[3] GT 0 and cat_sex[gv_sex_match].flux_aper[4] GT 0., n_gv)

						mag_diff=cat_ref[gv_ref_match[gv]].mag[3] + 2.5*alog10( cat_sex[gv_sex_match[gv]].flux_aper[4]/im_exptime)
						pipe_zp=biweight_mean(mag_diff, pipe_zperr)
						plot, cat_ref[gv_ref_match].mag[3], mag_diff, psym=1, xrange=sex_mag_range+[0,2], yrange=median(mag_diff)+[-1.,1.]
						oplot, [0,100], pipe_zp*[1,1], color=200, line=2
						print, 'Pipeline zero point ', string(pipe_zp, pipe_zperr, FORMAT='(F0.2,4X,F0.2)')

						for k=0L, n_gv-1 do begin
							printf, lun_zp_full, j+1, cat_sex[gv_sex_match[gv[k]]].number, cat_sex[gv_sex_match[gv[k]]].alpha_j2000, cat_sex[gv_sex_match[gv[k]]].delta_j2000, cat_sex[gv_sex_match[gv[k]]].flux_aper[4], cat_sex[gv_sex_match[gv[k]]].flux_auto, cat_ref[gv_ref_match[gv[k]]].mag[3], mag_diff[k], pipe_airmass, FORMAT='(I4,4X,I4,4X,F0.6,4X,F0.6,4X,F15.3,4X,F15.3,4X,F0.3,4X,F0.3,4X,F0.3)'
						endfor
						end
					'z': begin
						gv=where(cat_ref[gv_ref_match].mag[4] GT 0 and cat_sex[gv_sex_match].flux_aper[4] GT 0., n_gv)

						mag_diff=cat_ref[gv_ref_match[gv]].mag[4] + 2.5*alog10( cat_sex[gv_sex_match[gv]].flux_aper[4]/im_exptime)
						pipe_zp=biweight_mean(mag_diff, pipe_zperr)
						plot, cat_ref[gv_ref_match].mag[4], mag_diff, psym=1, xrange=sex_mag_range+[0,2], yrange=median(mag_diff)+[-1.,1.]
						oplot, [0,100], pipe_zp*[1,1], color=200, line=2
						print, 'Pipeline zero point ', string(pipe_zp, pipe_zperr, FORMAT='(F0.2,4X,F0.2)')

						for k=0L, n_gv-1 do begin
							printf, lun_zp_full, j+1, cat_sex[gv_sex_match[gv[k]]].number, cat_sex[gv_sex_match[gv[k]]].alpha_j2000, cat_sex[gv_sex_match[gv[k]]].delta_j2000, cat_sex[gv_sex_match[gv[k]]].flux_aper[4], cat_sex[gv_sex_match[gv[k]]].flux_auto, cat_ref[gv_ref_match[gv[k]]].mag[4], mag_diff[k], pipe_airmass, FORMAT='(I4,4X,I4,4X,F0.6,4X,F0.6,4X,F15.3,4X,F15.3,4X,F0.3,4X,F0.3,4X,F0.3)'
						endfor
						end
					else: stop
				endcase

				printf, lun_zp, j+1, im_zp-2.5*alog10(im_exptime), pipe_airmass, pipe_zp, pipe_zperr, n_gv, im_fwhm, FORMAT='(I4,4X,F0.3,4X,F0.3,4X,F0.3,4X,F0.3,4X,I4,4X,F0.2)'
      endif
     
		endfor
		print, 'DECam zero point ', im_zp
		free_lun, lun_zp
		free_lun, lun_zp_full
		
	endfor

	gv_phot=where(input_calib.airmass LT 2. AND input_calib.mjd GT 0. AND input_calib.photflag EQ 'T', n_gv_phot)
	for i=0L, n_gv_phot-1 do begin
		print, 'Processing file ', input_calib[gv_phot[i]].sex_zp_file

		readcol, input_calib[gv_phot[i]].sex_zp_file, temp_ext, temp_decam_zp, temp_airmass, temp_zp, temp_zperr, temp_zp_nstars, temp_fwhm, FORMAT='I,F,F,F,F,I,F'
		create_struct, temp_struct, '', ['ext','zp','airmass','filter','mjd','photflag'], 'I,F,F,A,F,A', dim=n_elements(temp_ext)
		temp_struct.ext=temp_ext
		temp_struct.zp=temp_zp
		temp_struct.airmass=temp_airmass
		temp_struct.filter=input_calib[gv_phot[i]].filter
		temp_struct.mjd=input_calib[gv_phot[i]].mjd
		temp_struct.photflag=input_calib[gv_phot[i]].photflag
		if i EQ 0L then pipe_zp=temp_struct $
		else pipe_zp=[pipe_zp,temp_struct]

		readcol, input_calib[gv_phot[i]].sex_zp_full_file, temp_ext, temp_id, temp_ra, temp_dec, temp_flux_aper, temp_flux_auto, temp_mag_ref, temp_zp, temp_airmass, FORMAT='I,I,F,F,D,D,F,F,F'
		create_struct, temp_struct, '', ['ext','zp','airmass','filter','mjd'], 'I,F,F,A,F', dim=n_elements(temp_ext)
		temp_struct.ext=temp_ext
		temp_struct.zp=temp_zp
		temp_struct.airmass=temp_airmass
		temp_struct.filter=input_calib[gv_phot[i]].filter
		temp_struct.mjd=input_calib[gv_phot[i]].mjd
		if i EQ 0L then pipe_zp_full=temp_struct $
		else pipe_zp_full=[pipe_zp_full,temp_struct]

	endfor

	filter_uniq=pipe_zp[uniq(pipe_zp.filter, sort(pipe_zp.filter))].filter
	mjd_uniq=floor(pipe_zp[uniq(floor(pipe_zp.mjd+2./24), sort(floor(pipe_zp.mjd+2./24)))].mjd)

	openw, lun, 'survey_zp_'+do_program+'.dat', /get_lun
	printf, lun, '# MJD		filter		zp		k		zp_err		k_err		PHOTFLAG'
	for i=0L, n_elements(filter_uniq)-1 do begin
		for j=0L, n_elements(mjd_uniq)-1 do begin

			gv_pipe=where(pipe_zp.filter EQ filter_uniq[i] AND floor(pipe_zp.mjd+2./24) EQ mjd_uniq[j], n_gv_pipe)
			gv_pipe_phot=where(pipe_zp.filter EQ filter_uniq[i] AND floor(pipe_zp.mjd+2./24) EQ mjd_uniq[j] AND pipe_zp.photflag EQ 'T', n_gv_pipe_phot)
			coeff=robust_linefit(pipe_zp[gv_pipe_phot].airmass, pipe_zp[gv_pipe_phot].zp, temp_yfit, temp_sigma, coeff_error)
			photflag = coeff_error[0] LT 0.015 ? 'T' : 'F'
			printf, lun, mjd_uniq[j], filter_uniq[i], coeff[0], coeff[1], coeff_error[0], coeff_error[1], photflag, FORMAT='(F0.1,4X,A,4X,F0.2,4X,F0.2,4X,F0.2,4X,F0.2,4X,A)'

			cgloadct, 0
			cgwindow, wxsize=800, wysize=600
			cgplot, [0], [0], xrange=[1.,1.9], yrange=mean(pipe_zp[gv_pipe].zp)+[-1,1], /nodata, /window, xtitle='airmass', ytitle='mag_std - mag_inst', title='Filter '+filter_uniq[i]+' - Date '+date_conv(mjd_uniq[j]+2400000.5D,'S')
			cgplot, pipe_zp[gv_pipe].airmass, pipe_zp[gv_pipe].zp, psym=cgsymcat('filled circle'), symsize=0.4, color='blue', /over, /addcmd
			cgplot, pipe_zp[gv_pipe_phot].airmass, pipe_zp[gv_pipe_phot].zp, psym=cgsymcat('filled circle'), symsize=0.4, color='red', /over, /addcmd
			cgplot, [0.,100.], coeff[0]+coeff[1]*[0.,100.], color='black', line=2, /over, /addcmd
			cglegend, color='red', align=0, location=[0.14,0.86], length=0.02, titles=filter_uniq[i]+'-band coeff='+string(coeff,format='(F0.3,",",F0.3)')+' coeff_error='+string(coeff_error,format='(F0.3,",",F0.3)'), vspace=2., /addcmd
			cgcontrol, output=output_calib_dir+'/standard_zp_airmass_uncorrected_'+filter_uniq[i]+'_'+string(mjd_uniq[j],FORMAT='(I0)')+'.pdf'

			print, 'Airmass term compute - Filter ', filter_uniq[i], ' Date ', date_conv(mjd_uniq[j]+2400000.5D,'S')
			print, coeff[0], coeff_error[0], coeff[1], coeff_error[1], FORMAT='("zp = ",F0.3," +- ",F0.3," / k = ", F0.3, " +- ", F0.3)'

			gv_pipe=where(input_calib.filter EQ filter_uniq[i] AND floor(input_calib.mjd+2./24) EQ mjd_uniq[j], n_gv_pipe)
			gv_pipe_phot=where(input_calib.filter EQ filter_uniq[i] AND floor(input_calib.mjd+2./24) EQ mjd_uniq[j] AND input_calib.photflag EQ 'T', n_gv_pipe_phot)
			cgloadct, 0
			cgwindow, wxsize=800, wysize=600
			cgplot, [0], [0], xrange=mjd_uniq[j]+[-2,12]/24., yrange=[1.,1.8], /nodata, /window, xtitle='MJD', ytitle='airmass', title='Filter '+filter_uniq[i]+' - Date '+date_conv(mjd_uniq[j]+2400000.5D,'S'), XTICKFORMAT='(F0.1)'
			cgplot, input_calib[gv_pipe].mjd, input_calib[gv_pipe].airmass, psym=cgsymcat('filled circle'), symsize=1., color='blue', /over, /addcmd
			cgplot, input_calib[gv_pipe_phot].mjd, input_calib[gv_pipe_phot].airmass, psym=cgsymcat('filled circle'), symsize=1., color='red', /over, /addcmd
			cgcontrol, output=output_calib_dir+'/standard_airmass_mjd_'+filter_uniq[i]+'_'+string(mjd_uniq[j],FORMAT='(I0)')+'.pdf'

			cgdelete, /all

		endfor
	endfor
	free_lun, lun

	cgdelete, /all

	create_struct, pipe_phot, '', ['filter','zp','k'], 'A,F,F', dim=n_elements(filter_uniq)

endif else $
if recipe EQ 'sextractor' then begin

	if do_sky_method EQ 'median global' then begin
		input_target.im_file = repstr(input_target.ss_im_file,'.fits','.002.fits')
		input_target.sex_cat_file = repstr(input_target.ss_sex_cat_file,'.ldac','.002.ldac')
		input_target.sex_xml_file = repstr(input_target.ss_sex_xml_file,'.xml','.002.xml')
		input_target.scamp_head_file = repstr(input_target.ss_scamp_head_file,'.head','.002.head')
	endif else $
	if do_sky_method EQ 'gradient tps' then begin
		input_target.im_file = repstr(input_target.ss_im_file,'.fits','.003.fits')
		input_target.sex_cat_file = repstr(input_target.ss_sex_cat_file,'.ldac','.003.ldac')
		input_target.sex_xml_file = repstr(input_target.ss_sex_xml_file,'.xml','.003.xml')
		input_target.scamp_head_file = repstr(input_target.ss_scamp_head_file,'.head','.003.head')
	endif

	if do_overwrite then begin 
		file_delete, input_target.sex_cat_file, /noexpand, /allow_non, /quiet
	  file_delete, input_target.scamp_head_file, /noexpand, /allow_non, /quiet
	endif

	if n_elements(bridges) eq 0 then bridges = build_bridges(do_n_cpu)
	n_cpu = n_elements(bridges)

	for k=0L, n_cpu-1 do begin
		(bridges[k])->execute, '.r worker'
		(bridges[k])->execute, '.r callback'
		(bridges[k])->setproperty, callback='callback'
	endfor

	for i=0L, n_elements(input_target)-1 do begin
			
		print, 'Sexractor - Iteration '+strn(i+1)+' of '+strn(n_elements(input_target))

		if file_test(input_target[i].sex_cat_file, /regular) then fits_info, input_target[i].sex_cat_file, n_ext=sex_n_ext, /silent $
		else sex_n_ext=0

		if sex_n_ext NE input_target[i].n_chip*2 OR do_overwrite EQ 1 then begin

			case input_target[i].filter of
				'u': sex_satur_level='25000.'
				'g': sex_satur_level='25000.'
				'i': sex_satur_level='30000.'
				'z': sex_satur_level='30000.'
				else: stop
			endcase
		
;			if n_elements(input_zp) GT 0 then begin
;				gv=where(input_zp.filter EQ input_target[i].filter AND input_zp.mjd EQ input_target[i].mjd_floor, n_gv)
;				if n_gv GT 0 then sex_zp=input_zp[gv].zp $
;				else stop
;			endif else begin
;				sex_zp=input_target[i].zp
;			endelse
			sex_zp=input_target[i].zp

			command = 'sex '+input_target[i].im_file +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+input_target[i].sex_cat_file+' -WEIGHT_IMAGE '+input_target[i].weight_file+' -XML_NAME '+input_target[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(sex_zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_target[i].sex_check_file

;			ud = {xrange:temp_xrange[ii:ii+1], yrange:temp_yrange[jj:jj+1], pout:pout, im_size:im_size}
			bridge = get_idle_bridge(bridges)
;			bridge->setproperty, userdata=ud
			bridge->setvar, 'command', command

			print, 'Sexractor - Running command ', command
			bridge->execute, nowait=1, 'worker, command, out'

;			command = 'sex '+input_target[i].im_file +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+input_target[i].sex_cat_file+' -WEIGHT_IMAGE '+input_target[i].weight_file+' -XML_NAME '+input_target[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(sex_zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_target[i].sex_check_file
;			print, command
;			spawn, command

		endif

	endfor

	barrier_bridges, bridges

	burn_bridges, bridges
	temp=temporary(bridges)

endif else $
if recipe EQ 'iq' then begin

	window, 0, XSIZE=800, YSIZE=600
	window, 1, XSIZE=800, YSIZE=600


	for i=0L, n_elements(do_program)-1 do begin
		log_file='survey_iq_'+do_program[i]+'.dat'
		if file_test(log_file) EQ 0 OR do_overwrite then begin
			openw, lun, log_file, /get_lun
			printf, lun, '#  im_file        filter  tile  dither  sky_level   FWHM    TYPE'
			free_lun, lun
		endif
	endfor

	for i=0L, n_elements(input_target)-1 do begin
		print, 'IQ - Iteration '+strn(i+1)+' of '+strn(n_elements(input_target))
		print, 'IQ - Analyzing image quality of file '+input_target[i].im_file
		
		log_file='survey_iq_'+input_target[i].program+'.dat'
		readcol, log_file, temp_im_file, temp_filter, temp_tile, temp_dither, temp_sky, temp_fwhm, FORMAT='A,A,A,A,F,F'
		if n_elements(temp_im_file) GT 0 then begin
			gv_im_file=where(temp_im_file EQ input_target[i].im_file, n_gv_im_file)
			if n_gv_im_file GT 0 then continue
		endif

		vig_diam=101
		sex_mag_bin=0.5
		sex_radius_bin=0.2
		sex_flux_radius_min=1.4
		sex_ellipticity_min=0.8
		sex_flag_max=1
	
		case input_target[i].filter of
			'u': begin
				sex_mag_range=[10.,18.]
				plot_mag_range=[21,10]
			end
			'g': begin
				case input_target[i].type of
					'short': begin 
						sex_mag_range=[14.,19.]
						plot_mag_range=[22,13]
						end
					'long':	begin
						sex_mag_range=[15.,21.]
						plot_mag_range=[25,15]
						end
				endcase
			end
			'i': begin
				case input_target[i].type of
					'short': begin 
						sex_mag_range=[14.,19.]
						plot_mag_range=[22,13]
						end
					'long':	begin
						sex_mag_range=[15.,21.]
						plot_mag_range=[25,15]
						end
				endcase
			end
			else: stop
		endcase

		if file_test(input_target[i].sex_cat_file) EQ 0 then begin
			print, "IQ error - Please run decam_pipeline, 'sextractor'"
			stop
		endif

		iq_sky=list()
		iq_fwhm=list()

		foreach j, input_chip_iq do begin

			im_data=readfits(input_target[i].im_file, im_h, ext=j)
			im_size=size(im_data, /dim)
			im_gain=fxpar(im_h, 'GAINA')
			im_ron=fxpar(im_h, 'RDNOISEA')

			cat_sex=mrdfits(input_target[i].sex_cat_file, 2*j, cat_sex_h, COLUMNS=['NUMBER','X_IMAGE','Y_IMAGE','FLUX_RADIUS','MAG_AUTO','MAGERR_AUTO','FLAGS','A_IMAGE','B_IMAGE'], /silent)

			wset, 0
			plot, cat_sex.flux_radius, cat_sex.mag_auto, psym=1, xrange=[1,6], yrange=plot_mag_range
			oplot, [0,100], sex_mag_range[0]*[1,1], line=2, color=100
			oplot, [0,100], sex_mag_range[1]*[1,1], line=2, color=100
			gv_stars=where(cat_sex.mag_auto GT sex_mag_range[0] and cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_flux_radius_min AND cat_sex.flags LE 3 AND cat_sex.b_image/cat_sex.a_image GT sex_ellipticity_min, n_gv_stars)
			if n_gv_stars LT 10 then begin
				print, 'IQ - Error, there is not enough number of stars'
				continue
			endif
			oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=100


			plothist, cat_sex[gv_stars].flux_radius, temp_xhist, temp_yhist, bin=sex_radius_bin, /noplot
			temp=max(temp_yhist, gv) & sex_radius=temp_xhist[gv]
			sex_radius=median((cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 , n_gv)]).flux_radius)
			print, j+1, sex_radius, FORMAT='("Chip ",I2,"  Flux_radius ",F0.1)'

			oplot, sex_radius*[1,1], [0,100], color=200
			oplot, sex_radius*[0.9,0.9], [0,100], line=2, color=200
			oplot, sex_radius*[1.1,1.1], [0,100], line=2, color=200

			gv_stars=where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1]-1 AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 AND cat_sex.flags LE sex_flag_max AND cat_sex.x_image GT vig_diam/2. AND cat_sex.x_image LT (im_size[0]-vig_diam/2.) AND cat_sex.y_image GT vig_diam/2. AND cat_sex.y_image LT (im_size[1]-vig_diam/2.), n_gv_stars)
			oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=200
	
		  if n_gv_stars GT 0 then begin
				gv_stars=gv_stars[sort(cat_sex[gv_stars].mag_auto)]
	
				temp_n=n_gv_stars<5
	  	  temp_x=cat_sex[gv_stars].x_image
	   		temp_y=cat_sex[gv_stars].y_image

				im_fwhm=list()	
				im_sky=list()	

	  	  for k=0L, temp_n-1 do begin
		      x_range=[ floor(temp_x[k]-(vig_diam-1.)/2)>0, ceil(temp_x[k]+(vig_diam-1.)/2)<(im_size[0]-1) ]
	        y_range=[ floor(temp_y[k]-(vig_diam-1.)/2)>0, ceil(temp_y[k]+(vig_diam-1.)/2)<(im_size[1]-1) ]
	     	  vig_data = im_data[x_range[0]:x_range[1],y_range[0]:y_range[1]]; - im_sky
					vig_size=size(vig_data, /dim)
	
		      x_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
	    	  y_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
					vig_center_data=vig_data[x_center_range[0]:x_center_range[1],y_center_range[0]:y_center_range[1]]
					vig_center_size=size(vig_center_data, /dim)
	
		      im_max=max(vig_center_data, gv_max)
	       	im_c= [gv_max mod vig_center_size[0], gv_max/vig_center_size[0]]
					im_c += [x_center_range[0], y_center_range[0]]

					dist_circle, vig_mask, vig_size, im_c[0], im_c[1]
					vig_sky=median(vig_data[where(vig_mask GT 20., n_sky)])
					vig_mag=25.-2.5*alog10(total(vig_data[where(vig_mask LE 20., n_star)]) - n_star*vig_sky )

					vig_res=abs( vig_data-vig_sky - max(vig_center_data-vig_sky)/2. )
					gv=sort(vig_res*vig_mask^2)
					vig_fwhm=2*median(vig_mask[gv[1:5]])
					print, 'Radius for computing FWHM ', vig_mask[gv[1:5]]
					wset, 1
					plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
					oplot, [0,100], max(vig_data-vig_sky)/2.*[1,1], line=2, color=200
					oplot, vig_fwhm/2.*[1,1], [-1e5,1e5], line=2, color=200

					vig_skyrad=4*vig_fwhm < 40.
					vig_psfrad=3*vig_fwhm < 50.
					vig_fitrad=vig_fwhm < 50.
	
					gcntrd, vig_data, im_c[0], im_c[1], im_cx, im_cy, vig_fwhm
					dist_circle, vig_mask, vig_size, im_cx, im_cy
					vig_sky=median(vig_data[where(vig_mask GT vig_skyrad, n_sky)])
					vig_mag=25.-2.5*alog10(total(vig_data[where(vig_mask LE vig_skyrad, n_star)]) - n_star*vig_sky )
					plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
					oplot, [0,100], max(vig_data-vig_sky)/2.*[1,1], line=2, color=200
					oplot, vig_fwhm/2.*[1,1], [-1e5,1e5], line=2, color=200

		      fit_x_range=[ floor(im_cx- 2*vig_fwhm), ceil(im_cx+ 2*vig_fwhm) ]
	    	  fit_y_range=[ floor(im_cy- 2*vig_fwhm), ceil(im_cy+ 2*vig_fwhm) ]
					fit_vig_data=vig_data[fit_x_range[0]:fit_x_range[1],fit_y_range[0]:fit_y_range[1]]
					fit_vig_size=size(fit_vig_data, /dim)
		      fit_peak=max(fit_vig_data-vig_sky, gv_max)
	       	im_c= [gv_max mod fit_vig_size[0], gv_max/fit_vig_size[0]]
					gcntrd, fit_vig_data, im_c[0], im_c[1], fit_im_cx, fit_im_cy, vig_fwhm

					fit_im_cx=fit_im_cx + (floor(im_cx- 2*vig_fwhm)-floor(im_cx- 1.2*vig_fwhm))
					fit_im_cy=fit_im_cy + (floor(im_cy- 2*vig_fwhm)-floor(im_cy- 1.2*vig_fwhm))
		      fit_x_range=[ floor(im_cx- 1.2*vig_fwhm), ceil(im_cx+ 1.2*vig_fwhm) ]
	    	  fit_y_range=[ floor(im_cy- 1.2*vig_fwhm), ceil(im_cy+ 1.2*vig_fwhm) ]
					fit_vig_data=vig_data[fit_x_range[0]:fit_x_range[1],fit_y_range[0]:fit_y_range[1]]

;					fit_sigma_range=vig_fwhm/2.35482*[0.8,1.4]
;					fit_ellipticity_range=[0.95,1.05]
;					fit_peak_range=[0.9,1.1]
;					fit_niter=30
;					fit_res=fltarr(fit_niter*[1,1,1])
;
;					fit_peak=(fit_peak_range[0] + findgen(fit_niter)/(fit_niter-1.)*(fit_peak_range[1]-fit_peak_range[0]))*max(fit_vig_data-vig_sky)
;					fit_sigma=fit_sigma_range[0] + findgen(fit_niter)/(fit_niter-1.)*(fit_sigma_range[1]-fit_sigma_range[0])
;					fit_ellipticity=fit_ellipticity_range[0] + findgen(fit_niter)/(fit_niter-1.)*(fit_ellipticity_range[1]-fit_ellipticity_range[0])
;
;					for ii=0L, fit_niter-1 do begin
;						for jj=0L, fit_niter-1 do begin
;							for kk=0L, fit_niter-1 do begin
;								psf_param=[ [fit_peak[ii]*[1,1]],[[fit_im_cx,fit_im_cy]],[fit_sigma[jj]*[fit_ellipticity[kk],1]] ]
;								fit_res[ii,jj,kk]=total( abs(fit_vig_data-vig_sky - psf_gaussian(psf_param, npix=fit_vig_size)) )
;							endfor
;						endfor
;					endfor
;					temp=min(fit_res, gv_min)
;					gv_min=array_indices(fit_res, gv_min)
;					psf_param1=[ fit_peak[gv_min[0]], fit_im_cx, fit_im_cy, fit_sigma[gv_min[1]]*[fit_ellipticity[gv_min[2]],1] ]

          parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 7)
          parinfo[*].value = [0.,fit_peak, vig_fwhm, vig_fwhm, fit_im_cx, fit_im_cy, 0.]
          parinfo[0].fixed = 1
          temp=mpfit2dpeak(fit_vig_data-vig_sky, psf_param1, /GAUSSIAN, /TILT, ESTIMATES=parinfo.value, PARINFO=parinfo)

		      im_fwhm.add, 2*sqrt(2*alog(2))*sqrt( (psf_param1[2]^2 + psf_param1[3]^2)/2. )
					im_sky.add, vig_sky
					print, 'FWHM and SKY: ', im_fwhm[-1], ' pixels', im_sky[-1], ' ADU'

;					getpsf, vig_data, im_cx, im_cy, vig_mag, vig_sky, im_ron, im_gain, psf_param2, psf_residuals, [0], vig_psfrad, vig_fitrad, output_sex_dir+'/im_sextractor_psf.fits' 
;		      im_fwhm.add, 2*sqrt(2*alog(2))*sqrt( (psf_param2[3]^2 + psf_param2[4]^2)/2. )
;					im_sky.add, vig_sky
;					print, 'FWHM and SKY: ', im_fwhm[-1], ' pixels', im_sky[-1], ' ADU'
	
					dist_circle, vig_mask, vig_size, psf_param1[4]+fit_x_range[0], psf_param1[5]+fit_y_range[0]
					plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
					x=findgen(100)/10.
					y1=gaussian(x,[psf_param1[1], 0., mean(psf_param1[2:3])])
;					y2=gaussian(x,[psf_param2[0], 0., mean(psf_param2[3:4])])
					oplot, x, y1, line=2, color=100
;					oplot, x, y2, line=2, color=200
					oplot, im_fwhm[-1]/2.*[1,1], [-1e5,1e5], line=2, color=200
					wait, 0.2
						
	    	endfor
				im_sky=median(im_sky.toarray(type='float'))
				im_fwhm=median(im_fwhm.toarray(type='float'))

				iq_sky.add, im_sky
				iq_fwhm.add, im_fwhm

			endif	

		endforeach

		print, 'Summary of SKY and FWHM'
		print, iq_sky
		print, iq_fwhm
		iq_sky=mean(iq_sky.toarray(type='float'))
		iq_fwhm=mean(iq_fwhm.toarray(type='float'))

		openw, lun, log_file, /get_lun, /append
		printf, lun, input_target[i].im_file, input_target[i].filter, input_target[i].tile, input_target[i].dither, iq_sky, iq_fwhm, input_target[i].type, FORMAT='(A,4X,A,4X,A,4X,A,4X,F0.1,4X,F0.1,4X,A)'
		free_lun, lun

	endfor

endif else $
if recipe EQ 'scamp' then begin

	if do_sky_method EQ 'median global' then begin
		input_target.im_file = repstr(input_target.ss_im_file,'.fits','.002.fits')
		input_target.sex_cat_file = repstr(input_target.ss_sex_cat_file,'.ldac','.002.ldac')
		input_target.scamp_cat_file = repstr(input_target.ss_scamp_cat_file,'.ldac','.002.ldac')
		input_target.scamp_head_file = repstr(input_target.ss_scamp_head_file,'.head','.002.head')
		input_target.scamp_ahead_file = repstr(input_target.ss_scamp_ahead_file,'.ahead','.002.ahead')
	endif else $
	if do_sky_method EQ 'gradient tps' then begin
		input_target.im_file = repstr(input_target.ss_im_file,'.fits','.003.fits')
		input_target.sex_cat_file = repstr(input_target.ss_sex_cat_file,'.ldac','.003.ldac')
		input_target.scamp_cat_file = repstr(input_target.ss_scamp_cat_file,'.ldac','.003.ldac')
		input_target.scamp_head_file = repstr(input_target.ss_scamp_head_file,'.head','.003.head')
		input_target.scamp_ahead_file = repstr(input_target.ss_scamp_ahead_file,'.ahead','.003.ahead')
	endif

  scamp_pos_error='1.2'
  scamp_scale_error='1.02'
  scamp_angle_error='0.02'
  scamp_sn_thresholds='40.,80.0'
	scamp_fwhm_thresholds='0.,100.'
	scamp_crossid_radius='4.'
	scamp_distort_degrees='4'
	scamp_astref_catalog = '2MASS' ;'USNO-B1' ;'2MASS'
	scamp_astref_band = 'Ks';'Rf' ;'DEFAULT' ;'Ks'
	scamp_astrefmag_limits='6.,20.' ;'-99.,18.' 
	scamp_match_resol='0.'

  scamp_refcat_dir = output_scamp_dir[0]+'/refcat'
	scamp_dir_out=output_scamp_dir[0]+'/PROGRAM_'+do_program_plain+'_FILTER_'+do_filter_plain+'_ORDER_'+scamp_distort_degrees+'_REFCAT_'+scamp_astref_catalog
	scamp_list_file=scamp_dir_out+'/scamp_decam_'+do_program_plain+'_'+do_filter_plain+'.dat'
	scamp_cat_file_out = scamp_dir_out+'/scamp_decam_'+do_program_plain+'_'+do_filter_plain+'.ldac'
  scamp_cat_type_out = 'FITS_LDAC'
  scamp_xml_file = scamp_dir_out+'/scamp_decam_'+do_program_plain+'_'+do_filter_plain+'.xml'
  scamp_check_type = 'SKY_ALL,FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR'
  scamp_check_file = strjoin(scamp_dir_out+'/'+['sky_all','fgroups','distort','astr_interror2d','astr_interror1d','astr_referror2d','astr_referror1d','astr_chi2','psphot_error'],',')


  if do_overwrite then begin
		file_delete, input_target.scamp_cat_file, /noexpand, /allow_non, /quiet
  	file_delete, input_target.scamp_head_file, /noexpand, /allow_non, /quiet
	endif
	if not file_test(scamp_refcat_dir, /directory) then file_mkdir, scamp_refcat_dir, /noexpand_path
	if not file_test(scamp_dir_out, /directory) then file_mkdir, scamp_dir_out, /noexpand_path
	loadct, 12

	for i=0L, n_elements(input_target)-1 do begin
		print, 'SCAMP - Creating file '+input_target[i].scamp_cat_file

		scamp_mag_bin=0.5
		scamp_radius_bin=0.2
		scamp_flux_radius_min=1.4
		scamp_ellipticity_min=0.8
		scamp_flag_max=1
		case input_target[i].filter of
			'u': begin
				scamp_mag_range=[10.,18.]
				plot_mag_range=[21,10]
  			scamp_sn_thresholds='20.,40.0'
				scamp_distort_degrees='4'
			end
			'g': begin
				scamp_mag_range=[15.,21.]
				plot_mag_range=[25,15]
  			scamp_sn_thresholds='20.,40.0'
				scamp_distort_degrees='4'
			end
			'i': begin
				scamp_mag_range=[15.,21.]
				plot_mag_range=[25,15]
  			scamp_sn_thresholds='20.,40.0'
				scamp_distort_degrees='4'
			end
			else: stop
		endcase
	
		if file_test(input_target[i].sex_cat_file) EQ 0 then stop

		if file_test(input_target[i].scamp_cat_file, /regular) then fits_info, input_target[i].scamp_cat_file, n_ext=scamp_n_ext, /silent $
		else scamp_n_ext=0
		if scamp_n_ext EQ input_target[i].n_chip*2 AND do_overwrite EQ 0 then continue

		fits_open, input_target[i].sex_cat_file, fcb_in  ; Lee la tabla fits sin intervenerla           
	  fits_read, fcb_in, cat_data0, cat_h0, exten=0 
  	writefits, input_target[i].scamp_cat_file, cat_data0, cat_h0

;		rdfits_struct, input_target[i].sex_cat_file, cat_sex, /silent
		for j=0L, input_target[i].n_chip-1 do begin
			cat_sex=mrdfits(input_target[i].sex_cat_file, 2*(j+1), cat_sex_h, COLUMNS=['NUMBER','FLUX_RADIUS','MAG_AUTO','MAGERR_AUTO','FLAGS','A_IMAGE','B_IMAGE'], /silent)
			plot, cat_sex.flux_radius, cat_sex.mag_auto, psym=1, xrange=[1,6], yrange=plot_mag_range
			oplot, [0,100], scamp_mag_range[0]*[1,1], line=2, color=100
			oplot, [0,100], scamp_mag_range[1]*[1,1], line=2, color=100

			gv_stars=where(cat_sex.mag_auto GT scamp_mag_range[0] and cat_sex.mag_auto LT scamp_mag_range[1] AND cat_sex.flux_radius GT scamp_flux_radius_min AND cat_sex.flags LE 3 AND cat_sex.b_image/cat_sex.a_image GT scamp_ellipticity_min, n_gv_stars)
			if n_gv_stars LT 10 then begin
				print, 'SCAMP - Error, there is not enough number of stars'
				stop
			endif
			oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=100

			plothist, cat_sex[gv_stars].flux_radius, temp_xhist, temp_yhist, bin=scamp_radius_bin, /noplot
			temp_yhist_weight=fltarr(n_elements(temp_yhist))
			for k=0L, n_elements(temp_xhist)-1 do begin
				gv=where( abs(cat_sex[gv_stars].flux_radius-temp_xhist[k]) LE scamp_radius_bin/2., n_gv)
				if n_gv GT 0 then temp_yhist_weight[k] = max(cat_sex[gv_stars[gv]].mag_auto)-min(cat_sex[gv_stars[gv]].mag_auto)
			endfor

			temp=max(temp_yhist*temp_yhist_weight, gv) & scamp_radius=temp_xhist[gv]
			scamp_radius=median((cat_sex[where(cat_sex.mag_auto GT scamp_mag_range[0] AND cat_sex.mag_auto LT scamp_mag_range[1] AND cat_sex.flux_radius GT scamp_radius*0.9 AND cat_sex.flux_radius LT scamp_radius*1.1 , n_gv)]).flux_radius)
			print, j+1, scamp_radius, FORMAT='("Chip ",I2,"  Flux_radius ",F0.1)'

			oplot, scamp_radius*[1,1], [0,100], color=200
			oplot, scamp_radius*[0.9,0.9], [0,100], line=2, color=200
			oplot, scamp_radius*[1.1,1.1], [0,100], line=2, color=200
			wait, 0.2
			if j EQ 0 AND do_debug then begin
				forprint, temp_xhist, temp_yhist, text=2
				print, 'DEBUG - CHECK mag range and flux radius are ok'
				stop
			endif

			gv_stars=where(cat_sex.mag_auto GT scamp_mag_range[0] AND cat_sex.mag_auto LT scamp_mag_range[1] AND cat_sex.flux_radius GT scamp_radius*0.9 AND cat_sex.flux_radius LT scamp_radius*1.1 AND cat_sex.flags LE scamp_flag_max AND cat_sex.b_image/cat_sex.a_image GT scamp_ellipticity_min, n_gv_stars)
			oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=200
			wait, 0.2

  		fits_read, fcb_in, cat_data1, cat_h1, exten=2*j+1                                   
			fits_read, fcb_in, cat_data2, cat_h2, exten=2*j+2

			cat_data1=reform(cat_data1, [n_elements(cat_data1),1])
			cat_data2=cat_data2[*,[gv_stars]]
			fxaddpar, cat_h2, 'NAXIS2', (size(cat_data2, /dim))[1]

  		fits_open, input_target[i].scamp_cat_file, fcb_out, /append
		  fits_write, fcb_out, cat_data1, cat_h1
			fits_write, fcb_out, cat_data2, cat_h2
  		fits_close, fcb_out

		endfor
  	fits_close, fcb_in
	endfor

	forprint, input_target.scamp_cat_file, textout=scamp_list_file, FORMAT='(A)', /NOCOMMENT

	for i=0L, n_elements(input_target)-1 do begin

		if file_test(input_target[i].scamp_ahead_file) EQ 0 OR do_overwrite OR do_ahead then begin

			im_h=headfits(input_target[i].im_file)
			im_mjd=fxpar(im_h, 'MJD-OBS')
			input_target[i].mjd=im_mjd
			print, input_target[i].im_file, input_target[i].filter, im_mjd, FORMAT='(A,2X,A,2X,F0.4)'

			openw, lun, input_target[i].scamp_ahead_file, /get_lun
			for j=0L, input_target[i].n_chip-1 do begin
				im_h=headfits(input_target[i].im_file, ext=j+1)
				extast, im_h, im_ast
				xy2ad, im_ast.naxis[0]/2., im_ast.naxis[1]/2., im_ast, im_ra, im_dec
				im_airmass=tai2airmass(im_ra, im_dec, 2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)

				printf, lun, "PROGRAM = '"+input_target[i].program+"'"
				printf, lun, "TYPE    = '"+input_target[i].type+"'"
				printf, lun, "FILTER  = '"+input_target[i].filter+"'"
				printf, lun, 'AIRMASS =     '+string(im_airmass,FORMAT='(F0.3)')
				printf, lun, 'EXPTIME =     '+string(input_target[i].exptime,FORMAT='(F0.2)')  ;string(600.,FORMAT='(F0.2)')
				if n_elements(input_zp) GT 0 then begin
					gv=where(input_zp.filter EQ input_target[i].filter AND input_zp.mjd EQ floor(input_target[i].mjd+2./24), n_gv)
					if n_gv EQ 1 then begin
						printf, lun, 'PHOTFLAG=     '+input_zp[gv].photflag
						printf, lun, 'PHOT_ZP =     '+string(input_zp[gv].zp,FORMAT='(F0.2)')
						printf, lun, 'PHOT_K  =     '+string(input_zp[gv].k,FORMAT='(F0.2)')
					endif	else begin
						gv=where(input_zp.filter EQ input_target[i].filter, n_gv)
						printf, lun, 'PHOTFLAG=     F'
						printf, lun, 'PHOT_ZP =     '+string(input_target[i].zp,FORMAT='(F0.3)')
						printf, lun, 'PHOT_K  =     '+string(input_zp[gv[0]].k,FORMAT='(F0.2)')
					endelse
				endif else begin
					printf, lun, 'PHOTFLAG=     T'
					printf, lun, 'PHOT_ZP =     '+string(input_target[i].zp,FORMAT='(F0.3)')
					printf, lun, 'PHOT_K  =      '+string(0.,FORMAT='(F0.3)')
				endelse
				printf, lun, 'END'
			endfor
			free_lun, lun

		endif
	endfor

	command='scamp @'+scamp_list_file+' -c scamp_config/ctio_decam.scamp'+' -MERGEDOUTCAT_TYPE '+scamp_cat_type_out+' -MERGEDOUTCAT_NAME '+scamp_cat_file_out+' -MATCH Y -WRITE_XML Y -XML_NAME '+scamp_xml_file+' -SAVE_REFCATALOG Y -REFOUT_CATPATH '+scamp_refcat_dir+' -CHECKPLOT_DEV PSC -CHECKPLOT_ANTIALIAS Y -CHECKPLOT_TYPE '+scamp_check_type + ' -CHECKPLOT_NAME '+scamp_check_file + ' -ASTREF_CATALOG '+scamp_astref_catalog+' -ASTREF_BAND '+scamp_astref_band+' -ASTREFMAG_LIMITS '+scamp_astrefmag_limits+' -DISTORT_DEGREES '+scamp_distort_degrees+' -PHOTCLIP_NSIGMA 2. -SOLVE_ASTROM Y -SOLVE_PHOTOM Y -POSITION_MAXERR '+scamp_pos_error+' -PIXSCALE_MAXERR '+scamp_scale_error+' -POSANGLE_MAXERR '+scamp_angle_error+' -SN_THRESHOLDS '+scamp_sn_thresholds+' -FWHM_THRESHOLDS '+scamp_fwhm_thresholds+' -CROSSID_RADIUS '+scamp_crossid_radius+' -MATCH_RESOL '+scamp_match_resol

	print, command
	spawn, command

endif else $
if recipe EQ 'swarp' then begin

	if do_sky_method EQ 'median global' then begin
		input_target.ss_weight_file = repstr(input_target.ss_weight_file,'.WEIGHT.fits','.002.WEIGHT.fits')
		input_target.scamp_head_file = repstr(input_target.ss_scamp_head_file,'.head','.002.head')
		input_target.swarp_im_file = repstr(input_target.ss_swarp_im_file,'.fits','.002.fits')
		input_target.swarp_head_file = repstr(input_target.ss_swarp_head_file,'.head','.002.head')
	endif else $
	if do_sky_method EQ 'gradient tps' then begin
		input_target.ss_weight_file = repstr(input_target.ss_weight_file,'.WEIGHT.fits','.003.WEIGHT.fits')
		input_target.scamp_head_file = repstr(input_target.ss_scamp_head_file,'.head','.003.head')
		input_target.swarp_im_file = repstr(input_target.ss_swarp_im_file,'.fits','.003.fits')
		input_target.swarp_head_file = repstr(input_target.ss_swarp_head_file,'.head','.003.head')
	endif

	tile_uniq=input_target[uniq(input_target.tile, sort(input_target.tile))].tile
	filter_uniq=input_target[uniq(input_target.filter, sort(input_target.filter))].filter
	type_uniq=input_target[uniq(input_target.type, sort(input_target.type))].type


	for i=0L, n_elements(tile_uniq)-1 do begin
		for k=0L, n_elements(type_uniq)-1 do begin

			if do_align_filter NE '' then begin
	
				swarp_list_file=output_stack_swarp_dir+'/swarp_header_tile'+tile_uniq[i]+'_'+do_align_filter+'_'+type_uniq[k]+'.lst'
				swarp_im_out=output_stack_swarp_dir+'/swarp_header_tile'+tile_uniq[i]+'_'+do_align_filter+'_'+type_uniq[k]+'.fits'
				swarp_weight_out=output_stack_swarp_dir+'/swarp_header_tile'+tile_uniq[i]+'_'+do_align_filter+'_'+type_uniq[k]+'.WEIGHT.fits'
				swarp_combine_type='CLIPPED'
				swarp_resampling_type='LANCZOS2'

				gv=where(input_target.tile EQ tile_uniq[i] AND input_target.type EQ type_uniq[k] AND input_target.filter EQ do_align_filter, n_gv)
				if n_gv GT 0 AND ( file_test(swarp_im_out) EQ 0 OR do_overwrite ) then begin

					for l=0L, n_gv-1 do begin
						file_copy, input_target[gv[l]].scamp_head_file, input_target[gv[l]].swarp_head_file, /OVERWRITE

 						file_delete, input_target[gv[l]].ss_weight_file, /noexpand, /allow_non, /quiet
						command='ln -s '+ input_target[gv[l]].weight_file+' '+ input_target[gv[l]].ss_weight_file
						print, command
						spawn, command
					endfor
	
					forprint, input_target[gv].swarp_im_file, textout=swarp_list_file, FORMAT='(A)', /NOCOMMENT
					command='swarp @'+swarp_list_file+' -c swarp_config/ctio_decam.swarp'+' -IMAGEOUT_NAME '+swarp_im_out+' -WEIGHTOUT_NAME '+swarp_weight_out+' -COMBINE_TYPE '+swarp_combine_type+' -RESAMPLE Y -VERBOSE_TYPE LOG -HEADER_ONLY Y'
					print, command
					spawn, command

					im_h=headfits(swarp_im_out)
				endif else $
				if n_gv GT 0 then begin
					print, 'SWARP - Reading header from align filter image'
					im_h=headfits(swarp_im_out)
				endif else $	
					print, 'SWARP - There are no images to create align filter image'
	
			endif
	
			for j=0L, n_elements(filter_uniq)-1 do begin
	
				swarp_combine_type='CLIPPED'
				swarp_resampling_type='LANCZOS2'
				swarp_resample_do='Y'
				swarp_resample_dir=output_stack_swarp_dir+'/resample'
				swarp_weight_suffix='.WEIGHT.fits'
				swarp_verbose='NORMAL'

				if do_align_filter NE '' then begin
					if filter_uniq[j] NE do_align_filter then do_align='_ALIGN'+do_align_filter $
					else do_align=''
				endif else $
					do_align=''

				if do_ndither NE 0 then do_dither='_NDITHER'+string(do_ndither, format='(I0)') $
				else do_dither=''
	
				if do_sky_method EQ 'median global' then begin
					swarp_list_file=output_stack_swarp_dir+'/swarp_ss_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.dat'
					swarp_xml_file=output_stack_swarp_dir+'/swarp_ss_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.xml'
					swarp_im_out=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.fits'
					swarp_weight_out=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.WEIGHT.fits'
					sex_stack_cat_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.ldac'
					sex_stack_xml_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.xml'
;					sex_stack_checkimage_file=strjoin(output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+['.002.CHECK_BACK.fits','.002.CHECK_SEGMENTATION.fits'],',')
					sex_stack_checkimage_type='SEGMENTATION'
					sex_stack_checkimage_file=output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.CHECK_SEGMENTATION.fits'
				endif else $
				if do_sky_method EQ 'gradient tps' then begin
					swarp_list_file=output_stack_swarp_dir+'/swarp_ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+'.003.dat'
					swarp_xml_file=output_stack_swarp_dir+'/swarp_ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+'.003.xml'
					swarp_im_out=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+'.003.fits'
					swarp_weight_out=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+'.003.WEIGHT.fits'
					swarp_im_head=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+'.003.head'
					sex_stack_cat_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+'.003.ldac'
					sex_stack_xml_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+'.003.xml'
;					sex_stack_checkimage_file=strjoin(output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+['.003.CHECK_BACK.fits','.003.CHECK_SEGMENTATION.fits'],',')
					sex_stack_checkimage_type='SEGMENTATION'
					sex_stack_checkimage_file=output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+'.003.CHECK_SEGMENTATION.fits'
					sex_stack_mask_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+do_dither+'.003.MASK.fits'
				endif else begin
					swarp_list_file=output_stack_swarp_dir+'/swarp_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.dat'
					swarp_xml_file=output_stack_swarp_dir+'/swarp_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.xml'
					swarp_im_out=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.fits'
					swarp_weight_out=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.WEIGHT.fits'
					sex_stack_cat_file=output_stack_sex_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.ldac'
					sex_stack_xml_file=output_stack_sex_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.xml'
;					sex_stack_checkimage_file=strjoin(output_stack_check_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+['.CHECK_BACK.fits','.CHECK_SEGMENTATION.fits'],',')
					sex_stack_checkimage_type='SEGMENTATION'
					sex_stack_checkimage_file=output_stack_check_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.CHECK_SEGMENTATION.fits'
				endelse		

				sex_stack_im_file=swarp_im_out
				sex_stack_weight_file=swarp_weight_out

				if do_align NE '' then begin
					if n_elements(im_h) GT 0 then begin
						openw, lun, swarp_im_head, /get_lun
						printf, lun, im_h, FORMAT='(A)'
						close, lun, /all
					endif
				endif

				gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.type EQ type_uniq[k], n_gv)
				if do_ndither NE 0 then begin
					n_gv=n_gv<do_ndither
					gv=gv[0:n_gv-1]
				endif

				if n_gv GT 0 then begin

					if file_test(swarp_im_out) EQ 0 OR do_overwrite EQ 1 then begin

						for l=0L, n_gv-1 do begin
							file_copy, input_target[gv[l]].scamp_head_file, input_target[gv[l]].swarp_head_file, /OVERWRITE

  						file_delete, input_target[gv[l]].ss_weight_file, /noexpand, /allow_non, /quiet
							command='ln -s '+ input_target[gv[l]].weight_file+' '+ input_target[gv[l]].ss_weight_file
							print, command
							spawn, command
						endfor
						if not file_test(swarp_resample_dir, /directory) then file_mkdir, swarp_resample_dir, /noexpand_path

						forprint, input_target[gv].swarp_im_file, textout=swarp_list_file, FORMAT='(A)', /NOCOMMENT
	
						command='swarp @'+swarp_list_file+' -c swarp_config/ctio_decam.swarp'+' -IMAGEOUT_NAME '+swarp_im_out+' -WEIGHTOUT_NAME '+swarp_weight_out+' -COMBINE_TYPE '+swarp_combine_type+' -RESAMPLE '+swarp_resample_do+' -RESAMPLE_DIR '+swarp_resample_dir+' -SATLEV_DEFAULT 35000 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX '+swarp_weight_suffix+' -WEIGHT_THRESH 0. -RESCALE_WEIGHTS N -BLANK_BADPIXELS Y -WRITE_XML Y -XML_NAME '+swarp_xml_file+' -VERBOSE_TYPE ' +swarp_verbose+' -RESAMPLING_TYPE '+swarp_resampling_type+' -SUBTRACT_BACK '+ (do_sky_method EQ '' ? 'Y':'N') + ' -BACK_SIZE 384'
						print, command
						spawn, command
					endif

					if file_test(swarp_im_head) then file_delete, swarp_im_head

					if file_test(sex_stack_cat_file) EQ 0 OR file_test(sex_stack_checkimage_file) EQ 0 OR do_overwrite EQ 1 then begin
						command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file+' -BACK_SIZE 384'
						print, command
						spawn, command
					endif

					if file_test(sex_stack_mask_file) EQ 0 OR do_overwrite EQ 1 then begin
						im_data=readfits_big(sex_stack_checkimage_file, im_h)
  	  			writefits, sex_stack_mask_file, (im_data GT 0), im_h
					endif

				endif

			endfor
		endfor
	endfor

endif else $
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
					stack_im_fix_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_FIX.003.fits' ;output_stack_check_dir+'/ss_fornax_FIX.003.fits'
					stack_weight_fix_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_FIX.003.WEIGHT.fits'  ;output_stack_check_dir+'/ss_fornax_FIX.003.WEIGHT.fits'
					stack_im_convol_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_FIX_CONVOL.003.fits'  ;output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_CONVOL.003.fits'

          sex_stack_cat_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.ldac'
          sex_stack_xml_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.xml'
          sex_stack_checkimage_file=output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.CHECK_SEGMENTATION.fits'

					stack_im_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.fits'
					stack_weight_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.WEIGHT.fits'
          stack_cat_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.ldac'
					stack_xml_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.xml'
          stack_checkimage_sex_file=output_stack_check_dir+'/ss_fornax_SEX.003.CHECK_SEGMENTATION.fits'
					stack_checkimage_sex_type='SEGMENTATION'

					sex_stack_cat_fix_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_FIX.003.ldac'
					sex_stack_xml_fix_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_FIX.003.xml'
					sex_stack_checkimage_fix_file=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_FIX.003.CHECK_SEGMENTATION.fits'
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
	
					command='sex '+stack_im_fix_file+' -c sex_config/ctio_decam_stack.sex -CATALOG_NAME '+sex_stack_cat_fix_file+' -WEIGHT_IMAGE '+stack_weight_fix_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_fix_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_fix_file+' -BACK_TYPE MANUAL -BACK_VALUE 0. -FILTER_NAME sex_config/gauss_3.0_7x7.conv -DETECT_THRESH 1.5 -ANALYSIS_THRESH 2. -DEBLEND_NTHRESH 16 -DEBLEND_MINCONT 0.01'
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


endif else $
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

				stack_mask_long_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_long.MASK.fits'

				stack_im_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+'.fits'
				stack_weight_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+'.WEIGHT.fits'
				stack_check_file=output_stack_check_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+'.CHECK_SEGMENTATION.fits'
				stack_mask_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+'.MASK.fits'

				if file_test(stack_mask_long_file, /regular) EQ 1 then begin
					command='ln -s '+stack_mask_long_file+' '+stack_mask_file
					print, command
					spawn, command
				endif

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
			stack_master_mask_long_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_long.MASTER_MASK.fits'

			stack_mask_file=file_search(output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_?_'+type_uniq[k]+'.MASK.fits', COUNT=n_stack_mask)
			stack_master_mask_file=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+type_uniq[k]+'.MASTER_MASK.fits'
			swarp_list_file=output_stack_swarp_dir+'/swarp_fornax_tile'+tile_uniq[i]+'_'+type_uniq[k]+'.MASTER_MASK.lst'
		
			if file_test(stack_master_mask_long_file, /regular) EQ 1 then begin
				command='ln -s '+stack_master_mask_long_file+' '+stack_master_mask_file
				print, command
				spawn, command
			endif

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
			'gradient tps': begin
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

			plot_window=0

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
		print, 'Saving file ', ss_im_file
		mwrfits, 3, ss_im_file, im_h, /create
		print, 'Saving file ', sky_im_file
		mwrfits, 3, sky_im_file, im_h, /create
	
		for l=0L, input_target[i].n_chip-1 do begin
			im_h=headfits(input_target[i].im_file, exten=l+1)
;			mwrfits, (*pout_im_masked)[*,*,l], masked_im_file, im_h, /silent
;			mwrfits, (*pout_im_masked_res)[*,*,l], masked_res_im_file, im_h, /silent
;			mwrfits, (*pout_im_masked_median)[*,*,l], masked_median_im_file, im_h, /silent
			mwrfits, (*pout_im)[*,*,l], ss_im_file, im_h, /silent
			mwrfits, (*pout_sky)[*,*,l], sky_im_file, im_h, /silent
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
if recipe EQ 'psfex' then begin

  vig_diam=101
  sex_mag_bin=0.5
  sex_radius_bin=0.2
  sex_flux_radius_min=1.5
  sex_ellipticity_min=0.8
  sex_flag_max=0

  psfex_psfvar_nsnap='20'
	psfex_psfvar_keys='X_IMAGE,Y_IMAGE';'X_IMAGE,Y_IMAGE';'X_IMAGE,Y_IMAGE,MAG_AUTO'
	psfex_psfvar_groups='1,1';'1,1';'1,1,2'
  psfex_psfvar_degrees='5';'5';'5,1'
  psfex_basis_type='PIXEL'
  psfex_basis_number='30'
  psfex_psf_size='43,43';'55,55'
  psfex_psf_sampling='1.0';'0.8'

	tile_uniq=input_target[uniq(input_target.tile, sort(input_target.tile))].tile
	filter_uniq=input_target[uniq(input_target.filter, sort(input_target.filter))].filter
	type_uniq=input_target[uniq(input_target.type, sort(input_target.type))].type

	for i=0L, n_elements(tile_uniq)-1 do begin
		for k=0L, n_elements(type_uniq)-1 do begin
			for j=0L, n_elements(filter_uniq)-1 do begin

				if do_align_filter NE '' then begin
					if filter_uniq[j] NE do_align_filter then do_align='_ALIGN'+do_align_filter $
					else do_align=''
				endif else $
					do_align=''
	
				if do_sky_method EQ 'median global' then begin
					swarp_list_file=output_stack_swarp_dir+'/swarp_ss_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.dat'
					swarp_xml_file=output_stack_swarp_dir+'/swarp_ss_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.xml'
					swarp_im_out=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.fits'
					swarp_weight_out=output_stack_swarp_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.WEIGHT.fits'
					sex_stack_cat_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.ldac'
					sex_stack_xml_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.xml'
;					sex_stack_checkimage_file=strjoin(output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+['.002.CHECK_BACK.fits','.002.CHECK_SEGMENTATION.fits'],',')
					sex_stack_checkimage_file=output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.002.CHECK_SEGMENTATION.fits'
				endif else $
				if do_sky_method EQ 'gradient tps' then begin
					sex_stack_im_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.fits'
					sex_stack_weight_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.WEIGHT.fits'
					sex_stack_cat_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_psf.003.ldac'
					sex_stack_xml_file=output_stack_sex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_psf.003.xml'
  	      sex_stack_cat_psfex_file=output_stack_sex_dir+'/ss_fornax_sex_psfex.ldac'
	        sex_stack_psfex_xml_file=output_stack_sex_dir+'/ss_fornax_sex_psfex.xml'
					sex_stack_checkimage_type='NONE'
					sex_stack_checkimage_file=output_stack_check_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.CHECK_SEGMENTATION.fits'
					sex_stack_magsize_file=output_stack_check_dir+'/sex_MAG_SIZE_ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.003.pdf'

        	psfex_stack_cat_file=output_stack_psfex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_stars.003.ldac'
        	psfex_stack_xml_file=output_stack_psfex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_stars.003.xml'
        	psfex_stack_psf_file=output_stack_psfex_dir+'/ss_fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'_stars.003.psf'
        	psfex_check_type = 'RESIDUALS,SNAPSHOTS,SAMPLES'
       		psfex_check_file = strjoin(output_stack_check_dir+'/psfex_'+['CHECK_RESIDUALS.fits','CHECK_SNAPSHOTS.fits','CHECK_SAMPLES.fits'],',')
        	psfex_checkplot_type = 'FWHM,ELLIPTICITY,COUNTS'
        	psfex_checkplot_file = strjoin(output_stack_check_dir+'/psfex_'+['CHECK_FWHM','CHECK_ELLIPTICITY','CHECK_COUNTS'],',')
				endif else begin
					swarp_list_file=output_stack_swarp_dir+'/swarp_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.dat'
					swarp_xml_file=output_stack_swarp_dir+'/swarp_fornax_'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.xml'
					swarp_im_out=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.fits'
					swarp_weight_out=output_stack_swarp_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.WEIGHT.fits'
					sex_stack_cat_file=output_stack_sex_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.ldac'
					sex_stack_xml_file=output_stack_sex_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.xml'
;					sex_stack_checkimage_file=strjoin(output_stack_check_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+['.CHECK_BACK.fits','.CHECK_SEGMENTATION.fits'],',')
					sex_stack_checkimage_file=output_stack_check_dir+'/fornax_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_'+type_uniq[k]+do_align+'.CHECK_SEGMENTATION.fits'
				endelse		

        case filter_uniq[j] of
          'u': begin
            sex_mag_range=[17.,21.]
            psfex_mag_range=[18.,20.]
            plot_mag_range=[24,14]
          end
          'g': begin
            sex_mag_range=[17.,21.]
            psfex_mag_range=[18.,20.];[17.,19.5]
            plot_mag_range=[25,15]
          end
          'i': begin
            sex_mag_range=[17.5,21.]
            psfex_mag_range=[18.,20.];[17.3,19.5]
            plot_mag_range=[25,15]
          end
          'z': begin
            sex_mag_range=[17.,21.]
            psfex_mag_range=[17.,19.]
            plot_mag_range=[25,15]
          end
          else: stop
        endcase

        print, sex_stack_im_file, file_test(sex_stack_im_file)
        print, psfex_stack_cat_file, file_test(psfex_stack_cat_file)
        print

        if file_test(sex_stack_im_file) AND (file_test(psfex_stack_cat_file) EQ 0 OR do_overwrite) then begin

          im_h=headfits(sex_stack_im_file)
          im_size=long([fxpar(im_h,'NAXIS1'),fxpar(im_h,'NAXIS2')])

;          if file_test(sex_stack_cat_psfex_file) EQ 0 OR do_overwrite then begin
          command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_psfex.sex -PARAMETERS_NAME sex_config/ctio_decam_psfex.param -CATALOG_NAME '+sex_stack_cat_psfex_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_psfex_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file+' -BACK_SIZE 384 -DETECT_THRESH 2. -ANALYSIS_THRESH 3. -PHOT_APERTURES '+ string(2*3*(1.2/0.2637),FORMAT='(F0.2)')
          print, command
          spawn, command
;          endif

          cat_data=mrdfits(sex_stack_cat_psfex_file, 2, cat_h, COLUMNS=['X_IMAGE','Y_IMAGE','ALPHA_J2000','DELTA_J2000','MAG_AUTO','FLUX_RADIUS','FLUX_APER','FLAGS','ELONGATION'])
          gv_stars=where(cat_data.flags LE 3 AND cat_data.mag_auto GE sex_mag_range[0] AND cat_data.mag_auto LT sex_mag_range[1], n_gv_stars)

					cgdelete, /all
					cgwindow, wxsize=800, wysize=600
          cgplot, cat_data.flux_radius, cat_data.mag_auto, xrange=[1,6], yrange=plot_mag_range, title='TILE'+tile_uniq[i]+'_FILTER'+filter_uniq[j], xtitle='Flux_radius (pixel)', ytitle='mag_auto', psym=cgsymcat('DOT'), symsize=0.6, /window
          cgplot, cat_data[gv_stars].flux_radius, cat_data[gv_stars].mag_auto, color=100, /over, /addcmd, psym=cgsymcat('DOT'), symsize=0.6
          cgplot, [0,100], sex_mag_range[0]*[1,1], line=2, color=100, /over, /addcmd
          cgplot, [0,100], sex_mag_range[1]*[1,1], line=2, color=100, /over, /addcmd

          plothist, cat_data[gv_stars].flux_radius, temp_xhist, temp_yhist, bin=sex_radius_bin, /noplot
          temp=max(temp_yhist, gv) & sex_radius=temp_xhist[gv]
          sex_radius=median((cat_data[where(cat_data.mag_auto GT sex_mag_range[0] AND cat_data.mag_auto LT sex_mag_range[1] AND cat_data.flux_radius GT sex_radius*0.9 AND cat_data.flux_radius LT sex_radius*1.1 , n_gv)]).flux_radius)

          cgplot, sex_radius*[1,1], [0,100], color=200, /over, /addcmd
          cgplot, sex_radius*[0.9,0.9], [0,100], line=2, color=200, /over, /addcmd
          cgplot, sex_radius*[1.1,1.1], [0,100], line=2, color=200, /over, /addcmd

          gv_stars=where(cat_data.mag_auto GT psfex_mag_range[0] AND cat_data.mag_auto LT psfex_mag_range[1] AND cat_data.flux_radius GT sex_radius*0.9 AND cat_data.flux_radius LT sex_radius*1.1 AND cat_data.flags LE sex_flag_max AND cat_data.x_image GT vig_diam/2. AND cat_data.x_image LT (im_size[0]-vig_diam/2.) AND cat_data.y_image GT vig_diam/2. AND cat_data.y_image LT (im_size[1]-vig_diam/2.), n_gv_stars)
          cgplot, [0,100], psfex_mag_range[0]*[1,1], line=2, color=200, /over, /addcmd
          cgplot, [0,100], psfex_mag_range[1]*[1,1], line=2, color=200, /over, /addcmd
          cgplot, cat_data[gv_stars].flux_radius, cat_data[gv_stars].mag_auto, color=200, /over, /addcmd, psym=cgsymcat('FILLEDCIRCLE'), symsize=0.5
					cgcontrol, output=sex_stack_magsize_file

          print, 'Number of selected stars ', n_gv_stars

          fits_open, sex_stack_cat_psfex_file, cat_fcb
          fits_read, cat_fcb, cat_data0, cat_h0, exten=0
          fits_read, cat_fcb, cat_data1, cat_h1, exten=1

          cat_data2=make_array(cat_fcb.axis[0:1,2], value=0, /byte)
          temp_lines=10000L
          temp_max=double(product(cat_fcb.axis[0:1,2])-1.)
          l=0L
          repeat begin
            print, 'Reading line ', strn(l*temp_lines), ' of ', strn(cat_fcb.axis[1,2])
            l1=double(l*temp_lines*cat_fcb.axis[0,2])
            l2=double((l+1)*temp_lines*cat_fcb.axis[0,2] - 1) < temp_max
            fits_read, cat_fcb, temp_data, cat_h2, first=l1, last=l2, exten_no=2
            cat_data2[*,l1/cat_fcb.axis[0,2]:(l2+1)/cat_fcb.axis[0,2]-1]=reform(temp_data, [cat_fcb.axis[0,2],(l2-l1+1)/cat_fcb.axis[0,2]])
            l++
          endrep until double(l*temp_lines*cat_fcb.axis[0,2]) GT temp_max
          temp_data=0
          fits_close, cat_fcb

          cat_data1=reform(cat_data1, [n_elements(cat_data1),1])
          cat_data2=cat_data2[*,[gv_stars]]
          fxaddpar, cat_h2, 'NAXIS2', (size(cat_data2, /dim))[1]

          if file_test(psfex_stack_cat_file, /regular, /noexpand) then  file_delete, psfex_stack_cat_file, /noexpand
          writefits, psfex_stack_cat_file, cat_data0, cat_h0
          fits_open, psfex_stack_cat_file, cat_fcb, /append
          fits_write, cat_fcb, cat_data1, cat_h1
          fits_write, cat_fcb, cat_data2, cat_h2
          fits_close, cat_fcb
          cat_data0=0 & cat_data1=0 & cat_data2=0

				endif

        if file_test(psfex_stack_cat_file) AND ((file_test(psfex_stack_psf_file) EQ 0 OR do_overwrite)) then begin
          command='psfex '+psfex_stack_cat_file+' -c psfex_config/ctio_decam.psfex'+' -PSF_DIR '+output_stack_psfex_dir+' -WRITE_XML Y -XML_NAME '+psfex_stack_xml_file+' -CHECKIMAGE_TYPE '+psfex_check_type+' -CHECKIMAGE_NAME '+psfex_check_file+' -CHECKPLOT_DEV PNG -CHECKPLOT_TYPE '+psfex_checkplot_type+' -CHECKPLOT_NAME '+psfex_checkplot_file+' -PSFVAR_NSNAP '+psfex_psfvar_nsnap+' -PSFVAR_DEGREES '+psfex_psfvar_degrees+' -PSFVAR_KEYS '+psfex_psfvar_keys+' -PSFVAR_GROUPS '+psfex_psfvar_groups+' -BASIS_TYPE '+psfex_basis_type+' -BASIS_NUMBER '+psfex_basis_number+' -SAMPLE_VARIABILITY 1.,1. -SAMPLE_MAXELLIP 1. -NEWBASIS_TYPE NONE -NEWBASIS_NUMBER 10 -SAMPLEVAR_TYPE NONE -STABILITY_TYPE EXPOSURE -SAMPLE_MINSN 1. -SAMPLE_FWHMRANGE 0.1,10. -SAMPLE_AUTOSELECT N -PSF_SIZE '+psfex_psf_size+' -PSF_SAMPLING '+psfex_psf_sampling+' -PSF_ACCURACY 0.01 -BADPIXEL_FILTER Y -BADPIXEL_NMAX 50'
          print, command
          spawn, command
        endif

        if file_test(psfex_stack_psf_file) AND ((file_test(sex_stack_cat_file) EQ 0 OR do_overwrite)) then begin
          command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack_psf.param -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file + ' -PSF_NAME '+psfex_stack_psf_file+' -PSF_NMAX 1'
;          print, command
;          spawn, command
        endif

      endfor
    endfor
  endfor

endif else $
if recipe EQ 'sex psf' then begin

	sex_stack_im_file=swarp_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.fits'
	sex_stack_weight_file=swarp_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.WEIGHT.fits'
	sex_stack_cat_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_psf.ldac'
	sex_stack_xml_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_psf.xml'
	sex_stack_checkimage_type='BACKGROUND,SEGMENTATION'
	sex_stack_checkimage_file=strjoin(sextractor_check_dir+'/fornax_tile'+do_tile+'_'+do_filter+['.CHECK_BACK.fits','.CHECK_SEGMENTATION.fits'],',')
	psfex_psf_file=psfex_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_psfex.psf'

	command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack_psf.param -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file + ' -PSF_NAME '+psfex_psf_file+' -PSF_NMAX 1'
	print, command
	spawn, command

endif else $
if recipe EQ 'cutout' then begin

	type_uniq=input_stack[uniq(input_stack.type, sort(input_stack.type))].type
	if file_test(output_stack_swarp_dir+'/cutout', /directory) EQ 0 then file_mkdir, output_stack_swarp_dir+'/cutout'

	im_file=file_search(output_stack_swarp_dir+'/*.fits')
	for i=0L, n_elements(type_uniq)-1 do begin

	for j=0L, n_elements(survey_cutout)-1 do begin

		gv=where(input_stack.type EQ type_uniq[i] AND input_stack.tile EQ (survey_cutout[j]).tile AND input_stack.filter EQ (survey_cutout[j]).filter, n_gv)
		if n_gv EQ 1 then begin
			print, 'CUTOUT - Processing image for galaxy ', (survey_cutout[j]).name
			
			im_h=headfits(input_stack[gv].im_file)
			wim_h=headfits(input_stack[gv].weight_file)

		  fits_open, input_stack[gv].im_file, fcb
 			im_data=make_array(fcb.axis[0:1], value=0., /float)
		  temp_lines=4000L
		  temp_max=long(product(fcb.axis[0:1])-1.)
		  k=0L
		  repeat begin
		  	print, 'Reading line ', strn(k*temp_lines), ' of ', strn(fcb.axis[1])
		    k1=long(k*temp_lines*fcb.axis[0])
		    k2=long((k+1)*temp_lines*fcb.axis[0] - 1) < temp_max
		    fits_read, fcb, temp_data, temp_h, first=k1, last=k2, exten_no=0
		    im_data[*,k1/fcb.axis[0]:(k2+1)/fcb.axis[0]-1]=reform(temp_data, [fcb.axis[0],(k2-k1+1)/fcb.axis[0]])
		    k++
		  endrep until long(k*temp_lines*fcb.axis[0]) GT temp_max
		  fits_close, fcb

		  fits_open, input_stack[gv].weight_file, fcb
 			wim_data=make_array(fcb.axis[0:1], value=0., /float)
		  temp_lines=4000L
		  temp_max=long(product(fcb.axis[0:1])-1.)
		  k=0L
		  repeat begin
		  	print, 'Reading line ', strn(k*temp_lines), ' of ', strn(fcb.axis[1])
		    k1=long(k*temp_lines*fcb.axis[0])
		    k2=long((k+1)*temp_lines*fcb.axis[0] - 1) < temp_max
		    fits_read, fcb, temp_data, temp_h, first=k1, last=k2, exten_no=0
		    wim_data[*,k1/fcb.axis[0]:(k2+1)/fcb.axis[0]-1]=reform(temp_data, [fcb.axis[0],(k2-k1+1)/fcb.axis[0]])
		    k++
		  endrep until long(k*temp_lines*fcb.axis[0]) GT temp_max
		  fits_close, fcb

			extast, im_h, im_ast

			ad2xy, (survey_cutout[j]).ra+4/60./cos((survey_cutout[j]).dec*!DTOR), (survey_cutout[j]).dec-4/60., im_ast, temp_x1, temp_y1
			ad2xy, (survey_cutout[j]).ra-4/60./cos((survey_cutout[j]).dec*!DTOR), (survey_cutout[j]).dec+4/60., im_ast, temp_x2, temp_y2
		  print, FORMAT='(%"Triming image --> X = [ %6i , %6i ]   Y = [ %6i , %6i ]")', (x1=(temp_x1<temp_x2>0L)), (x2=(temp_x1>temp_x2<(im_ast.naxis[0]-1))), (y1=(temp_y1<temp_y2>0)), (y2=(temp_y1>temp_y2<(im_ast.naxis[1]-1)))

		  hextract, im_data, im_h, nim_data, nim_h, x1, x2, y1, y2
		  hextract, wim_data, wim_h, nwim_data, nwim_h, x1, x2, y1, y2

		  writefits, output_stack_swarp_dir+'/cutout/'+(survey_cutout[j]).name+'_'+(survey_cutout[j]).filter+'_'+type_uniq[i]+'.fits', nim_data, nim_h
		  writefits, output_stack_swarp_dir+'/cutout/'+(survey_cutout[j]).name+'_'+(survey_cutout[j]).filter+'_'+type_uniq[i]+'.WEIGHT.fits', nwim_data, nwim_h

		endif else $
		if n_gv EQ 0 then begin
			print, 'There is no DECAm imaging for galaxy ', (survey_cutout[j]).name
		endif
		
	endfor
	endfor

endif else $
if recipe EQ 'sky variation' then begin

	plot_extname=['N15','N18','S15','S18']

	mjd_uniq=input_target[uniq(input_target.mjd_floor, sort(input_target.mjd_floor))].mjd_floor
	filter_uniq=input_target[uniq(input_target.filter, sort(input_target.filter))].filter
	type_uniq=input_target[uniq(input_target.type, sort(input_target.type))].type

	for i=0L, n_elements(do_program)-1 do begin
		log_file='survey_sky_'+do_program[i]+'.dat'
		if file_test(log_file) EQ 0 OR do_overwrite then begin
			openw, lun, log_file, /get_lun
			printf, lun, '#			im_file				mask_file			mjd			filter		  n_chip		chip  extname		sky_mean   sky_sigma'
			free_lun, lun
		endif
	endfor

	for i=0L, n_elements(mjd_uniq)-1 do begin
		print, 'Date ', date_conv(mjd_uniq[i]+2400000.5D, 'string')
		for j=0L, n_elements(filter_uniq)-1 do begin
			print, 'Filter ', filter_uniq[j]
			for k=0L, n_elements(type_uniq)-1 do begin
				print, 'Type ', type_uniq[k]
				
				gv_im=where(input_target.mjd_floor EQ mjd_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.type EQ type_uniq[k], n_gv_im)
				if n_gv_im EQ 0 then continue

				create_struct, im_stats, '', ['im_file','mask_file','mjd','filter','n_chip','chip','extname','sky_mean','sky_sigma'], 'A,A,D,A,I,I(62),A(62),F(62),F(62)', dim=n_gv_im
				
				for ii=0L, n_gv_im-1 do begin
					print, 'SKY SUBTRACTION - Processing image ', input_target[gv_im[ii]].im_file

					log_file='survey_sky_'+input_target[gv_im[ii]].program+'.dat'

					im_stats[ii].im_file=input_target[gv_im[ii]].im_file
					im_stats[ii].mask_file=input_target[gv_im[ii]].mask_file
					im_stats[ii].mjd=input_target[gv_im[ii]].mjd
					im_stats[ii].filter=input_target[gv_im[ii]].filter
					im_stats[ii].n_chip=input_target[gv_im[ii]].n_chip

					rdfits_struct, input_target[gv_im[ii]].im_file, im_struct, /silent
					rdfits_struct, input_target[gv_im[ii]].mask_file, mask_struct, /silent
					for jj=0L, input_target[gv_im[ii]].n_chip-1 do begin
						im_h=im_struct.(2*(jj+1))
						im_data=im_struct.(2*(jj+1)+1)
						mask_data=mask_struct.(2*(jj+1)+1)

						gv_sky=where(mask_data EQ 0, n_gv_sky)
						sky_mean=biweight_mean(im_data[gv_sky], sky_sigma)

						im_stats[ii].chip[jj]=jj+1
						im_stats[ii].extname[jj]=strtrim(fxpar(im_h, 'EXTNAME'),2)
						im_stats[ii].sky_mean[jj]=sky_mean
						im_stats[ii].sky_sigma[jj]=sky_sigma

					endfor

					n_chip=strn(im_stats[ii].n_chip-1)
					openw, lun, log_file, /get_lun, /append
					printf, lun, im_stats[ii].im_file, im_stats[ii].mask_file, im_stats[ii].mjd, im_stats[ii].filter, im_stats[ii].n_chip, im_stats[ii].chip[0:im_stats[ii].n_chip-1], im_stats[ii].extname[0:im_stats[ii].n_chip-1], im_stats[ii].sky_mean[0:im_stats[ii].n_chip-1], im_stats[ii].sky_sigma[0:im_stats[ii].n_chip-1], FORMAT='(A,4X,A,4X,F0.6,4X,A,4X,I,4X,'+n_chip+'(I0,","),I0,4X,'+n_chip+'(A,","),A,4X,'+n_chip+'(F0.2,","),F0.2,4X,'+n_chip+'(F0.2,","),F0.2)'
					free_lun, lun

				endfor

				!x.margin=[14,4]
				!y.margin=[5,4]
				plotsym, 0, 0.2, /fill 

				gv_chip=where(im_stats[0].extname EQ plot_extname[0], n_gv_chip)
				plot_xrange=[floor(min(im_stats.mjd)*1e2)/1d2, ceil(max(im_stats.mjd)*1e2)/1d2]
				plot_yrange=[floor(min(im_stats.sky_mean[gv_chip])/2e2)*2d2, ceil(max(im_stats.sky_mean[gv_chip])/2e2)*2d2]

				cgdelete, /all
				cgloadct, 0
				cgwindow, wxsize=800, wysize=1200;, WMulti=[0,1,n_elements(plot_extname)], WOXMargin=[2,6], WOYMargin=[2,2]
				plot_pos = cgLayout([1,n_elements(plot_extname)], OXMargin=[8,6], OYMargin=[2,2], XGap=1, YGap=4)

				for jj=0L, n_elements(plot_extname)-1 do begin
					gv_chip=where(im_stats[0].extname EQ plot_extname[jj], n_gv_chip)

					cgplot, [0], [0], position=plot_pos[*,jj], xrange=plot_xrange, yrange=plot_yrange, xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, xtickformat='(F0.3)', xtitle=(jj EQ n_elements(plot_extname)-1 ? 'MJD (days)':''), ytitle='Sky (ADU)', /nodata, /addcmd, /noerase, title='EXTNAME '+plot_extname[jj], charsize=1.
					cgplot, im_stats.mjd, transpose(im_stats.sky_mean[gv_chip]), position=plot_pos[*,jj] , psym=cgsymcat('filled circle'), symsize=0.8, color='red', /over, /addcmd
;					cgtext, plot_xrange[1]-(plot_xrange[1]-plot_xrange[0])/20, plot_yrange[1]-(plot_yrange[1]-plot_yrange[0])/10, 'EXTNAME '+plot_extname[jj], /data, align=1, charsize=1., /addcmd

				endfor

				cgcontrol, output='results/fornax_sky_level_'+string(mjd_uniq[i],FORMAT='(I0)')+'_'+filter_uniq[j]+'_'+type_uniq[k]+'.pdf'

			endfor
		endfor
	endfor
	stop

endif



end
