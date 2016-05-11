forward_function decam_data

pro decam_pipeline, recipe, PROGRAM=program, FILTER=filter, OVERWRITE=overwrite, TILE=tile, TYPE=type, DEBUG=debug, STANDARD=standard, AHEAD=ahead, SKY_METHOD=sky_method, SKY_NIM=sky_nim, SKY_NCHIP=sky_nchip, SKY_USE_TARGET=sky_use_target

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
do_sky_nchip= (n_elements(sky_nchip) GT 0) ? sky_nchip : 31
do_sky_nim= (n_elements(sky_nim) GT 0) ? sky_nim : 5
do_sky_method = (n_elements(sky_method) GT 0) ? sky_method : ''
do_sky_tile_exclude =  (n_elements(sky_tile_exclude) GT 0) ? sky_tile_exclude : '1'
do_sky_use_target = (n_elements(sky_use_target) GT 0) ? keyword_set(sky_use_target) : 0
do_sky_region = (n_elements(sky_region) GT 0) ? sky_region : [[0,0,1022,4093],[1023,0,2045,4093]]

input_dir='/Volumes/Q6/NGFS/DECam'
output_dir='/Volumes/Q6/NGFS/DECam'

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
output_stack_check_dir=output_dir+'/stacks/check';output_dir+'/'+do_program+'/pipeline/swarp'
output_stack_psfex_dir=output_dir+'/stacks' ;output_dir+'/'+do_program+'/pipeline/psfex'

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
		gv=where( total(rebin(byte(input_calib.filter),[2,n_elements(input_calib.filter)]) EQ rebin(transpose(byte(do_filter_split)),[2,n_elements(input_calib.filter)]),1) GT 0, n_gv)
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

