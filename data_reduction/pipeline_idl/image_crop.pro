im_file='/Volumes/Q6/rmunoz/Research/NGFS/DECam/'+['2014B-0609/data_reduction/pipeline/swarp/fornax_tile1_u.fits','2013B-0613/data_reduction/pipeline/swarp/fornax_tile1_g.fits','2013B-0613/data_reduction/pipeline/swarp/fornax_tile1_i.fits']
nim_file='/Volumes/Q6/rmunoz/Research/NGFS/DECam/color_images/'+['fornax_tile1_u_crop.fits','fornax_tile1_g_crop.fits','fornax_tile1_i_crop.fits']

im_coo= { ra:tenv(['03:40:53.46','03:35:51.78'])*360./24, dec:tenv(['-35:56:58.4','-34:54:56.6']) }

for i=0L, n_elements(im_file)-1 do begin

	print, 'Processing image ', im_file[i]

	im_h=headfits(im_file[i])
  fits_open, im_file[i], fcb
  im_data=make_array(fcb.axis[0:1], value=0., /float)
  temp_lines=4000L
  temp_max=long(product(fcb.axis[0:1])-1.)
  j=0L
  repeat begin
  	print, 'Reading line ', strn(j*temp_lines), ' of ', strn(fcb.axis[1])
    j1=long(j*temp_lines*fcb.axis[0])
    j2=long((j+1)*temp_lines*fcb.axis[0] - 1) < temp_max
    fits_read, fcb, temp_data, temp_h, first=j1, last=j2, exten_no=0
    im_data[*,j1/fcb.axis[0]:(j2+1)/fcb.axis[0]-1]=reform(temp_data, [fcb.axis[0],(j2-j1+1)/fcb.axis[0]])
    j++
  endrep until long(j*temp_lines*fcb.axis[0]) GT temp_max
  fits_close, fcb

	extast, im_h, im_ast

	ad2xy, im_coo.ra[0], im_coo.dec[0], im_ast, temp_x1, temp_y1
	ad2xy, im_coo.ra[1], im_coo.dec[1], im_ast, temp_x2, temp_y2
  print, FORMAT='(%"Triming image --> X = [ %6i , %6i ]   Y = [ %6i , %6i ]")', (x1=(temp_x1<temp_x2>0L)), (x2=(temp_x1>temp_x2<(im_ast.naxis[0]-1))), (y1=(temp_y1<temp_y2>0)), (y2=(temp_y1>temp_y2<(im_ast.naxis[1]-1)))
  hextract, im_data, im_h, nim_data, nim_h, x1, x2, y1, y2
  writefits, nim_file[i], nim_data, nim_h

endfor

end
