readcol, 'decam_footprint.dat', chip, ra_bottom, dec_bottom, ra_top, dec_top, FORMAT='I,A,A,A,A'

im_file='/Users/rmunoz/Research/NGFS/DECam/data_reduction/pipeline/swarp/fornax_tile1_i.fits'
im_h=headfits(im_file)
extast, im_h, im_ast
ra=im_ast.CRVAL[0]
dec=im_ast.CRVAL[1]

ra_bottom=tenv(ra_bottom)*360./24
dec_bottom=tenv(dec_bottom)
ra_top=tenv(ra_top)*360./24
dec_top=tenv(dec_top)

ad2xy, ra, dec, im_ast, x_center, y_center
ad2xy, ra_bottom, dec_bottom, im_ast, x_bottom, y_bottom
ad2xy, ra_top, dec_top, im_ast, x_top, y_top

delta_ra_bottom=(x_bottom-x_center)*im_ast.CD[0,0]
delta_dec_bottom=(y_bottom-y_center)*im_ast.CD[1,1]
delta_ra_top=(x_top-x_center)*im_ast.CD[0,0]
delta_dec_top=(y_top-y_center)*im_ast.CD[1,1]

;delta_ra_bottom=(ra_bottom-ra)*cos(dec_bottom*!DTOR)
;delta_dec_bottom=dec_bottom-dec
;delta_ra_top=(ra_top-ra)*cos(dec_top*!DTOR)
;delta_dec_top=dec_top-dec

openw, lun, 'aladin_CTIO_DECam.xml', /get_lun
printf, lun, '<?xml version="1.0" encoding="UTF-8"?'
printf, lun, '<VOTABLE xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" version="1.1" xmlns="http://www.ivoa.net/xml/VOTable/v1.1" xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.1 http://www.ivoa.net/xml/VOTable/v1.1">'
printf, lun, '<DESCRIPTION>Created by Roberto Pablo Munoz, Pontificia Universidad Catolica de Chile</DESCRIPTION>'
printf, lun, ''
printf, lun, '<RESOURCE ID="Blanco_DECam" utype="dal:footprint.geom">'
printf, lun, '    <PARAM datatype="char" arraysize="*" ID="TelescopeName" value="Blanco" />'
printf, lun, '    <PARAM datatype="char" arraysize="*" ID="InstrumentName" value="DECam" />'
printf, lun, '    <PARAM datatype="char" arraysize="*" ID="InstrumentDescription" value="Dark Energy Camera FOV" />'
printf, lun, '    <PARAM datatype="char" arraysize="*" ID="Origin" value="Roberto Pablo Munoz" />'
printf, lun, ''
printf, lun, '    <PARAM datatype="char" arraysize="*" utype="stc:AstroCoordSystem.CoordFrame.CARTESIAN" name="reference frame" value="*"/>'
printf, lun, '    <PARAM name="projection" utype="stc:AstroCoordSystem.CoordFrame.Cart2DRefFrame.projection" datatype="char" arraysize="*" value="TAN"/>'
printf, lun, ''

for i=0L, n_elements(chip)-1 do begin
	printf, lun, ''
	printf, lun, '    <RESOURCE ID="shape_'+strn(i+1)+'" name="chip_'+strn(chip[i])+'">'
	printf, lun, '    <TABLE utype="dal:footprint.geom.segment">'
	printf, lun, '    <PARAM datatype="char" name="Shape" arraysize="*" value="Polygon" utype="dal:footprint.geom.segment.shape" />'
	printf, lun, '    <FIELD unit="arcsec" datatype="double" name="xPtPosition" utype="stc:AstroCoordArea.Polygon.Vertex.Position.C1" />'
	printf, lun, '    <FIELD unit="arcsec" datatype="double" name="yPtPosition" utype="stc:AstroCoordArea.Polygon.Vertex.Position.C2" />'
	printf, lun, '    <DATA>'
	printf, lun, '    <TABLEDATA>'
	printf, lun, delta_ra_bottom[i]*3600., delta_dec_bottom[i]*3600., FORMAT='("        <TR><TD>",F0.1,"</TD><TD>",F0.1,"</TD></TR>")'
	printf, lun, delta_ra_top[i]*3600., delta_dec_bottom[i]*3600.,  FORMAT='("        <TR><TD>",F0.1,"</TD><TD>",F0.1,"</TD></TR>")'
	printf, lun, delta_ra_top[i]*3600., delta_dec_top[i]*3600.,  FORMAT='("        <TR><TD>",F0.1,"</TD><TD>",F0.1,"</TD></TR>")'
	printf, lun, delta_ra_bottom[i]*3600., delta_dec_top[i]*3600.,  FORMAT='("        <TR><TD>",F0.1,"</TD><TD>",F0.1,"</TD></TR>")'
	printf, lun, '    </TABLEDATA>'
	printf, lun, '    </DATA>'
	printf, lun, '    </TABLE>'
	printf, lun, '    </RESOURCE>'
endfor

printf, lun, '</RESOURCE>'
printf, lun, '</VOTABLE>'

free_lun, lun

end
