pro callback_sky_surface, status, error, bridge, ud
	out_im = bridge->getvar('out_im')
	out_sky = bridge->getvar('out_sky')

	ud=*ud

	if size(*(ud.pout_im), /n_dim) EQ 3 then (*(ud.pout_im))[*,*,ud.chip_offset+ud.chip] = out_im
	if size(*(ud.pout_sky), /n_dim) EQ 3 then (*(ud.pout_sky))[*,*,ud.chip_offset+ud.chip] = out_sky
	
end
