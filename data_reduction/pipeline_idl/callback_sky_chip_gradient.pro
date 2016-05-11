pro callback_sky_median_global, status, error, bridge, ud
	out_im = bridge->getvar('out_im')
	out_sky = bridge->getvar('out_sky')

	(*(ud.pout_im))[*,*,ud.chip_offset+ud.chip] = out_im
	(*(ud.pout_sky))[*,*,ud.chip_offset+ud.chip] = out_sky
end
