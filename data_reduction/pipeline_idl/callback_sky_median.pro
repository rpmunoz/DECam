pro callback_sky_median, status, error, bridge, ud
	out = bridge->getvar('out')

	(*(ud.pout))[*,*,ud.chip_offset+ud.chip] = out
end
