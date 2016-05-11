pro callback_convolve, status, error, bridge, ud
	out = bridge->getvar('out')

	(*(ud.pout))[ud.xrange[0]:ud.xrange[1]-1,ud.yrange[0]:ud.yrange[1]-1] = out
end
