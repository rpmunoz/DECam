pro callback_sky_surface, status, error, bridge, ud
	out_im = bridge->getvar('out_im')
	out_im_masked = bridge->getvar('out_im_masked')
;	out_im_masked_median = bridge->getvar('out_im_masked_median')
	out_im_masked_res = bridge->getvar('out_im_masked_res')
	out_sky = bridge->getvar('out_sky')
;	plot_x = bridge->getvar('plot_x')
;	plot_y = bridge->getvar('plot_y')
;	plot_xrange = bridge->getvar('plot_xrange')
;	plot_yrange = bridge->getvar('plot_yrange')
;	plot_pos = bridge->getvar('plot_pos')

	ud=*ud

	if size(*(ud.pout_im), /n_dim) EQ 3 then (*(ud.pout_im))[*,*,ud.chip_offset+ud.chip] = out_im
	if size(*(ud.pout_sky), /n_dim) EQ 3 then (*(ud.pout_sky))[*,*,ud.chip_offset+ud.chip] = out_sky
;	if size(*(ud.pout_im_masked), /n_dim) EQ 3 then (*(ud.pout_im_masked))[*,*,ud.chip_offset+ud.chip] = out_im_masked
;	if size(*(ud.pout_im_masked_median), /n_dim) EQ 3 then (*(ud.pout_im_masked_median))[*,*,ud.chip_offset+ud.chip] = out_im_masked_median
;	if size(*(ud.pout_im_masked_res), /n_dim) EQ 3 then (*(ud.pout_im_masked_res))[*,*,ud.chip_offset+ud.chip] = out_im_masked_res

;	if ud.plot_do then begin
;		(ud.plot_window).select
;		for i=0, n_elements(plot_x)-1 do begin
;			print, plot_x.p1
;			print, plot_y.p1
;;			p=plot(plot_x[i], plot_y[i], position=plot_pos[i], xrange=plot_xrange[i], yrange=plot_yrange[i], xcharsize=0.8, ycharsize=0.8, xstyle=1, ystyle=1, /current)
;		endfor
;	endif
	
end
