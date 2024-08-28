pro plot_kinetic_ntb,HD_variables,beam,time_index,length=length,energy=energy,$
  filename=filename,title=title,xlinear=xlinear,cadence=cadence,step=step


  if not keyword_set(length) then length = get_loop_length()
  if not keyword_set(filename) then filename = 'hxr_with_position.eps'
  if not keyword_set(cadence) then cadence = 1.0
  time = cadence * (time_index - HD_variables.start_time) / HD_variables.time_step
  if not keyword_set(title) then title = 'Time = '+strtrim(string(time,format='(F6.2)'),1)+' s'
  if not keyword_set(step) then step=5
  if not keyword_set(energy) then energy = findgen(2001,start=0.,increment=0.1)

  ntb = kinetic_nonthermal_emission(time_index,time_index,1)
  total_ntb = dblarr(n_elements(energy))
  for i=0,n_elements(ntb[*,0,0])-1 do total_ntb += ntb[i,0,*]>0

  apex_index = locate_apex(HD_variables,time_index)

  for i=0,apex_index,step do begin
    if beam.energy[-1,i,time_index] ge 1.0 then break
  endfor
  bottom_index = i

  s_index = reverse(indgen( apex_index/step, start=bottom_index, increment=step ))

  !x.thick=6
  !y.thick=6
  !P.thick=6
  !p.multi = 0

  distinct_colors,n_colors=1

  set_plot,'PS'
  device,filename=filename,/COLOR, BITS=8, xsize=10, ysize=6,/encap,/inch

  if keyword_set(xlinear) then begin
    plot,energy,total_ntb,$
      xrange=[0,max(beam.energy)],xstyle=1,xtitle='Photon Energy [keV]', $
      yrange=[max(total_ntb)/1e10, max(total_ntb)],/ylog,ytitle='Photon Flux [photons s!E-1!N cm!E-2!N keV!E-1!N]', $
      charsize=2., charth=5, thick=10, title=title
  endif else begin
    plot,energy,total_ntb,$
      xrange=[1.0,max(beam.energy)],/xlog,xtitle='Photon Energy [keV]', $
      yrange=[max(total_ntb)/1e10, max(total_ntb)],/ylog,ytitle='Photon Flux [photons s!E-1!N cm!E-2!N keV!E-1!N]', $
      charsize=2., charth=5, thick=10, title=title
  endelse

  if n_elements(s_index) le 0 then begin
    distinct_colors,n_colors=n_elements(s_index)
    for i=0,n_elements(s_index)-1 do begin
      oplot, energy, ntb[s_index[i], 0, *], $
        linestyle=0, thick=6, color=i
    endfor
  endif else begin
    loadcv,80   ;; the first 13 elements of this table have dummy colors
    col_increment = (255-13.)/n_elements(s_index)
    for i=0,n_elements(s_index)-1 do begin
      oplot, energy, ntb[s_index[i], 0, *], $
        linestyle=0, thick=6, color=13+i*col_increment
    endfor
  endelse


  device,/close_file



end