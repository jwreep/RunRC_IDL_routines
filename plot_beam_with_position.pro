pro plot_beam_with_position,HD_variables,beam,time_index,length=length,$
    filename=filename,title=title,xlinear=xlinear,cadence=cadence,step=step


  if not keyword_set(length) then length = get_loop_length()
  if not keyword_set(filename) then filename = 'beam_with_position.eps'
  if not keyword_set(cadence) then cadence = 1.0
  time = cadence * (time_index - HD_variables.start_time) / HD_variables.time_step
  if not keyword_set(title) then title = 'Time = '+strtrim(string(time,format='(F6.2)'),1)+' s'
  if not keyword_set(step) then step=5

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
    plot,beam.energy[*,apex_index,time_index],beam.flux[*,s_index,time_index], $
      xrange=[0,max(beam.energy)],xstyle=1,xtitle='Electron Energy [keV]', $
      yrange=[max(beam.flux)/1e5,max(beam.flux)],/ylog,ytitle='Energy Flux [erg s!E-1!N cm!E-2!N keV!E-1!N]', $
      charsize=2.,charth=5,thick=10,/nodata,title=title
  endif else begin
    plot,beam.energy[*,apex_index,time_index],beam.flux[*,s_index,time_index], $
      xrange=[1.0,max(beam.energy)],/xlog,xtitle='Electron Energy [keV]', $
      yrange=[max(beam.flux)/1e5,max(beam.flux)],/ylog,ytitle='Energy Flux [erg s!E-1!N cm!E-2!N keV!E-1!N]', $
      charsize=2.,charth=5,thick=10,/nodata,title=title
  endelse

  if n_elements(s_index) le 0 then begin
    distinct_colors,n_colors=n_elements(s_index)
    for i=0,n_elements(s_index)-1 do begin
      oplot,beam.energy[*,s_index[i],time_index],beam.flux[*,s_index[i],time_index], $
        linestyle=0,thick=6,color=i
    endfor
  endif else begin
    loadcv,80   ;; the first 13 elements of this table have dummy colors
    col_increment = (255-13.)/n_elements(s_index)
    for i=0,n_elements(s_index)-1 do begin
      oplot,beam.energy[*,s_index[i],time_index],beam.flux[*,s_index[i],time_index], $
        linestyle=0,thick=6,color=13+i*col_increment
    endfor
  endelse


  device,/close_file



end