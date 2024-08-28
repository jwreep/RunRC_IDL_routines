function read_bm_files,start_time,end_time,time_step,N_NT_ENERGY=N_NT_ENERGY
;; Reads in the kinetic non-thermal beam files ('.bm') which contain the
;; beam parameters as they propagate down the loop

  if not keyword_set(N_NT_ENERGY) then N_NT_ENERGY=200

  maximum_number_of_rows = find_maximum_number_of_rows(start_time,end_time,time_step)

  input_array = dblarr(N_NT_ENERGY, /nozero)
  energy_array = make_array(N_NT_ENERGY, maximum_number_of_rows, $
                            (end_time-start_time)/time_step+1,/double,value=!Values.F_NAN)
  flux_array = make_array(N_NT_ENERGY, maximum_number_of_rows, $
                            (end_time-start_time)/time_step+1,/double,value=!Values.F_NAN)

  for time=start_time,end_time,time_step do begin

    filen='profile'+strtrim(string(time),1)+'.bm'
    openr,lun,filen,/get_lun                ;open
    count=0

    while( not eof(lun)) do begin

      readf,lun,pos
      readf,lun,input_array
      energy_array[*,count,(time-start_time)/time_step] = input_array / 1.602e-9
      readf,lun,input_array
      flux_array[*,count,(time-start_time)/time_step] = input_array

      count = count+1

    endwhile

    close,/all
    free_lun,lun

  endfor

  beam = create_struct( 'energy',energy_array, 'flux',flux_array, $
                          'energy_units','keV','flux_units','erg/s/cm^2')

  return, beam


end