;; Synthesizes the nonthermal emission for flare simulations
;; photons sec^-1 cm^-2 kev^-1


;; Calculate the relativistic Bethe-Heitler cross-section with Elwert correction factor:
function cross_section,photon_energy, electron_energy

if electron_energy le photon_energy or electron_energy le 0.0 or photon_energy le 0.0 then begin
  return,0.0
endif else begin

  m_e = 510.998950d   ;; electron mass in keV
  c1 = 1.1289762657670394d-21 ;; = 16 <Z^2> r_0^2 alpha m_e^2 / 3
  c2 = 0.045848400982856752 ;; 2 pi alpha (alpha = fine structure constant)
  
  pe = photon_energy
  ee = electron_energy
  
  ;; Relativistic Bethe-Heitler cross-section with Elwert correction factor:
  cs = c1 / (pe * ee * (ee + 2*m_e)) * $
      alog( (1.0 + sqrt(((ee - pe)*(ee - pe + 2*m_e))/(ee * (ee + 2*m_e)) ) ) $
      / (1.0 - sqrt(((ee-pe)*(ee-pe+2*m_e))/(ee*(ee+2*m_e))) ) ) $
      * sqrt(ee*(ee+2*m_e))*(ee-pe+m_e) / (sqrt((ee-pe)*(ee-pe+2*m_e))*(ee+m_e)) $
      * (1.0 - exp(-c2*(ee+m_e)/sqrt(ee*(ee+2*m_e)))) $
      / (1.0 - exp(-c2*(ee-pe+m_e)/sqrt((ee-pe)*(ee-pe+2*m_e))) )

  return,cs 

endelse

end



function kinetic_nonthermal_emission,start_time,end_time,time_step,$
  area=area,use_hstate=use_hstate,energy=energy

  h = 6.626070150d-27   ; erg second
  h_kev = 4.135668d-18   ; keV second
  c = 2.99792458d18  ; Angstrom per second
  pi = 3.14159265358979d
  r_sun = 6.957d10 ; cm
  AU = 1.496d13 ; cm

  if not keyword_set(area) then area = 1.0d17  ;; cross-sectional area of the loop

  ;; Read in simulation data:
  variables = read_phy_files(start_time,end_time,time_step)
  beam = read_bm_files(start_time,end_time,time_step)

  ;; Get the H II fractions
  if keyword_set(use_hstate) then begin
    Hstate = read_Hstate_files(start_time,end_time,time_step)
  endif else begin
    ; Use CHIANTI's ionization equilibrium values if not read in
    restore,'/Users/reep/Desktop/IDL/HYDRAD_IDL_Routines/chianti_ioneq.sav'  ;; v10.1 CHIANTI ionization equilibria
    hstate_temp = 10.^(findgen(101,start=4.,increment=0.01))
    hstate = interpol(ioneq[*,0,1],T,alog10(hstate_temp) )
    ;hstate = get_ieq(hstate_temp,'h_2')   ;; Get the LTE fraction of H II from Chianti
  endelse

  loop_length = get_loop_length()

  ;; Set up arrays:
  time = findgen(1+(end_time-start_time)/time_step,start=start_time,increment=time_step)

  if not keyword_set(energy) then energy = findgen(2001,start=0.,increment=0.1)
  
  thick_brems = make_array(variables.num_columns, n_elements(time), n_elements(energy), value=0.0d)
  
  total_brems = make_array(n_elements(time), n_elements(energy), value=0.0d)

  ;; Constants:
  m_e = 510.998950d   ;; electron mass in keV
  ;; Average energy of an electron in the distribution:
  ;avg_energy = cutoff * (delta - 1.0)/(delta - 2.0)
  ;L2 = 25.1 + alog(avg_energy)


  for t=0,n_elements(time)-1 do begin

    ;; Loop over the grid cells to calculate the non-thermal bremsstrahlung in each cell
    for n = 0,variables.num_columns-1 do begin

      if not finite(variables.position[n,t]) then break

      if keyword_set(use_hstate) then begin
        ion_frac = Hstate.ionized[n,t]
      endif else begin   ;; Use LTE if we don't have the Hstate files from HYDRAD
        if variables.electron_temperature[n,t] gt 10.0^4.9 then begin
          ion_frac = 1.0
        endif else if variables.electron_temperature[n,t] lt 10.0^4.0 then begin
          ion_frac = 0.0
        endif else begin
          ion_frac = interpol(hstate, hstate_temp, variables.electron_temperature[n,t])
        endelse
      endelse

      ;L1 = 66. + 1.5*alog(avg_energy) - 0.5 *alog(variables.electron_density[n,t])
      L1 = 20.0

      ;; No emission if the beam is thermalized:
      if min(beam.energy[*,n,t]) eq max(beam.energy[*,n,t]) then continue
      
      dE = beam.energy[-1,n,t] - beam.energy[-2,n,t]
      count = 0
      emin = beam.energy[0,n,t]
      while( emin-dE gt 0.0 ) do begin
        count += 1
        emin -= dE
      endwhile
      if count gt 0 then begin
        energy_array = [dindgen(count, start=emin, increment=dE), beam.energy[*,n,t]]
        flux_array = [dblarr(count), beam.flux[*,n,t]]/1.602e-9
                                                      ;; convert from erg to keV
      endif else begin
        energy_array = beam.energy[*,n,t]
        flux_array = beam.flux[*,n,t]/1.60e2-9
                                      ;; convert from erg to keV
      endelse
      brems_array = make_array(n_elements(energy_array), value=0.0d)
      
      for k=0,n_elements(energy_array)-1 do begin
        
        n_i = 20  ;; number of integration points for inner integral -- must be even
        n_o = 20  ;; number of integration points for outer integral -- must be even
        if n_i mod 2 ne 0 then n_i += 1
        if n_o mod 2 ne 0 then n_o += 1
       
        d = max(energy_array) ;; upper limit of the integral
        c = energy_array[k] ;; lower limit of the integral

        if c ge d then break

        prefactor = 2.0*(d-c)/(3.0*n_o)
        outer_sum1 = 0.0
        outer_sum2 = 0.0
        
        for j=1,n_o/2-1 do begin
          electron_energy = c + 2.0*j*(d-c)/n_o
          if electron_energy gt d then break
          
          cs = cross_section(c, electron_energy)
          inner_sum1 = 0.0
          inner_sum2 = 0.0
          
          dummy = min(abs(energy_array - electron_energy), index1)
          
          for i=1,n_i/2-1 do begin
            ee = electron_energy + 2.0*i*(d - electron_energy)/n_i
            dummy = min(abs(energy_array - ee), index2)
            inner_sum1 += flux_array[index2]
          endfor
          
          for i=1,n_i/2 do begin
            ee = electron_energy + (2.0*i-1.0)*(d - electron_energy)/n_i
            dummy = min(abs(energy_array - ee), index3)
            inner_sum2 += flux_array[index3]
          endfor
            
          outer_sum1 += electron_energy * cs * (flux_array[-1] + flux_array[index1] + 2.0*inner_sum1 + 4.0 * inner_sum2)
          
        endfor
        
        for j=1,n_o/2 do begin
          electron_energy = c + (2.0*j-1.0)*(d-c)/n_o
          if electron_energy gt d then break
          
          cs = cross_section(c, electron_energy)
          inner_sum1 = 0.0
          inner_sum2 = 0.0

          dummy = min(abs(energy_array - electron_energy), index1)
          
          for i=1,n_i/2-1 do begin
            ee = electron_energy + 2.0*i*(d - electron_energy)/n_i
            dummy = min(abs(energy_array - ee), index2)
            inner_sum1 += flux_array[index2]
          endfor

          for i=1,n_i/2 do begin
            ee = electron_energy + (2.0*i-1.0)*(d - electron_energy)/n_i
            dummy = min(abs(energy_array - ee), index3)
            inner_sum2 += flux_array[index3]
          endfor

          outer_sum2 += electron_energy * cs * (flux_array[-1] + flux_array[index1] + 2.0*inner_sum1 + 4.0 * inner_sum2)

        endfor

        brems_array[k] = prefactor * (outer_sum1 + 2.0 * outer_sum2)
        
      endfor

      brems_array *= 2.7294774653513575d-9 * area / (ion_frac * L1)
      ;brems_array *= 2.7294774653513575d-9 * area / (ion_frac * L1 + (1.0d - ion_frac) * L2)
      ;; 2.7294774653513575d-9 = 1 / (8 pi^2 AU^2 e^4 (1.602d-9)^2)
      ;; 1.602d-9 erg to keV conversion

      thick_brems[n,t,*] = interpol(brems_array, energy_array, energy)

    endfor ;; position loop

  endfor ;;time loop


;  for t=0,n_elements(time)-1 do begin
;    for i=0,n_elements(energy)-1 do begin
;      total_brems[t,i] = total( thick_brems[*,t,i], /nan )
;      if keyword_set(ntr) then total_brems[t,i] += total_recomb[t,i]
;    endfor
;  endfor

  ;return, total_brems   ;; photon/s/cm^2/keV
  return,thick_brems

end