function hist_yubo, lst, GAUSS=gauss, PLOT=plot, SIGM=sigm
compile_opt idl2, HIDDEN
; saves some time, just inputs some reasonable params for hist_wrapper
; uses bins of width  [[ stddev(list) / n_elements(list)^(1/3) ]]
list = lst
if keyword_set(sigm) then begin
    list = list[where(list lt mean(list) + 5 * sigm)]
    list = list[where(list gt mean(list) - 5 * sigm)]
endif else begin
    list = list[where(list lt mean(list) + 5 * stddev(list))]
    list = list[where(list gt mean(list) - 5 * stddev(list))]
endelse

if keyword_set(PLOT) then begin
    if keyword_set(GAUSS) then begin
        RETURN, hist_wrapper(list, stddev(list)/ n_elements(list)^(1.0/3), min(list), max(list), /plot, /gauss, /ylog) 
    endif else begin
        return, hist_wrapper(list, stddev(list)/ n_elements(list)^(1.0/3), min(list), max(list), /plot)
    endelse
endif else begin
    if keyword_set(GAUSS) then begin
        RETURN, hist_wrapper(list, stddev(list)/ n_elements(list)^(1.0/3), min(list), max(list), /gauss)
    endif else begin
        return, hist_wrapper(list, stddev(list)/ n_elements(list)^(1.0/3), min(list), max(list))
    endelse
endelse

end
