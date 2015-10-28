function addfunc_multi, in_sig, in_a, in_sigm, num_beams, MINSEP=minsep, emissivity=emissivity, tdust=tdust
compile_opt idl2, HIDDEN
    ; Input
    ;   in_sig       - Existing map (usual struct with signal, freqs)
    ;   in_a         - height of Gaussians to be added
    ;   in_sigm      - sigma of gaussians to be added
    ;   num_beams    - num of Gaussians to add
    ;
    ; Output - struct
    ;   signal       - Usual struct containing noise + added Gaussians
    ;   x0           - x coords of peaks
    ;   y0           - y coords of peaks
    ;
    ; wrapper around addGauss to add random Gaussians


edgeRange = 3 * in_sigm        ; do not make Gaussians too close to edge
if keyword_set(minsep) then minDist = minsep else minDist = 0 * in_sigm
    ; do not make Gaussians too close together
range = double(sqrt(n_elements(in_sig.signal[*,*,0])))


; set up some stuff to compute Gaussians
x = RANDOMU(seed) * (range - 2 * edgeRange - 1) + edgeRange
y = RANDOMU(seed) * (range - 2 * edgeRange - 1) + edgeRange 
x0 = [x] ; store generated coords
y0 = [y]
signal = addgauss_multi(in_a, in_sigm, x, y, in_sig)

for i=2, num_beams do begin
    repeat begin ; generate new coords at least once
        x = RANDOMU(seed) * (range - 2 * edgeRange - 1) + edgeRange
        y = RANDOMU(seed) * (range - 2 * edgeRange - 1) + edgeRange 
    endrep until min((x0 - x)^2 + (y0 - y)^2) gt minDist^2 ; until mindist exceeded
    x0 = [x0,x]
    y0 = [y0,y]
    if keyword_set(tdust) then begin
        if keyword_set(emissivity) then signal = addgauss_multi(in_a, in_sigm, x, y, signal, emissivity=emissivity, tdust=tdust) else signal = addgauss_multi(in_a, in_sigm, x, y, signal, tdust=tdust)
    endif else begin
        if keyword_set(emissivity) then signal = addgauss_multi(in_a, in_sigm, x, y, signal, emissivity=emissivity) else signal = addgauss_multi(in_a, in_sigm, x, y, signal)
    endelse
endfor

return, {sig:signal, x0:x0, y0:y0}
end
