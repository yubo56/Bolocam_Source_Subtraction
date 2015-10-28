function addfunc, noise, gaussHeight, sigm, numGaussians, MINSEP=minsep
compile_opt idl2, HIDDEN
    ; Input
    ;   noise        - Existing noise map
    ;   gaussHeight  - height of Gaussians to be added
    ;   sigm         - sigma of gaussians to be added
    ;   numGaussians - num of Gaussians to add
    ;
    ; Output - struct
    ;   signal       - Map containing noise + added Gaussians
    ;   x0           - x coords of peaks
    ;   y0           - y coords of peaks
    ;
    ; wrapper around addGauss to add random Gaussians


edgeRange = 3 * sigm        ; do not make Gaussians too close to edge
if keyword_set(minsep) then minDist = minsep else minDist = 0 * sigm
    ; do not make Gaussians too close together
range = double(sqrt(n_elements(noise)))


; set up some stuff to compute Gaussians
x = RANDOMU(seed) * (range - 2 * edgeRange - 1) + edgeRange
y = RANDOMU(seed) * (range - 2 * edgeRange - 1) + edgeRange 
x0 = [x] ; store generated coords
y0 = [y]
signal = addGauss(gaussHeight, sigm, x, y, noise)

for i=2, numGaussians do begin
    repeat begin ; generate new coords at least once
        x = RANDOMU(seed) * (range - 2 * edgeRange - 1) + edgeRange
        y = RANDOMU(seed) * (range - 2 * edgeRange - 1) + edgeRange 
    endrep until min((x0 - x)^2 + (y0 - y)^2) gt minDist^2 ; until mindist exceeded
    x0 = [x0,x]
    y0 = [y0,y]
    signal = addGauss(gaussHeight, sigm, x, y, signal)
endfor

return, {signal:signal, x0:x0, y0:y0}
end
