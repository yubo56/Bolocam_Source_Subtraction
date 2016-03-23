function bbody_seekmin, sig, specdens, sigm, binwidth, binit, tinit, bstep=bstep, tstep=tstep
; sig, specdens, sigm, binwidth the usual args to subtractmax
; performs step-by-step search to brute force minimum chi2 surface
; searches a parametergrid that is bstep x tstep until it finds itself at minimum
; default: bstep = 0.0025, tstep = 0.025

if ~keyword_set(bstep) then bstep = 0.0003
if ~keyword_set(tstep) then tstep = 0.01

grid = [[1, 1], [1, 0], [1, -1], [0, 1], [0, 0], [0, -1], [-1, 1], [-1, 0], [-1, -1]]
steps = [50, 20, 5, 1]
ngrid = 9
chi2s = findgen(ngrid)

fromindex = 0
tmin = tinit
bmin = binit
for s=0, n_elements(steps) - 1 do begin
    while fromindex ne 4 do begin ; while 4 isn't the minimum
        ; iterate over parameter grid
        for i=0, ngrid - 1 do begin
            emis = bmin + bstep * grid[0, i] * steps[s]
            tdust = tmin + tstep * grid[1, i] * steps[s]

            ; subtract, get chi2s, populate grid, make step
            ret = subtractmax_multi(sig, specdens, sigm, binwidth, emis, tdust, /bbody)
            chi2s[i] = ret.chi2
        endfor
        min_chi2 = min(chi2s, fromindex)
        print, min_chi2, format='(d20.8)'
        print, bmin, tmin
        if min_chi2 eq chi2s[4] then break
        bmin += bstep * grid[0, fromindex] * steps[s]
        tmin = tmin + tstep * grid[1, fromindex] * steps[s]
    endwhile
endfor
return, {emis:bmin, tdust:tmin}
end
