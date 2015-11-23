function tnmin_wrap, p, dp, signal=sig, specdens=specdens, sigm=sigm, binwidth=binwidth, bbody=bbody, real_x=real_x, real_amp=real_amp, real_y=real_y
compile_opt idl2, HIDDEN
; wrapper
;   p = beta
;   dp = stores derivative [optional]
;   functargs = {sig:sig, specdens:specdens, sigm:sigm, binwidth:binwidth}
;   /bbody = keyword_set(bbody) means we are doing double fit for tdust, emissivity

if keyword_set(bbody) then begin
    if keyword_set(real_x) then begin
        if keyword_set(real_y) then begin
            if keyword_set(real_amp) then begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p[0], tdust=p[1], /bbody, real_x=real_x, real_y=real_y, real_amp=real_amp)
                dp = [ret.dbeta, ret.dt_dust]
                return, ret.chi2
            endif else begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p[0], tdust=p[1], /bbody, real_x=real_x, real_y=real_y)
                dp = [ret.dbeta, ret.dt_dust]
                return, ret.chi2
            endelse
        endif else begin
            if keyword_set(real_amp) then begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p[0], tdust=p[1], /bbody, real_x=real_x, real_amp=real_amp)
                dp = [ret.dbeta, ret.dt_dust]
                return, ret.chi2
            endif else begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p[0], tdust=p[1], /bbody, real_x=real_x)
                dp = [ret.dbeta, ret.dt_dust]
                return, ret.chi2
            endelse
        endelse
    endif else begin
        if keyword_set(real_y) then begin
            if keyword_set(real_amp) then begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p[0], tdust=p[1], /bbody, real_amp=real_amp, real_y=real_y)
                dp = [ret.dbeta, ret.dt_dust]
                return, ret.chi2
            endif else begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p[0], tdust=p[1], /bbody, real_y=real_y)
                dp = [ret.dbeta, ret.dt_dust]
                return, ret.chi2
            endelse
        endif else begin
            if keyword_set(real_amp) then begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p[0], tdust=p[1], /bbody, real_amp=real_amp)
                dp = [ret.dbeta, ret.dt_dust]
                return, ret.chi2
            endif else begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p[0], tdust=p[1], /bbody)
                dp = [ret.dbeta, ret.dt_dust]
                return, ret.chi2
            endelse
        endelse
    endelse
endif else begin
    if keyword_set(real_x) then begin
        if keyword_set(real_y) then begin
            if keyword_set(real_amp) then begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p, real_x=real_x, real_y=real_y, real_amp=real_amp)
                dp = ret.dbeta
                return, ret.chi2
            endif else begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p, real_x=real_x, real_y=real_y)
                dp = ret.dbeta
                return, ret.chi2
            endelse
        endif else begin
            if keyword_set(real_amp) then begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p, real_x=real_x, real_amp=real_amp)
                dp = ret.dbeta
                return, ret.chi2
            endif else begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p, real_x=real_x)
                dp = ret.dbeta
                return, ret.chi2
            endelse
        endelse
    endif else begin
        if keyword_set(real_y) then begin
            if keyword_set(real_amp) then begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p, real_amp=real_amp, real_y=real_y)
                dp = ret.dbeta
                return, ret.chi2
            endif else begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p, real_y=real_y)
                dp = ret.dbeta
                return, ret.chi2
            endelse
        endif else begin
            if keyword_set(real_amp) then begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p, real_amp=real_amp)
                dp = ret.dbeta
                return, ret.chi2
            endif else begin
                ret = subtractmax_multi(sig, specdens, sigm, binwidth, emissivity=p)
                dp = ret.dbeta
                return, ret.chi2
            endelse
        endelse
    endelse
endelse
end
