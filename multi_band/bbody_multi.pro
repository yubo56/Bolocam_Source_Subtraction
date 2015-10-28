function bbody_multi, sig, specdens0, sigm, binwidth, tdust = tdust, emissivity=emissivity
compile_opt idl2, HIDDEN
; usual inputs
; computes emissivity + tdust by calling tnmin

if ~ keyword_set(emissivity) then emissivity = 1.5D else emissivity = double(emissivity) ; init values for tnmin search
if ~ keyword_set(tdust) then tdust = 40 else tdust = double(tdust)
limits_em = [0,3] ; limits on emissivity; 1.5 +- 1.5
limits_tdust = [5, 100 ] ; limits on tdust

num_iters = 0
params = tnmin('tnmin_wrap', [emissivity, tdust], functargs={signal:sig, specdens:specdens0, sigm:sigm, binwidth:binwidth, bbody:1}, autoderivative=0, parinfo=[{limited:[1,1], limits:limits_em}, {limited:[1,1], limits:limits_tdust}], niter=num_iters, /quiet)
if n_elements(params) eq 2 then begin ; if params is not effed up
    ret = subtractmax_multi(sig, specdens0, sigm, binwidth, emissivity = params[0], tdust=params[1], /bbody)
    return, {xparam:ret.xparam,$
        yparam:ret.yparam,$
        Aest:ret.Aest,$
        sig:ret.sig,$
        sigm_x0:ret.sigm_x0,$
        sigm_a:ret.sigm_a,$
        chi2:ret.chi2,$
        emissivity:params[0],$
        tdust:params[1],$
        num_iters:num_iters,$
        sigm_beta:ret.sigm_beta,$
        dbeta:ret.dbeta,$
        dt_dust:ret.dt_dust}
endif else return, {tdust:-1}

end
