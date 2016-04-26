function bbody_multi, sig, specdens0, sigm, binwidth, e_in, t_in, normband=normband, autoder=autoder, real_x=real_x, real_y=real_y, real_amp=real_amp,emis_fix=emis_fix, td_fix=td_fix
compile_opt idl2, HIDDEN
; usual inputs
; computes emissivity + tdust by calling tnmin

if e_in eq !NULL then emissivity = 1.5D else emissivity = double(e_in) ; init values for tnmin search
if t_in eq !NULL then tdust = 15 else tdust = double(t_in)
if ~keyword_set(autoder) then autoder=0
limits_em = [0,4] ; limits on emissivity
limits_tdust = [3, 150 ] ; limits on tdust
if ~keyword_set(normband) then normband = 0 ; normalize to leading band by default

; hackish... for passing to tnmin
if ~keyword_set(real_x) then real_x=0
if ~keyword_set(real_y) then real_y=0
if ~keyword_set(real_amp) then real_amp=0
if ~keyword_set(emis_fix) then emis_fix=fix(0)
if ~keyword_set(td_fix) then td_fix=fix(0)

num_iters = 0
params = tnmin('tnmin_wrap', [emissivity, tdust],$
    functargs={signal:sig, specdens:specdens0, sigm:sigm, binwidth:binwidth, bbody:1, real_x:real_x, real_y:real_y, real_amp:real_amp, normband:normband},$
    parinfo=[{limited:[1,1], limits:limits_em, fixed:emis_fix}, {limited:[1,1], limits:limits_tdust, fixed:td_fix}],$
    autoderivative=autoder, niter=num_iters, /quiet)
if n_elements(params) eq 2 then begin ; if params is not effed up
    ret = subtractmax_multi(sig, specdens0, sigm, binwidth, params[0], params[1], real_x=real_x, real_y=real_y, real_amp=real_amp, normband=normband, /bbody)
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
        sigm_tdust:ret.sigm_tdust,$
        dbeta:ret.dbeta,$
        dt_dust:ret.dt_dust}
endif else return, {tdust:-1}

end
