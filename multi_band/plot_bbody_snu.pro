function plot_bbody_snu, emis, tdust, in_amps, normbands, fn=fn, large=large, title=title, name=name
; plots SEDs for list of emis, tdusts
; plots data points using parameters from [0]-th entry in arrays
; 0th is reference
range = 256D
binwidth = 480 / range
freqs = [400, 352.94, 272.73, 230.77, 150, 100] ; frequencies in GHz
num_bands = n_elements(freqs)
specdens0 = spdensgen_multi(range, $
    [0.181, 0.137, 0.112, 0.0947, 0.049, 0.009] / 2, $
    replicate(0,num_bands), replicate(8.0/3, num_bands), binwidth) ; PSD in units of 0.01 mJy^2
noise = gennoise_multi(specdens0, binwidth, freqs)
sigm = 15D / (2 * sqrt(2 * alog(2))) ; sigm in arcsec, in terms of FWHM
sigm /= binwidth ; sigm in bins
if ~keyword_set(title) then title = 'Title'

; plotting, source properties
if ~keyword_set(large) then begin
    plotfreqs = findgen(8) * 50 + 100; [100, 1500]
endif else begin
    plotfreqs = findgen(28) * 50 + 100; [100, 1500]
endelse
plotfreqs = [plotfreqs, freqs]
plotfreqs = plotfreqs[sort(plotfreqs)]
bin400 = where(abs(plotfreqs - 400) eq 0) ; 400 bin, for normalization
colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'black', 'yellow']
amp = in_amps[0]

for i=0, n_elements(emis) - 1 do begin

    ; full 
    normbin = where(abs(plotfreqs - freqs[normbands[i]]) eq 0)
    normbin = normbin[0] ; instead of array of indicies, just want one index
    plotamps = (plotfreqs / 1500)^(3.0 + emis[i]) / (exp(0.04799 * plotfreqs / tdust[i]) - 1)
    plotamps = plotamps / plotamps[normbin]
    plotamps = plotamps * in_amps[i] * sqrt(!PI) * sigm

    if i eq 0 then begin
        if ~keyword_set(name) then begin
            p = plot(plotfreqs, plotamps, COLOR=colors[i MOD n_elements(colors)],$
               YTITLE='Flux (mJy)', XTITLE='Frequency (GHz)',TITLE=title,$
               NAME=strcompress(string([emis[i], tdust[i]])), /ylog, yrange=[min(plotamps),max(plotamps)])
        endif else begin
            p = plot(plotfreqs, plotamps, COLOR=colors[i MOD n_elements(colors)],$
               YTITLE='Flux (mJy)', XTITLE='Frequency (GHz)',TITLE=title,$
               NAME=name[i])
        endelse
        p.XRANGE = [0, max(plotfreqs) + 100]
    endif else begin
        if ~keyword_set(name) then begin
            p1 = plot(plotfreqs, plotamps, /overplot,$
                NAME=strcompress(string([emis[i], tdust[i]])), COLOR=colors[i MOD n_elements(colors)])
        endif else begin
            p1 = plot(plotfreqs, plotamps, /overplot,$
                NAME=name[i], COLOR=colors[i MOD n_elements(colors)])
        endelse
    endelse
endfor

; overplot points in signal
sig_bbody = addgauss_multi(amp, sigm, range/2 + 0.5, range/2 + 0.5, noise, emis[0], tdust[0], /bbody)
amps = amps_multi(amp, freqs, emis[0], tdust[0], /bbody)
ret = getSEDpts(sig_bbody, specdens0, sigm, binwidth, range/2 + 0.5 * [1,1])
if ~keyword_set(name) then begin
    q1 = errorplot(freqs, ret.amps * sqrt(!PI) * sigm, ret.sigms * sqrt(!PI) * sigm, SYMBOL='D', /overplot)
endif else begin
    q1 = errorplot(freqs, ret.amps * sqrt(!PI) * sigm, ret.sigms * sqrt(!PI) * sigm, SYMBOL='D', name=name[i], /overplot)
endelse
stop
q1.SYM_COLOR = 'black'
q1.ERRORBAR_COLOR = 'black'
q1.SYM_FILLED = 1
q1.SYM_SIZE = 2.0
q1.LINESTYLE='none'
; p.YRANGE = [0, max(amps + ret.sigms) * 1.1]
if ~keyword_set(fn) then q1.save, 'SEDcomp.png' else q1.save, fn
return, q1
end
; SEDcomp.png
; plot_bbody_snu, [1.7, 1.807, 0.765, 1.318], [13, 10.563, 150, 25.619], [1.71, 2.524, 0.0722, 0.651]
