pro hist_wrap, data, ps_file=ps_file,  title=title, xtitle=xtitle, ytitle=ytitle, xstep=xstep, ystep=ystep, theory=theory
; wrapper around contour to write to filename
;   data: data to plot
;   ps_file: filename to write to
;   x_bin: x-axis bins
;   y_bin: y-axis bins
;   title: title of plot
;   xtitle: xaxis title
;   ytitle: yaxis title
;   xstep: x-axis steps, auto-generate bins centered at 0
;   ystep: y-axis steps, auto-generate bins centered at 0
;   theory: a theoretical value to display on the graph (passed as double)



; check data dimensions and fill in un-supplied parameters
data = reform(data) ; so that when receiving 1D slices of multi-D data, automatically reform
size_dat = size(data)
if size_dat[0] ne 1 then begin
    print, "Error: Wrong dimensionality of Data, should be 1D"
    return
endif



; get histogram + gaussian fit
hist = histogram(data, binsize=max([stddev(data)/10, 1e-8]), locations = bins)
if n_elements(hist) lt 2 then begin
    print, "Mean: " + string(mean(data)) + " and zero deviation, not plotting..."
    return
endif
gfit = gaussfit(bins, hist, coeff, nterms=3)

; if plotting to file, then:
if keyword_set(ps_file) then begin
    ; get dimensions
    set_plot, 'X'
    device, decomposed=0
    tek_color
    window, 0, xsize=800, ysize=800
    pageInfo = pswindow()
    cleanplot, /silent
    wdelete

    set_plot, 'PS'
    device, _Extra = pageInfo, /color, filename=ps_file, language=2

    ; plot contour
    !P.MULTI = [0, 1, 1] ; grids window into top, bottom [win_number, cols, rows]
    ; if levels is set, lpot
    plot, bins, hist, charsize=1.5, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4, psym=10
    oplot, bins, gfit, color=2, thick=6
    plot_xpos = !X.CRANGE[0] + 0.05 * (!X.CRANGE[1] - !X.CRANGE[0]); coordinates for top left corner
    plot_ypos = !Y.CRANGE[0] + 0.9 * (!Y.CRANGE[1] - !Y.CRANGE[0])
    str = 'mean = ' + string(format = '(G8.2,"!C")', coeff[1])
    str = str + 'rms = ' + string(format = '(E11.4, "!C")', coeff[2])
    if keyword_set(theory) then begin
        str = str + 'theory = ' + string(format = '(E11.4)', theory)
    endif
    xyouts, /data, plot_xpos, plot_ypos, str, charsize=1.5

    ; close file
    device, /close_file

    ; restore plotting method
    set_plot, 'X'
endif else begin
    plot, bins, hist, charsize=1.5, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4, psym=10
    oplot, bins, gfit, color=2, thick=6
    plot_xpos = !X.CRANGE[0] + 0.05 * (!X.CRANGE[1] - !X.CRANGE[0]); coordinates for top left corner
    plot_ypos = !Y.CRANGE[0] + 0.9 * (!Y.CRANGE[1] - !Y.CRANGE[0])
    str = 'mean = ' + string(format = '(E11.4,"!C")', coeff[1])
    str = str + 'rms = ' + string(format = '(E11.4, "!C")', coeff[2])
    if keyword_set(theory) then begin
        str = str + 'theory = ' + string(format = '(E11.4)', theory)
    endif
    xyouts, /data, plot_xpos, plot_ypos, str, charsize=1.5
endelse

print, coeff[1], coeff[2]

end
