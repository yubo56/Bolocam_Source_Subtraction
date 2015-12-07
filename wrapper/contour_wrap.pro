pro contour_wrap, data, FN=FN, x_bin=x, y_bin=y, title=title, xtitle=xtitle, ytitle=ytitle, xstep=xstep, ystep=ystep, gauss=gauss, quad=quad
; wrapper around contour to write to filename
;   data: data to plot
;   FN: filename to write to
;   x_bin: x-axis bins
;   y_bin: y-axis bins
;   title: title of plot
;   xtitle: xaxis title
;   ytitle: yaxis title
;   xstep: x-axis steps, auto-generate bins centered at 0
;   ystep: y-axis steps, auto-generate bins centered at 0
;   /gauss: Assume contour is of gaussian, choose levels at 0-5 sigma
;       NOTE: overrides /quad
;   /quad: Assume contour is quadratic, choose levels [1.5, 2, 2.5...]^2

size_dat = size(data)
if size_dat[0] ne 2 then begin
    print, "Error: Wrong dimensionality of Data"
    return
endif
if ~keyword_set(x) and ~keyword_set(xstep) then x = indgen(size_dat[1])
if ~keyword_set(y) and ~keyword_set(xstep) then y = indgen(size_dat[2])
if keyword_set(xstep) then x = xstep * (indgen(size_dat[1]) - size_dat[1] / 2)
if keyword_set(ystep) then y = ystep * (indgen(size_dat[2]) - size_dat[2] / 2)
if ~keyword_set(title) then title=''
if ~keyword_set(xtitle) then xtitle=''
if ~keyword_set(ytitle) then ytitle=''

; compute levels if gauss, quad are set
if keyword_set(gauss) then begin
    ; 1, 2, 3, 4 sigma
    levels = max(data) * reverse([0.606531, 0.13534, 0.0111090, 0.0003355])
endif else if keyword_set(quad) then begin
    levels = min(data) * [1.5^2, 2.0^2, 2.5^2, 3.0^2, 3.5^2]
endif

; if plotting to file, then:
if keyword_set(FN) then begin
    ; get dimensions
    set_plot, 'X'
    device, decomposed=0
    tek_color
    window, 0, xsize=800, ysize=800
    pageInfo = pswindow()
    cleanplot, /silent
    wdelete

    set_plot, 'PS'
    device, _Extra = pageInfo, /color, filename=FN, language=2

    ; plot contour
    !P.MULTI = [0, 1, 1] ; grids window into top, bottom [win_number, cols, rows]
    ; if levels is set, lpot
    if levels eq !NULL then begin
        contour, data, x, y, charsize=2.0, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4
    endif else begin
        contour, data, x, y, charsize=2.0, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4, levels=levels
    endelse

    ; close file
    device, /close_file

    ; restore plotting method
    set_plot, 'X'
endif else begin
    if levels eq !NULL then begin
        contour, data, x, y, charsize=2.0, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4
    endif else begin
        contour, data, x, y, charsize=2.0, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4, levels=levels
    endelse
endelse

end
