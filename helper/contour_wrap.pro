function contour_wrap, data, ps_file=ps_file, x_bin=x, y_bin=y, title=title, xtitle=xtitle, ytitle=ytitle, xstep=xstep, ystep=ystep, gauss=gauss, quad=quad, levels=levels, c_annotation=c_annotation
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
;   levels: levels for contour plot
;       NOTE: overrides both /gauss, /quad
;   /gauss: Assume contour is of gaussian, choose levels at 0-5 sigma
;       NOTE: overrides /quad
;   /quad: Assume contour is quadratic, choose levels [1.5, 2, 2.5...]^2
;   c_annotation: annotations of contour levels, just passed on to contour procedure

!P.FONT = 1
loadct, 1
size_dat = size(data)
if size_dat[0] ne 2 then begin
    print, "Error: Wrong dimensionality of Data; should be 2D"
    return, -1
endif
; if no step or coordinates are provided, auto-generate [-n/2 - 1, ... n/2] for domain
if ~keyword_set(x) and ~keyword_set(xstep) then x = indgen(size_dat[1]) - size_dat[1] / 2
if ~keyword_set(y) and ~keyword_set(xstep) then y = indgen(size_dat[2]) - size_dat[2] / 2
; if step is provided and domain not provided, equidistant points centered on zero
if keyword_set(xstep) and ~keyword_set(x) then x = xstep * (indgen(size_dat[1]) - size_dat[1] / 2)
if keyword_set(ystep) and ~keyword_set(y) then y = ystep * (indgen(size_dat[2]) - size_dat[2] / 2)
if ~keyword_set(title) then title=''
if ~keyword_set(xtitle) then xtitle=''
if ~keyword_set(ytitle) then ytitle=''

; compute levels if gauss or quad are set and levels are not set
if ~keyword_set(levels) then begin
    if keyword_set(gauss) then begin
        ; 1, 2, 3, 4 sigma
        levels = max(data) * reverse([0.606531, 0.13534, 0.0111090, 0.0003355])
    endif else if keyword_set(quad) then begin
        levels = min(data) * [1.02^2, 1.04^2, 1.06^2, 1.08^2, 1.1^2]
    endif
endif

; if plotting to file, then:
if keyword_set(ps_file) then begin
    ; get dimensions
    set_plot, 'X'
    device, decomposed=0
    tek_color
    window, 0, xsize=900, ysize=800
    pageInfo = pswindow()
    cleanplot, /silent
    wdelete

    set_plot, 'PS'
    device, _Extra = pageInfo, /color, filename=ps_file, language=2, SET_FONT='Helvetica Bold', /TT_font

    ; plot contour
    !P.MULTI = [0, 1, 1] ; grids window into top, bottom [win_number, cols, rows]
    ; if levels is set, lpot
    if levels eq !NULL then begin
        contour, data, x, y, charsize=1.5, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, font=1, ythick=4, /fill, c_colors=[230:130:-20]
    endif else begin
        contour, data, x, y, charsize=1.5, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4, font=1, levels=levels, /fill, c_colors=[230:130:-100/(n_elements(levels) - 1)]
    endelse

    ; close file
    device, /close_file

    ; restore plotting method
    set_plot, 'X'
endif else begin
    if levels eq !NULL then begin
        if ~keyword_set(c_annocation) then begin
            contour, data, x, y, charsize=1.5, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, font=1, ythick=4, /fill, c_colors=[230:130:-20]
        endif else begin
            contour, data, x, y, charsize=1.5, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4, font=1, c_annotation=c_annotation, /fill, c_colors=[230:130:-20]
        endelse
    endif else begin
        if ~keyword_set(c_annocation) then begin
            contour, data, x, y, charsize=1.5, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4, font=1, levels=levels, /fill, c_colors=[230:130:-100/(n_elements(levels) - 1)]
        endif else begin
            contour, data, x, y, charsize=1.5, title=title, ytitle=ytitle, xtitle=xtitle, thick=4, xthick=4, ythick=4, levels=levels, font=1, c_annotation=c_annotation, /fill, c_colors=[230:130:-100/(n_elements(levels) - 1)]
        endelse
    endelse
endelse

if levels ne !NULL then return, levels else return, -1

end
