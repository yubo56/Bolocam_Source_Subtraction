function getcovars, x2, sigm_x2, xy, sigm_xy, y2, sigm_y2, xstep, ystep
; takes coefficients of fit to chi^2 and spits out covariance matrix terms
; x2, xy, y2 - coefficients of chi2 fit
; sigm_x2,... - uncertainties on fit coeffs
; xstep, ystep - step sizes of the chi2 surfaces

x2 = double(x2)
sigm_x2 = double(sigm_x2)
xy = double(xy)
sigm_xy = double(sigm_xy)
y2 = double(y2)
sigm_y2 = double(sigm_y2)

; correct for step size and factors, see A.3 in 2D.pdf
x2 /= xstep^2
sigm_x2 /= xstep^2
xy /= xstep * ystep * 2
sigm_xy /= xstep * ystep * 2
y2 /= ystep^2
sigm_y2 /= ystep^2

mat = [[x2, xy], [xy, y2]]
sigm_mat = [[sigm_x2, sigm_xy], [sigm_xy, sigm_y2]]

temp = invert_matuncert(mat, sigm_mat)
return, temp
end
