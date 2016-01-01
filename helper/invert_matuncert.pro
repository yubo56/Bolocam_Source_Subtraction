function invert_matuncert, mat, uncerts
; performs 2x2 matrix inversion with uncertainties
; article is at http://arxiv.org/pdf/hep-ex/9909031v1.pdf
; 
; Specifically, we assume no covariance in matrix elements, and so use
; Equation B-7 from the aforementioned article

inv = invert(mat) ; matrix inverse, i.e. inv # mat = mat # inv = [[1,0],[0,1]]
inv_uncerts = make_array(2,2) ; initialize it

; use the formula
inv_uncerts[0,0] = sqrt(inv[0,0]^4 * uncerts[0,0]^2 + $
                        inv[0,0]^2 * inv[1,0]^2 * uncerts[0,1]^2 + $
                        inv[0,0]^2 * inv[0,1]^2 * uncerts[1,0]^2 + $
                        inv[1,0]^2 * inv[0,1]^2 * uncerts[1,1]^2)
inv_uncerts[1,0] = sqrt(inv[1,0]^4 * uncerts[1,0]^2 + $
                        inv[1,0]^2 * inv[0,0]^2 * uncerts[1,1]^2 + $
                        inv[1,0]^2 * inv[1,1]^2 * uncerts[0,0]^2 + $
                        inv[0,0]^2 * inv[1,1]^2 * uncerts[0,1]^2)
inv_uncerts[0,1] = sqrt(inv[0,1]^4 * uncerts[0,1]^2 + $
                        inv[0,1]^2 * inv[1,1]^2 * uncerts[0,0]^2 + $
                        inv[0,1]^2 * inv[0,0]^2 * uncerts[1,1]^2 + $
                        inv[1,1]^2 * inv[0,0]^2 * uncerts[1,0]^2)
inv_uncerts[1,1] = sqrt(inv[1,1]^4 * uncerts[1,1]^2 + $
                        inv[1,1]^2 * inv[0,1]^2 * uncerts[1,0]^2 + $
                        inv[1,1]^2 * inv[1,0]^2 * uncerts[0,1]^2 + $
                        inv[0,1]^2 * inv[1,0]^2 * uncerts[0,0]^2)

return, {inv:inv, uncerts:inv_uncerts}
end
