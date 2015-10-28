function fitquad, signal
compile_opt idl2, HIDDEN
    ; helper to fit log(gaussian) = quad
range = SQRT(N_ELEMENTS(signal))

x = DINDGEN(range) ## REPLICATE(1D, range)
y = REPLICATE(1D, range) ## DINDGEN(range) 
xy = x^2 + y^2 ; x^2 + y^2 term

; matrix equation Ax = b
b = dblarr(4)
A = dblarr(4, 4)

b[0] = TOTAL(signal)
b[1] = -TOTAL(signal * xy)
b[2] = TOTAL(signal * x)
b[3] = TOTAL(signal * y)

A[0,0] = range^2
A[0,1] = -TOTAL(xy)
A[0,2] = TOTAL(x)
A[0,3] = TOTAL(y)
A[1,0] = A[0,1]
A[1,1] = TOTAL(xy^2)
A[1,2] = -TOTAL(x * xy)
A[1,3] = -TOTAL(y * xy)
A[2,0] = A[0,2]
A[2,1] = A[1,2]
A[2,2] = TOTAL(x^2)
A[2,3] = TOTAL(x * y)
A[3,0] = A[0,3]
A[3,1] = A[1,3]
A[3,2] = A[2,3]
A[3,3] = TOTAL(y^2)

return, INVERT(A) # b
end
