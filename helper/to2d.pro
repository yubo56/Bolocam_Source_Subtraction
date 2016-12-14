function to2d, range, index
; for a (range)x(range) map, returns the 2D coordinates corresponding to the 1D index index
compile_opt idl2, HIDDEN
    ; simply returns the 2D index given the side length and the 1d index
    range2 = long(range) ; just in case, change to fix since relies on int div
    index2 = long(index)
    return, [index2 MOD range2, index2 / range2]
end
