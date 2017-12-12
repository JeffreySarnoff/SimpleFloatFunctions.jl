module SimpleFloatFunctions

export square, cube, invsquare, invcube, spread, tld, sld

"""
    square(x)

squares x
"""
function square(x::T) where T<:SysFloat
    return x * x
end

function square(x::T) where T<:Number
    return x * x
end

"""
    cube(x)

cubes x
"""
function cube(x::T) where T<:SysFloat
    hi, lo = one_square(x)
    hi, md = two_prod(hi, x)
    md += lo * x
    hi += md
    return hi
end

function cube(x::T) where T<:Number
    result = x*x*x
    return result
end

"""
    fourth(x)
    
x^4
"""
function fourth(x::T) where T<:SysFloat
    hi1, lo1 = one_square(x)
    hi, lo = two_prod(hi1, hi1, lo1*lo1)
    

"""
    invsquare(x)

1 / x^2
"""
function invsquare(x::T) where T<:Number
    return inv(square(x))
end

"""
    invcube(x)

1 / x^3
"""
function invcube(x::T) where T<:Number
    return inv(cube(x))
end

"""
    spread
    
spread complements trunc()    
the nearest integer to x, away from zero
"""
function spread(x::T) where T<:AbstractFloat
    return signbit(x) ? floor(x) : ceil(x)
end

"""
    tld(x, y)

truncates @ref(trunc) the result of x/y

like cld @ref(cld), fld @ref(fld)
"""
function tld(x::T) where T<:AbstractFloat
    return trunc(x/y)
end

"""
    sld(x, y)
    
spreads @ref(spred) the result of x/y

like cld @ref(cld), fld @ref(fld)
"""
function sld(x::T) where T<:AbstractFloat
    return spread(x/y)
end




# errorfree transformations (internal use only)
const SysFloat = Union{Float64, Float32, Float16}

"""
    two_sum_hilo(a, b)

*unchecked* requirement `|a| ≥ |b|`

Computes `s = fl(a+b)` and `e = err(a+b)`.
"""
@inline function two_sum_hilo(a::T, b::T) where T<:SysFloat
    s = a + b
    e = b - (s - a)
    return s, e
end

"""
    two_sum(a, b)

Computes `s = fl(a+b)` and `e = err(a+b)`.
"""
@inline function two_sum(a::T, b::T) where T<:SysFloat
    s = a + b
    v = s - a
    e = (a - (s - v)) + (b - v)
"""

# errorfree transformations (internal use only)
const SysFloat = Union{Float64, Float32, Float16}

"""
    two_sum_hilo(a, b)

*unchecked* requirement `|a| ≥ |b|`

Computes `s = fl(a+b)` and `e = err(a+b)`.
"""
@inline function two_sum_hilo(a::T, b::T) where T<:SysFloat
    s = a + b
    e = b - (s - a)
    return s, e
end

"""
    two_sum(a, b)

Computes `s = fl(a+b)` and `e = err(a+b)`.
"""
@inline function two_sum(a::T, b::T) where T<:SysFloat
    s = a + b
    v = s - a
    e = (a - (s - v)) + (b - v)
"""
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end
    return s, e
end

"""
    two_diff_hilo(a, b)
    
*unchecked* requirement `|a| ≥ |b|`

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff_hilo(a::T, b::T) where T<:SysFloat
    s = a - b
    e = (a - s) - b
    s, e
end

"""
    two_diff(a, b)

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff(a::T, b::T) where T<:SysFloat
    s = a - b
    v = s - a
    e = (a - (s - v)) - (b + v)
    s, e
end

"""
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end

"""
    two_prod(a, b)

Computes `s = fl(a*b)` and `e = err(a*b)`.
"""
@inline function two_prod(a::T, b::T) where T<:SysFloat
    p = a * b
    e = fma(a, b, -p)
    p, e
end

"""
    one_square(a)
    
Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end

"""
    one_cube(a)
    
Computes `s = fl(a*a*a)` and `e = err(a*a*a)`.
"""
@inline function one_cube(a::T) where T<:SysFloat
    hi, lo = one_square(a)
    hihi, hilo = two_prod(hi, a)
    lohi, lolo = two_prod(lo, a)
    hilo, lohi = two_sum_hilo(hilo, lohi)
    hi, lo = two_sum_hilo(hihi, hilo)
    lo += lohi + lolo
    return hi, lo
end
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end
    return s, e
end

"""
    two_diff_hilo(a, b)
    
*unchecked* requirement `|a| ≥ |b|`

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff_hilo(a::T, b::T) where T<:SysFloat
    s = a - b
    e = (a - s) - b
    s, e
end

"""
    two_diff(a, b)

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff(a::T, b::T) where T<:SysFloat
    s = a - b
    v = s - a
    e = (a - (s - v)) - (b + v)
    s, e
end

"""
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end

    hihi, hilo = two_prod(hi, a)
    lohi, lolo = two_prod(lo, a)
    hilo, lohi = two_sum_hilo(hilo, lohi)
    hi, lo = two_sum_hilo(hihi, hilo)
    lo += lohi + lolo
    return hi, lo
end
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end
    return s, e
end

"""
    two_diff_hilo(a, b)
    
*unchecked* requirement `|a| ≥ |b|`

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff_hilo(a::T, b::T) where T<:SysFloat
    s = a - b
    e = (a - s) - b
    s, e
end

"""
    two_diff(a, b)

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff(a::T, b::T) where T<:SysFloat
    s = a - b
    v = s - a
    e = (a - (s - v)) - (b + v)
    s, e
end

"""
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end

"""
    two_prod(a, b)
"""
128
    two_prod(a, b)
129
    
130
Computes `s = fl(a*b)` and `e = err(a*b)`.
131
"""
132
@inline function two_prod(a::T, b::T) where T<:SysFloat
133
    p = a * b
134
    e = fma(a, b, -p)
135
    p, e
136

Computes `s = fl(a*b)` and `e = err(a*b)`.
"""
@inline function two_prod(a::T, b::T) where T<:SysFloat
    p = a * b
    e = fma(a, b, -p)
    p, e
end

""""""
128
    two_prod(a, b)
129
    
130
Computes `s = fl(a*b)` and `e = err(a*b)`.
131
"""
132
@inline function two_prod(a::T, b::T) where T<:SysFloat
133
    p = a * b
134
    e = fma(a, b, -p)
135
    p, e
136

    one_square(a)
    
Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end

"""
    one_cube(a)
    
Computes `s = fl(a*a*a)` and `e = err(a*a*a)`.
"""
@inline function one_cube(a::T) where T<:SysFloat
    hi, lo = one_square(a)
    hihi, hilo = two_prod(hi, a)
    hihi, hilo = two_prod(hi, a)
    lohi, lolo = two_prod(lo, a)
    hilo, lohi = two_sum_hilo(hilo, lohi)
    hi, lo = two_sum_hilo(hihi, hilo)
    lo += lohi + lolo
    return hi, lo
end
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end
    return s, e
end

"""
    two_diff_hilo(a, b)
    
*unchecked* requirement `|a| ≥ |b|`

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff_hilo(a::T, b::T) where T<:SysFloat
    s = a - b
    e = (a - s) - b
    s, e
end

"""
    two_diff(a, b)

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff(a::T, b::T) where T<:SysFloat
    s = a - b
    v = s - a
    e = (a - (s - v)) - (b + v)
    s, e
end

"""
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end

"""
    two_prod(a, b)
"""
128
    two_prod(a, b)
129
    
130
Computes `s = fl(a*b)` and `e = err(a*b)`.
131
"""
132
@inline function two_prod(a::T, b::T) where T<:SysFloat
133
    p = a * b
134
    e = fma(a, b, -p)
135
    p, e
136

Computes `s = fl(a*b)` and `e = err(a*b)`.
"""
@inline function two_prod(a::T, b::T) where T<:SysFloat
    p = a * b
    e = fma(a, b, -p)
    p, e
end

""""""
128
    two_prod(a, b)
129
    
130
Computes `s = fl(a*b)` and `e = err(a*b)`.
131
"""
132
@inline function two_prod(a::T, b::T) where T<:SysFloat
133
    p = a * b
134
    e = fma(a, b, -p)
135
    p, e
136

    one_square(a)
    
Computes `s 
    hihi, hilo = two_prod(hi, a)
    lohi, lolo = two_prod(lo, a)
    hilo, lohi = two_sum_hilo(hilo, lohi)
    hi, lo = two_sum_hilo(hihi, hilo)
    lo += lohi + lolo
    return hi, lo
end
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end
    return s, e
end

"""
    two_diff_hilo(a, b)
    
*unchecked* requirement `|a| ≥ |b|`

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff_hilo(a::T, b::T) where T<:SysFloat
    s = a - b
    e = (a - s) - b
    s, e
end

"""
    two_diff(a, b)

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff(a::T, b::T) where T<:SysFloat
    s = a - b
    v = s - a
    e = (a - (s - v)) - (b + v)
    s, e
end

"""
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end

"""
    two_prod(a, b)

Computes `s = fl(a*b)` and `e = err(a*b)`.
"""
@inline function two_prod(a::T, b::T) where T<:SysFloat
    p = a * b
    e = fma(a, b, -p)
    p, e
end

"""
    one_square(a)
    
Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    p, e
end

"""
    one_cube(a)
    
Computes `s = fl(a*a*a)` and `e = err(a*a*a)`.
"""
@inline function one_cube(a::T) where T<:SysFloat
    hi, lo = one_square(a)
    hihi, hilo = two_prod(hi, a)
    lohi, lolo = two_prod(lo, a)
    hilo, lohi = two_sum_hilo(hilo, lohi)
    hi, lo = two_sum_hilo(hihi, hilo)
    lo += lohi + lolo
    return hi, lo
end






end # module SimpleFloatFunctions
