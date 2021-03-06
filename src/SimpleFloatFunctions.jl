module SimpleFloatFunctions

export midpoint, # an artful midpoint
       modulo, nint,
       square, cube, invsquare, invcube, invsqrt, invcbrt, spread, tld, sld,

const SysFloat = Union{Float64, Float32, Float16}

include("hasnan.jl") # does a vector or an array contain any NaNs


function midpoint(x::T,y::T) where {T}
    x, y = ifelse(abs(y) > abs(x), (y, x), (x, y))
    return x - (x/2 - y/2)
end


#=
   This well-behaved bounded modulo implementation is from
   The pitfalls of verifying floating-point computations
   by David Monniaux, 2008 
   http://arxiv.org/abs/cs/0701192v5
=#
function modulo(a::T, lowerbound::T, upperbound::T) where {T}
    delta = upperbound - lowerbound
    a - (floor((a - lowerbound)/delta) * delta)
end


nint(a::T) where {T<:AbstractFloat} = (trunc(a + copysign(one(T)/2,a)))
nint(::Type{I}, a::T) where {T<AbstractFloat, I<:Integer} = I(nint(a))

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
    invsquare(x)

1 / x^2
"""
@inline function invsquare(x::T) where T<:Number
    return square(inv(x))
end

"""
    invcube(x)

1 / x^3
"""
@inline function invcube(x::T) where T<:Number
    return cube(inv(x))
end

"""
    invsqrt(x)

1 / x^(1/2)
"""
@inline function invsqrt(x::T) where T<:Number
    return sqrt(inv(x))
end

"""
    invcbrt(x)

1 / x^(1/3)
"""
@inline function invcbrt(x::T) where T<:Number
    return cbrt(inv(x))
end

"""
    spread(x)
    
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

#=
    errorfree transformations (internal use)
=#

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
    return s, e
end

"""
    two_diff_hilo(a, b)
    t = lo + t
247
    
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
    two_prod(a, b)
Computes `s = fl(a*b)` and `e = err(a*b)`.
"""
@inline function two_prod(a::T, b::T) where T<:SysFloat
    p = a * b
    e = fma(a, b, -p)
    p, e
end
cube
"""
    one_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`.
"""
@inline function one_square(a::T) where T<:SysFloat
    p = a * a
    e = fma(a, a, -p)
    return p, e
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

"""
    one_cube_hi(a)

Computes `s = fl(a*a*a)`.
"""
@inline function one_cube_hi(a::T) where T<:SysFloat
    hi, lo = one_square(a)
    hi, md = two_prod(hi, a)
    md += lo * a
    hi += md
    return hi
end


# like doubledouble

# Algorithm 6 from Tight and rigourous error bounds for basic building blocks of double-word arithmetic
function add_dd_dd(xhi::T, xlo::T, yhi::T, ylo::T) where T<:SysFloat
    hi, lo = two_sum(xhi, yhi)
    thi, tlo = two_sum(xlo, ylo)
    c = lo + thi
    hi, lo = two_sum_hilo(hi, c)
    c = tlo + lo
    hi, lo = two_sum_hilo(hi, c)
    return hi, lo
end

function prod_dd_fl(ahi, alo, b)
    hi, lo = two_prod(ahi, b)
    lo += alo * b
    hi, lo = two_sum_hilo(hi, lo)
    return hi, lo
end

function prod_dd_fl_hi(ahi, alo, b)
    hi, lo = two_prod(ahi, b)
    lo += alo * b
    hi += lo
    return hi
end

# Algorithm 12 from Tight and rigourous error bounds for basic building blocks of double-word arithmetic
function prod_dd_dd(xhi::T, xlo::T, yhi::T, ylo::T) where T<:SysFloat
    hi, lo = two_prod(xhi, yhi)
    t = xlo * ylo
    t = fma(xhi, ylo, t)
    t = lo + t
    t = fma(xlo, yhi, t)
    t = lo + t
    hi, lo = two_sum_hilo(hi, t)
    return hi, lo
end

function prod_dd_dd_hi(xhi::T, xlo::T, yhi::T, ylo::T) where T<:SysFloat
    hi, lo = two_prod(xhi, yhi)
    t = xlo * ylo
    t = fma(xhi, ylo, t)
    t = fma(xlo, yhi, t)
    hi += lo + t
    return hi
end


end # module SimpleFloatFunctions
