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



end # module SimpleFloatFunctions
