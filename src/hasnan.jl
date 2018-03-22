import Base: IEEEFloat

uint(x::Float64) = reinterpret(UInt64, x)
uint(x::Float32) = reinterpret(UInt32, x)
uint(x::Float16) = reinterpret(UInt16, x)

const nanbit64 = xor(uint(NaN), uint(Inf))
const nanbit32 = xor(uint(NaN32), uint(Inf32))
const nanbit16 = xor(uint(NaN16), uint(Inf16))

nanbit(::Type{Float64}) = nanbit64
nanbit(::Type{Float32}) = nanbit32
nanbit(::Type{Float16}) = nanbit16

@inline function isolate_nan(a::T, b::T) where {T<:IEEEFloat}
    return uint(a) | uint(b)
end

@inline function isolate_nan(a::U, b::T) where {U<:Unsigned, T<:IEEEFloat}
    return a | uint(b)
end

@inline function isolate_nan(vec::AbstractVector{T}) where {T<:IEEEFloat}
    reduce( isolate_nan, vec )
end

function hasnan(vec::AbstractVector{T}) where {T<:IEEEFloat}
     nanbit(T) == nanbit(T) & isolate_nan(vec)
end
