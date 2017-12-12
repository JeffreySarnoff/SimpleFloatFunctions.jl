# SimpleFloatFunctions.jl
Extra floating point functions

## Exports

square, cube, invsquare, invcube, invsqrt, invcbrt, spread, tld, sld

spread(x) complements trunc(x): nearest integer away from zero

tld(x,y) is like cld(x,y), fld(x,y) for trunc: trunc(x/h)

sld(x,y) is like cld(x,y), fld(x,y) for spread: spread(x/h)

## Accuracy

- cube and invcube are more accurate than their obvious implementations

- the other exports are at least as accurate as their obvious implementations
