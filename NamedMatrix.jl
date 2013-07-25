##require( "base.jl" ) ## required to overrequire 'size' and other funcs
#import Base.size, Base.ndims, Base.print_matrix, Base.isnan, Base.ref, Base.convert

## Let's make rnames/cnames possible! Use a Dict for the rnames -> int lookup
# rnames = Dict{ASCIIString,Int64}()
# for i=1:size(x,1) rnames[x_rnames[i]] = i; end
# cnames = Dict{ASCIIString,Int64}()
# for i=1:size(x,2) cnames[x_cnames[i]] = i; end

# ref(a::Matrix, i::ASCIIString, j::ASCIIString) = a[rnames[i],cnames[j]]
# ref(a::Matrix, I::Vector{ASCIIString}, J::Vector{ASCIIString}) = [ a[i,j] for i=I, j=J ]

## Now can use x[["a","b"],["c"]]
## Need to turn this into a class/object.
## let's try...
type NamedMatrix{T} <: AbstractArray{T}
    x::Matrix{T}
    rnames::Dict{ASCIIString,Int64}
    cnames::Dict{ASCIIString,Int64}
    
    function NamedMatrix{T}( xx::Matrix{T}, rn::Vector{ASCIIString}, cn::Vector{ASCIIString} )
        rnd = Dict{ASCIIString,Int64}()
        for i=1:size(rn,1) rnd[rn[i]] = i; end
        cnd = Dict{ASCIIString,Int64}()
        for i=1:size(cn,1) cnd[cn[i]] = i; end
        new( xx, rnd, cnd )
    end

    ## how do you fuckn convert Any -> ASCIIString? (see below for anno table)
#     function NamedMatrix( xx, rn::Array{Any,1}, cn::Array{Any,1} ) 
#         rnd = Dict{ASCIIString,Int64}()
#         for i=1:size(rn,1) rnd[rn[i]] = i; end
#         cnd = Dict{ASCIIString,Int64}()
#         for i=1:size(cn,1) cnd[cn[i]] = i; end
#         new( xx, rnd, cnd )
#     end
end

## But do still need to implement indexing x by the value output from findn()
# function ref(x::Matrix, inds::(Array{Int64,1},Array{Int64,1}))
#     o=Array(typeof(x[1]),size(inds[1],1))
#     for i=1:size(inds[1],1)
#         o[i]=x[inds[1][i],inds[2][i]]
#     end
#     o
# end

# size(x::NamedMatrix,ind::Int64) = size(x.x,ind)
# size(x::NamedMatrix) = size(x.x)
# nrow(x::NamedMatrix) = size(x,1)
# ncol(x::NamedMatrix) = size(x,2)
# ndims(x::NamedMatrix) = ndims(x.x)
# function print_matrix(IOStream,x::NamedMatrix) 
#     print_matrix(IOStream,x.x) ## Update it to print row/col names too
#     println()
#     println(IOStream,x.rnames)
#     println(IOStream,x.cnames)
# end

# ref(x::NamedMatrix, i::ASCIIString, j::ASCIIString) = x.x[x.rnames[i],x.cnames[j]]
# ref(x::NamedMatrix, I::Vector{ASCIIString}, J::Vector{ASCIIString}) =
#                    NamedMatrix{eltype(x.x)}( [ x[i,j] for i=I, j=J ], I, J )
# ref(x::NamedMatrix, i::ASCIIString, J::Vector{ASCIIString}) = x[[i],J]
# ref(x::NamedMatrix, I::Vector{ASCIIString}, j::ASCIIString) = x[I,[j]]
# ref(x::NamedMatrix, i::Range1{Int64}, j::ASCIIString) = x.x[i,x.cnames[j]]
# ref(x::NamedMatrix, i::ASCIIString, j::Range1{Int64}) = x.x[x.rnames[i],j]
# ref(x::NamedMatrix, i::Int64, j::ASCIIString) = x.x[i,x.cnames[j]]
# ref(x::NamedMatrix, i::ASCIIString, j::Int64) = x.x[x.rnames[i],j]
## How to allow combos of string vectors and int indexes? Seems complicated.
## How to allow combos of string vectors and ranges? Seems complicated.
## How to allow combos of int vectors and/or ranges and/or int indexes? Seems complicated.

function rownames{T<:NamedMatrix}( x::T ) ## Could be faster but this is easy as long as its not called too often
    #d = reverse_dict( x.rnames )
    #out = [ d[i] for i=1:length(d) ]
    #convert( Array{ASCIIString,1}, out )
    ##keys(x.rnames)[order(values(x.rnames))] ## faster!
    keyz = Array(ASCIIString, length(x.rnames))      ## faster still!!
    valz = Array(Int64, length(x.rnames))
    i = 0
    for (k,v) in x.rnames 
        keyz[i+=1] = k 
        valz[i] = v
    end
    keyz[order(valz)]
end

function colnames{T<:NamedMatrix}( x::T ) ## Could be faster but this is easy as long as its not called too often
    #d = reverse_dict( x.cnames )
    #out = [ d[i] for i=1:length(d) ]
    #convert( Array{ASCIIString,1}, out )
    #keys(x.cnames)[order(values(x.cnames))] ## faster!
    keyz = Array(ASCIIString, length(x.cnames))     ## faster still!!
    valz = Array(Int64, length(x.cnames))
    i = 0
    for (k,v) in x.cnames 
        keyz[i+=1] = k 
        valz[i] = v
    end
    keyz[order(valz)]
end
