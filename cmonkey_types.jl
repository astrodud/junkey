## CONSTANTS AND TYPES

## CONSTANTS only defined here if they are likely to be changed (organism specific)

type bicluster
    k::Int64
    rows::Vector{Int64}
    cols::Vector{Int64}
    var::Float32
    resid::Float32
    dens_string::Float32
    meanp_meme::Float32
    scores_r::Vector{Float32}
    scores_c::Vector{Float32}
    scores_n::Vector{Float32}
    scores_m::Vector{Float32}
    meme_out::Vector{ASCIIString}
    mast_out::DataFrame ##Array{Any,2}
    changed::BitArray{1} ## rows, cols
end

bicluster( k, rows::Vector, cols::Vector ) = bicluster( k, rows, cols,
                                 ##unique( int64( [rand(1:size(x,2)) for i=1:div(size(x,2),2)] ) ),
                                 typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), 
                                 Array(Float32,0), Array(Float32,0), Array(Float32,0), 
                                 Array(Float32,0), Array(ASCIIString,0), DataFrame(), trues(2) )

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

bicluster( k, rows::Vector, x::NamedMatrix ) = bicluster( k, rows, unique(randi(size(x,2),div(size(x,2),2))) )
