#Pkg.add( "Distributions" ); Pkg.add( "DataFrames" ); Pkg.add( "Clustering" ) ## only need once
using Distributions
#import Distributions.sample
using DataFrames
using GZip
##TODO: use NamedArrays package...
## n = NamedArray(rand(2,4),(["A","B"],["Q","R","S","T"]),("string","string")

try
    include( "NamedMatrix.jl" ) ## catch it so the rest of this file can be require()-ed without an error
    include( "Bicluster.jl" )
#    include( "ArrayView.jl" )
catch x
    println( "Not re-loading types" )
    println( x )
end

include( "util.jl" )
#include( "kmeans.jl" )
include( "funcs.jl" )
include( "scores.jl" )
include( "init.jl" )
#include( "kmeans2.jl" )
include( "sequence.jl" )
include( "floc.jl" )
#include( "gibbs-sampler-mult.jl" )
include( "meme.jl" )

gc()
