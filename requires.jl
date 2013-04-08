#Pkg.add( "Distributions" ) ## only need once
using Distributions
#import Distributions.sample
using DataFrames
#using Clustering  ## will use this after upgrade to 0.2.0
using GZip  

try
    reload( "./cmonkey/cmonkey_types.jl" ) ## at end so this file can be reload()-ed without an error
catch x
    println( "Not re-loading cmonkey_types.jl" )
    println( x )
end

require( "./cmonkey/ArrayView.jl" )
require( "./cmonkey/util.jl" )
require( "./cmonkey/kmeans.jl" )
require( "./cmonkey/cmonkey_funcs.jl" )
require( "./cmonkey/scores.jl" )
require( "./cmonkey/cmonkey_init.jl" )
#require( "./cmonkey/kmeans2.jl" )
require( "./cmonkey/sequence.jl" )
require( "./cmonkey/floc.jl" )
require( "./cmonkey/gibbs-sampler-mult.jl" )
require( "./cmonkey/meme.jl" ) ## run meme on most "flagellar-enriched" cluster

gc()
