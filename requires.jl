#Pkg.add( "Distributions" ) ## only need once
using Distributions
#import Distributions.sample
using DataFrames
#using Clustering  ## will use this after upgrade to 0.2.0
using GZip  

try
    reload( "./junkey/junkey_types.jl" ) ## at end so this file can be reload()-ed without an error
catch x
    println( "Not re-loading junkey_types.jl" )
    println( x )
end

require( "./junkey/ArrayView.jl" )
require( "./junkey/util.jl" )
require( "./junkey/kmeans.jl" )
require( "./junkey/junkey_funcs.jl" )
require( "./junkey/scores.jl" )
require( "./junkey/junkey_init.jl" )
#require( "./junkey/kmeans2.jl" )
require( "./junkey/sequence.jl" )
require( "./junkey/floc.jl" )
require( "./junkey/gibbs-sampler-mult.jl" )
require( "./junkey/meme.jl" ) ## run meme on most "flagellar-enriched" cluster

gc()
