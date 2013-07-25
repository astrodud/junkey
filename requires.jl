#Pkg.add( "Distributions" ) ## only need once
using Distributions
#import Distributions.sample
using DataFrames
using Clustering  ## will use this after upgrade to 0.2.0
using GZip

try
    require( "./junkey/NamedMatrix.jl" ) ## catch it so the rest of this file can be require()-ed without an error
    require( "./junkey/Bicluster.jl" )
    require( "./junkey/ArrayView.jl" )
catch x
    println( "Not re-loading types" )
    println( x )
end

require( "./junkey/util.jl" )
#require( "./junkey/kmeans.jl" )
require( "./junkey/funcs.jl" )
require( "./junkey/scores.jl" )
require( "./junkey/init.jl" )
#require( "./junkey/kmeans2.jl" )
require( "./junkey/sequence.jl" )
require( "./junkey/floc.jl" )
#require( "./junkey/gibbs-sampler-mult.jl" )
require( "./junkey/meme.jl" )

gc()
