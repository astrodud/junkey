junkey
======

Currently works on [Julia](http://julialang.org/) [v0.1.2](https://github.com/JuliaLang/julia/tree/v0.1.2), 
requires [v0.1.2](https://github.com/JuliaLang/julia/tree/v0.1.2)-compatible 
[DataFrames](https://github.com/HarlanH/DataFrames.jl) 
and [Distributions](https://github.com/JuliaStats/Distributions.jl) packages.

Fails with obscure errors on Julia 0.2.0.

<a name="Running"/>
## To run

Change to the directory **containing** the junkey directory and start julia

in Julia: 

    julia> load("junkey/main.jl")

will run the algorithm on the default (H. pylori), for 100 iterations. Will take ~an hour.
