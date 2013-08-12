## TODO: read in these paramters from cmd-line arguments

if ! isdefined(:organism)
    #organism = "Hpy"
    organism = "Eco"
    #organism = "Sce"
    #organism = "Mpn"
end

if ! isdefined(:ratios_file) ratios_file = "./junkey/$organism/ratios.tsv"; end

if ! isdefined(:k_clust) 
    if organism == "Hpy" k_clust = 75
    elseif organism == "Eco" k_clust = 471
    elseif organism == "Sce" k_clust = 500
    elseif organism == "Mpn" k_clust = 75
        ## DONE: set "default" k_clust if not already set:
    elseif isdefined(:ratios) k_clust = round( size( ratios, 1 ) / 10 )
    end
end

if ! isdefined(:n_iters) n_iters = 100; end

## these are defaults (for microbes):
if ! isdefined(:distance_search) distance_search = [-150,+20]; end 
if ! isdefined(:distance_scan)   distance_scan =   [-250,+50]; end
if ! isdefined(:motif_width_range) motif_width_range = [6,24]; end

if organism == "Sce"
    if ! isdefined(:distance_search) distance_search = [-250,+30]; end 
    if ! isdefined(:distance_scan)   distance_scan =   [-500,+50]; end
    if ! isdefined(:motif_width_range) motif_width_range = [6,18]; end
end
println( "$distance_search  $distance_scan  $motif_width_range" )

if ! isdefined(:max_network_weight) max_network_weight = 0.03; end 
if ! isdefined(:max_motif_weight)   max_motif_weight =  0.5; end 

if ! isdefined(:avg_genes_per_cluster) avg_genes_per_cluster = 31; end 
if ! isdefined(:avg_clusters_per_gene) avg_clusters_per_gene = 0.7000000000000002; end 

## allow more updates if there are more clusters??? Tuned to k_clust/2 for Hpy (where k_clust is 75) -- 
##     may need additional tuning; e.g. for eco (k_clust=450), k_clust/2 is too high
if ! isdefined(:max_improvements_per_iter) max_improvements_per_iter = convert(Int64, round(min(k_clust/2, 100))); end
