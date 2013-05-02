require( "./junkey/requires.jl" ) ## This code will be loaded onto all nodes

## MAIN PROGRAM

if ! isdefined(:organism)
    #organism = "Hpy"; k_clust = 75
    organism = "Eco"; k_clust = 450 
    #organism = "Sce"; k_clust = 500
    #organism = "Mpn"; k_clust = 75
end

iter = 1

if ! isdefined(:n_iters) n_iters = 100; end

## these are defaults (for microbes):
if ! isdefined(:distance_search) distance_search = [-150,+20]; end 
if ! isdefined(:distance_scan)   distance_scan =   [-250,+50]; end
if ! isdefined(:motif_width_range) motif_width_range = [6,24]; end
initialize_constants( organism ) ## update the three constants (above) for organism specific (e.g. yeast)

if ! isinteractive() stop(); end ## This stuff below should ONLY run on the head node

(ratios, genome_seqs, anno, op_table, string_net, allSeqs_fname, all_bgFreqs, all_genes) = junkey_init(organism, k_clust);

## DONE: set "default" k_clust if not already set:
if ! isdefined(:k_clust) k_clust = round( size( ratios, 1 ) / 10 ); end

if nprocs() > 1 pre_load_cluster_nodes(); end ## send ratios, string_net, k_clust, etc. over to children
gc()

startTime = time()
srand( 10 ) ## AFAIK random numbers are only used on the head node (during init/reseed and floc update)
clusters = init_biclusters( ratios, k_clust );
if nprocs() > 1 clusters = fill_all_cluster_scores_parallel( clusters, true, true );
else clusters = fill_all_cluster_scores( clusters, true, true ); end
println( @sprintf( "%.3f", (time() - startTime)/60 ), " minutes since initialization" )

stats_df = DataFrame()
run_junkey()

println( "DONE!" )
endTime = time()
println( @sprintf( "%.3f", (endTime - startTime)/60 ), " minutes since initialization" )

kInd = 1;
if organism == "Hpy" || organism == "Eco" kInd = clusters_w_func("flagell", clusters)[1]; end
b = clusters[kInd]
#b = re_meme_bicluster( b );
print(join(b.meme_out, '\n'))

## Save junkey code for safe keeping
junkey_code = load_junkey_code("junkey")

## Right now this fails on the clusters if we use gzopen (open works fine, though)
save_jld( "$(organism)_out.jldz", (organism, k_clust, ratios, genome_seqs, anno, op_table, string_net, 
                                   allSeqs_fname, all_bgFreqs, startTime, endTime,
                                   all_genes, iter, n_iters, distance_search, distance_scan, motif_width_range,
                                   clusters, stats_df, junkey_code) )

clusters_tab = clusters_to_dataFrame(clusters);
write_table("clusters.tsv", clusters_tab)

tmp = get_cluster_row_counts(clusters);
println( sum(tmp.==0), " genes in no clusters" )
println( sum(tmp.==max(tmp)), " genes in ", max(tmp), " clusters" )

println( @sprintf( "%.3f", (endTime - startTime)/60 ), " minutes since initialization" )


#genes = rownames(ratios)[clusters[kInd].rows] ##rows]
#seqs = get_sequences(genes);
#@time gibbs_out = gibbs_site_sampler(seqs[:,2])     ## run gibbs sampler on most "flagellar-enriched" cluster
#@time gibbs_out2 = gibbs_site_sampler(seqs, gibbs_out["pssm"])
