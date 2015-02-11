## converted anno into a DataFrame -- faster???
## string_net and op_table and anno was already converted to DataFrame
## TODO: Also try using Named package instead of NamedMatrix?
## TODO: Try using BioSeq package instead of ASCIIStrings for sequences (seq sequence.jl for details)
## TODO: Load and use other networks (with their respective scores)

function load_ratios(rats_file)
    ## TODO: make it work for any organism; download the sequence/anno data using HTTP package; 
    ## DONE: try on bigger data; multiple chrome's
    println(rats_file)
    x=readdlm(rats_file,'\t',ASCIIString)
    x_rnames=vec(x[2:end,1])
    x_cnames=vec(x[1,2:end])
    x=readdlm(rats_file,'\t',skipstart=1)[:,2:end]
    x[x.=="NA"] = NaN
    x = float(x)
    println(size(x))

## Examples for Eco:
## First get the taxinomy codes:
## wget -O taxdump.tar.gz 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
## tar -xvzf taxdump.tar.gz names.dmp
## wget -O genomeInfo.511145.tsv 'http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=511145;export=tab'
## wget -O genome.511145.txt 'http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=511145;export=genome'
## wget -O STRING_511145.tsv.gz 'http://baliga.systemsbiology.net/data/STRING/511145_STRING.tsv.gz'; gunzip STRING_511145.tsv.gz
## NOTE: extract first 200 columns from a file using cut:
##    gunzip -c TanayMSB_GEData.txt.gz | cut -f1-200 | more

    ##means=rowmeans(x)
    ##sds=rowsds(x)
    #inds=findn(isnan(x)) ## same as which(...,arr=T) in R
    #vals=x[inds]

    good_rows = find( sum(isnan(x), 2) .< size(x, 2)/2 )
    x = x[good_rows, :]
    x_rnames = x_rnames[good_rows]
    ##good_cols = find( sum(isnan(x), 1) .< size(x, 1)/2 )
    println(size(x))

    for i=1:size(x,1) x[i,:] = stdize_vector(x[i,:]); end ##(x[i,:]-nanmean(x[i,:]))/nansd(x[i,:]); end

    x = NamedMatrix{Float32}( x, x_rnames, x_cnames ) ## Works! Needs more work... or use Named package?
    x
end

function load_genome(organism)
    ## DONE: allow for multiple genome seqs (e.g. Halo, Sce)
    org_files = system(`ls ./$organism/`)
    genome_file = org_files[ findfirst( [ beginswith(i,"genome.") for i=org_files ] ) ]
    genome_seqs = readFastaDNA( "./$organism/$genome_file" );
    println(size(genome_seqs))
    genome_seqs
end

function load_annos(organism)
    ## Load the gene annotations
    org_files = system(`ls ./$organism/`)
    genomeInfo_file = org_files[ findfirst( [ beginswith(i,"genomeInfo.") for i=org_files ] ) ]
    ## Note current version of readtable() specifically parses '+' and '-' as integers (see bytestoint() in io.jl)
    df = DataFrames.readtable( "./$organism/$genomeInfo_file", separator='\t' )
    #io = open("./$organism/$genomeInfo_file", "r")
    #df = DataFrames.read_table( io, '\t', DataFrames.DEFAULT_QUOTATION_CHARACTER, DataFrames.DEFAULT_MISSINGNESS_INDICATORS, true,
    #                       DataFrames.determine_column_names(io, '\t', DataFrames.DEFAULT_QUOTATION_CHARACTER, true ),
    #                       DataFrames.determine_nrows( "./$organism/$genomeInfo_file", true ) ) ## pain in the ass but it works
    #close( io )
    for i=find([typeof(df[i]) for i=1:size(df,2)].==DataArray{UTF8String,1})
        df[i] = DataArray( convert( Vector{ASCIIString}, df[ i ].data ) )
    end
    anno = df
#     anno = readdlm("./$organism/$genomeInfo_file", '\t', Any) ##ASCIIString)
#     anno = NamedMatrix{Any}( anno[2:end,2:end], 
#                             Base.convert(Array{ASCIIString,1},vec(anno[2:end,8])), 
#                             Base.convert(Array{ASCIIString,1},vec(anno[1,2:end])) )
#     anno.x[:,2] = [ typeof(anno.x[i,2])==Float64 ? int32(anno.x[i,2]) : NA for i=1:size(anno,1) ]
#     anno.x[:,3] = [ typeof(anno.x[i,3])==Float64 ? int16(anno.x[i,3]) : NA for i=1:size(anno,1) ]
#     anno.x[:,4] = [ typeof(anno.x[i,4])==Float64 ? int32(anno.x[i,4]) : NA for i=1:size(anno,1) ]
#     anno.x[:,5] = [ typeof(anno.x[i,5])==Float64 ? int32(anno.x[i,5]) : NA for i=1:size(anno,1) ]
#     anno.x[:,6] = [ typeof(anno.x[i,6])==ASCIIString ? anno.x[i,6][1] : NA for i=1:size(anno,1) ]
    
#    if organism == "Mpn"
 #       println(anno[["MPN002","MPN003"],"start"]) ## test
 #       @time for i=1:100000 tmp = anno[["MPN002","MPN003"],"start"]; end
 #   end
 #   @time for i=1:100000 tmp = anno.x[[2,3],4]; end
    anno
end

function load_string_net(organism)
    ## Load the string network
    org_files = system(`ls ./$organism/`)
    string_file = org_files[ findfirst( [ beginswith(i,"STRING") for i=org_files ] ) ]
    ##df = readdlm("./$organism/$string_file", '\t')
    df = DataFrames.readtable( "./$organism/$string_file", separator='\t', header=false )
    #io = open("./$organism/$string_file", "r")
    #has_header = false
    #df = DataFrames.read_table( io, '\t', DataFrames.DEFAULT_QUOTATION_CHARACTER, DataFrames.DEFAULT_MISSINGNESS_INDICATORS, 
    #                           has_header,
    #                       DataFrames.determine_column_names(io, '\t', DataFrames.DEFAULT_QUOTATION_CHARACTER, has_header ),
    #                       DataFrames.determine_nrows( "./$organism/$string_file", has_header ) ) ## pain in the ass but it works
    #close( io )
    # df[1] = [ split( df[ i, 1 ], '.' )[ 2 ] for i=1:size( df, 1 ) ]
    # df[2] = [ split( df[ i, 2 ], '.' )[ 2 ] for i=1:size( df, 1 ) ]
    # if size(df,2) > 3 
    #     df = df[[1,2,15]]
    #     df[3] = [ float32( replace( split( df[ i, 3 ], '|' )[ 1 ], "score:", "" ) ) for i=1:size( df, 1 ) ]
    # end
    df[3] = df[3] / float32(maximum(df[3]))
    df[1] = PooledDataArray( convert(Vector{ASCIIString}, df[1].data) ) ## use PooledDataArray to use less memory; doesn't seem to slow it down
    df[2] = PooledDataArray( convert(Vector{ASCIIString}, df[2].data) )
    string_net = df

    ## Symmetrize it -- no, seems to already be done.
    # tmp = hcat( string_net[:,2], string_net[:,1], string_net[:,3] )
    # string_net = vcat( string_net, tmp )
    println(size(string_net))

    ## Faster for searching? Yes, 5x faster.
    #string_net[:,1] = [ all_genes[g] for g=string_net[:,1] ]
    #string_net[:,2] = [ all_genes[g] for g=string_net[:,2] ]
        
    # if organism == "Mpn"
    #     ## Whew, this is some work; must be a faster way: find edge between MPN391 and MPN390 only if score >= 0.9
    #     println( string_net[find(Base.convert(Array{Bool,1},string_net[:,1].=="MPN391")'&
    #                              Base.convert(Array{Bool,1},string_net[:,2].=="MPN390")'&
    #                              (float32(string_net[:,3]).>=0.9)'),:] )
    #     ## Try it this way (function defined in junkey_funcs.jl)
    #     ##@time for i=1:100 tmp = get_subnetwork( string_net, "MPN391", "MPN390", float32(0.9) ); end
    # end
    string_net
end

function load_op_table(organism)
    ## Now, require the operons table!
    ## DONE? (need to verify): Allow for no operons table (e.g. yeast)
    org_files = system(`ls ./$organism/`)
    op_file = findfirst( [ beginswith(i,"microbesonline_operons_") for i=org_files ] )
    op_table::DataFrame = DataFrame()
    if op_file > 0 ##filesize("./$organism/$op_file") > 0   ## eukaryotes - no op_file exists!
        op_file = org_files[ op_file ]
        df = DataFrames.readtable( "./$organism/$op_file", separator='\t' )
        df = df[ [ :SysName1, :SysName2, :bOp, :pOp ] ] ##3, 4, 7, 8 ] ];
        df[ :SysName1 ] = DataArray( convert( Vector{ASCIIString}, df[ :SysName1 ].data ) )
        df[ :SysName2 ] = DataArray( convert( Vector{ASCIIString}, df[ :SysName2 ].data ) )
        df[ :pOp ] = DataArray( convert( Vector{Float32}, df[ :pOp ].data ) )
        op_table = df
    else
       warn( "NO OPERONS TABLE!" )
    end
    op_table
end

## Load the junkey code as text so it can be stored in the run results
function load_junkey_code(path)
    files = system(`ls $path/`)
    files = files[ [ endswith( f, ".jl" ) for f in files ] ]
    files = files[ files .!= "nanopond-1.9C.jl" ]
    code = { f => readall( "./$f" ) for f in files }
    code
end

function junkey_init(organism, k_clust, force=true)
    global ratios_file, ratios, distance_scan, distance_search
    println(organism)

    ## TODO: shouldn't store ratios, k_clust in this file; this is just organism/genome data
    if ! force && filesize("./$organism/data.jldz") > 0 ## 0 if not exist
        warn( "Loading organism data from ./$organism/data.jldz" )
        (organism, k_clust, ratios, genome_seqs, anno, op_table, string_net, 
         allSeqs_fname, all_bgFreqs, all_genes) = load_jld("./$organism/data.jldz"); ##, all_rows

        ## Make sure we have written out the sequences file for MAST-ing
        if filesize(allSeqs_fname) <= 0 
            all_seqs_scan = get_sequences(anno["sysName"].data,anno,genome_seqs,true,op_table,distance_scan,false); 
            all_seqs_scan = all_seqs_scan[ find(all_seqs_scan[:,1].!=""), : ]
            all_seqs_scan = filter_sequences( all_seqs_scan, distance_scan )

            writeFasta( all_seqs_scan, allSeqs_fname )
        end
        return( (ratios, genome_seqs, anno, op_table, string_net, ##all_seqs, all_seqs_scan, 
                 allSeqs_fname, all_bgFreqs, all_genes) ) ##, all_rows) ##all_bgCounts, 
    end 

    if ! isdefined(:ratios) ratios = load_ratios(ratios_file); end
    genome_seqs = load_genome(organism);
    anno = load_annos(organism);
    string_net = load_string_net(organism);
    op_table = load_op_table(organism);

    ## build up a list of all genes in the expr. data + in the annotations + in the string network
    all_genes = Dict{ASCIIString,Int64}()
    for i=keys(ratios.rnames) all_genes[i] = ratios.rnames[i]; end
    ##for i=keys(anno.rnames) if get( all_genes, i, 0 ) == 0 all_genes[i] = length(all_genes)+1; end; end
    for i=anno[:sysName].data if get( all_genes, i, 0 ) == 0 all_genes[i] = length(all_genes)+1; end; end
    for i=1:size(string_net,1)
        if get( all_genes, string_net[i,1], 0 ) == 0 all_genes[string_net[i,1]] = length(all_genes)+1; end
        if get( all_genes, string_net[i,2], 0 ) == 0 all_genes[string_net[i,2]] = length(all_genes)+1; end
    end
    ##for i=keys(op_table.rnames) if get( all_genes, i, 0 ) == 0 all_genes[i] = length(all_genes)+1; end; end
    if size(op_table,1) > 0
        for i=op_table[:SysName1].data if get( all_genes, i, 0 ) == 0 all_genes[i] = length(all_genes)+1; end; end
        for i=op_table[:SysName2].data if get( all_genes, i, 0 ) == 0 all_genes[i] = length(all_genes)+1; end; end
    end
    println(length(all_genes))

    gene_regex = get_regex(collect(keys(all_genes)))
    println(gene_regex)
    
    ##all_rows::Vector{Int} = [ all_genes[ g ] for g in rownames(x) ] ##pre-define this (not used anymore?)
    
    ## Get all upstream seqs and their bgFreqs
    ## TODO: Need to add 0's for k-mers that are NOT in the sequences.
    ##   Use generate_all_kmers() for that.
    ## TODO: Need to include vague IUPAC symbols better

    all_seqs_scan = get_sequences(anno[:sysName].data,anno,genome_seqs,true,op_table,distance_scan,false); 
    all_seqs_scan = all_seqs_scan[ find(all_seqs_scan[:,1].!=""), : ]
    all_seqs_scan = filter_sequences( all_seqs_scan, distance_scan )

    allSeqs_fname = "./$(organism)/allSeqs.fst"
    writeFasta( all_seqs_scan, allSeqs_fname ) ## NOTE all_seqs_scan are not used from here on
    
    ##bgCounts = getBgCounts( [genome_seqs,revComp(genome_seqs)], [0:5], true );
    ##bgFreqs = getBgFreqs( bgCounts ); ## TODO: bgFreqs are currently not used in MEME-ing OR MAST-ing.

    ## TODO: remove 'N's from k-mers in counts
    all_bgCounts = getBgCounts( all_seqs_scan[:seq].data );
    all_bgFreqs = getBgFreqs( all_bgCounts );  ## TODO: bgFreqs are currently not used in MEME-ing OR MAST-ing.

    save_jld( "./$organism/data.jldz", (organism, k_clust, ratios, genome_seqs, anno, op_table, string_net, 
                                          allSeqs_fname, all_bgFreqs, all_genes) ) ##, all_rows
    
    (ratios, genome_seqs, anno, op_table, string_net, ##all_seqs, all_seqs_scan, 
     allSeqs_fname, all_bgFreqs, all_genes) ##, all_rows) ##all_bgCounts, 
end

function init_biclusters( x, k_clust, method="kmeans" ) 
    ## Init via kmeans -- TODO: other methods (including random)
    ## DONE: use k_means from Clustering package
    clusters = Dict{Int64,bicluster}()

    if method == "kmeans" || method == "kmeans+random"
        xx = x.x;
        is_nan = isnan(xx)
        xx[is_nan] = rand(sum(is_nan))*0.1 - 0.05; ## randomize a bit; note this works on the original (global) x!

        ##km1=kmeansclust( xx, 50, 20 ); ## original
        ##km1=kmeans( xx', k_clust, 20 ); ## my "optimized" version that pre-removes NAs
        ##assignments = km1[3]
        km1=Clustering.kmeans(convert(Matrix{Float64},xx'), k_clust)
        assignments = km1.assignments
        
        #     if nprocs() > 1 km2=kmeans2( x.x', k_clust, 20 ); ## parallelized version -- now seems to be broken
        #     else km1=kmeans( x.x', k_clust, 20 ); ## my "optimized" version that pre-removes NAs
        #     end
        
        ## seed each bicluster with rows=output from kmeans and cols=random (1/2 of all cols)
        for k=1:k_clust
            rows = findin(assignments, k)
            if method == "kmeans" 
                clusters[k] = bicluster( k, rows, x )
            elseif method == "kmeans+random" ##         ## add some random rows to each bicluster
                clusters[k] = bicluster( k, unique([ rows, randperm(size(xx,1))[1:length(rows)] ]), x )
            end
        end
        xx[is_nan] = NA ## reset the values to NA
    elseif method == "random"
        for k=1:k_clust
            clusters[k] = bicluster( k, unique( randperm(size(x.x,1))[1:20] ), x )
        end
    end
    clusters
end

#get_regex( strings::Vector{ASCIIString} ) = get_regex( strings, 2 )

function get_regex( strings::Vector{ASCIIString}, min_ignore::Int=2 )
    nchar=int32([length(strings[i]) for i=1:length(strings)])
    out = ""
    for i=1:maximum(nchar)
        d = Dict{Char,Int64}()
        for j=1:length(strings)
            if i > nchar[j] continue; end
            c::Char = strings[j][i]
	    d[c] = get(d,c,0) + 1
        end
        if isempty(d) break; end
        tmp = ""
        for c::Char=collect(keys(d)) if get(d,c,0) > min_ignore tmp = "$tmp$c"; end; end ## ignore values with <=2 occurences
        if length(tmp) > 1 tmp = "[$tmp]"; end
        if sum(collect(values(d))) < length(strings) * 0.95 tmp = "$tmp?"; end
        if tmp != "" out = "$out$tmp"; end
    end
    out
end
