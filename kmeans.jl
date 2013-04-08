##srand( 10 )

function l1dist{T}( a::Vector{T}, b::Vector{T} ) #= sum( abs( a - b ) );
    out::T = 0; ## this is about 2x faster than the one-liner above
    for i=1:size(a,1) 
        if a[i] != a[i] || b[i] != b[i] continue; end ## equiv to (but faster than) isnan(a[i])
        out += abs( a[i] - b[i] )
    end
    if out == 0 out = Inf; end
    out
end

function l2dist{T}( a::Vector{T}, b::Vector{T} ) ##= sum( abs2( a - b ) );
    out::T = 0; ## this is about 2x faster than the one-liner above
    for i=1:size(a,1) out += abs2( a[i] - b[i] ); end
    if out == 0 out = Inf; end
    out
end

function kmeans_getClass{T}( k::Int64, Xexp::Array{T,1}, centers::Array{T,2}, 
                         distfun::Function )
    dists::Array{T,1} = Array(T,k);
    for j=1:k dists[j] = distfun(centers[:,j],Xexp); end
    findmin(dists)[2];
end

function kmeans_getClass2{T}( k::Int64, Xexp::Array{T,1}, centers::Array{T,2}, 
                          distfun::Function ) 
    ## a little faster to put the distfun and findmin into a single loop?
    m::T = typemax(T) ## Plus we get min(dists) out as a bonus.
    mi::Int32 = 0
    for j=1:k 
        dist = distfun(centers[:,j],Xexp)
        if dist < m
            m = dist
            mi = j
        end
    end
    ( mi, m )
end

## Todo: multiple restarts (parallelized) and return the one with the best mean minDists
function kmeans{T}(X::Array{T,2}, k::Int, maxIterations::Int, distfun::Function)
    const nan_ = 0.0 / 0.0;
    classes::Array{Int,1} = Array(Int,size(X,2)); ##fill!(classes,Inf);
    centers::Array{T,2} = rand(size(X,1),k);
    minDists::Array{T,1} = Array(T,size(X,2)); ##zeros(Float32,size(X,2));
    oldMin::T = typemax(T);
    newMin::T = 0;

    if any(isnan(X)) X[isnan(X)] = 0; end
    while maxIterations > 0
        maxIterations -= 1;
        for i = 1:size(X,2) ## keep track of minDists to re-run with multiple starts to get best 
            classes[i], minDists[i] = kmeans_getClass2( k, X[:,i], centers, distfun );
        end 
        newMin = mean(minDists);
        println("iterations left: $maxIterations $newMin");
        if newMin >= oldMin break; end
        oldMin = newMin;
        for i = 1:k
            const inClass::Array{T,2} = X[:,find(classes .== i)];
            centers[:,i] = isempty(inClass) ? nan_ : mean(inClass,2);
        end
    end
    (centers, minDists, classes)
end

kmeans{T}(X::Array{T,2}, k::Int, maxIterations::Int) = kmeans(X, k, maxIterations, l1dist);

function kmeans{T}(X::Array{T,2}, k::Int, maxIterations::Int, nRestarts::Int, 
                distfun::Function) ## should use pmap() !!!
    out = @parallel ( collect_dict ) for i=1:nRestarts 
        kmeans( X, k, maxIterations, distfun )
    end
    out
end

kmeans{T}(X::Array{T,2}, k::Int, maxIterations::Int, nRestarts::Int) =
    kmeans(X, k, maxIterations, nRestarts, l1dist);

function kmeansclust (X, k, maxIterations) 
    nan_ = 0.0 / 0.0;
    n = size(X,2);
    classes = zeros(Int32,1,n);
    centers = rand(size(X,1),k);
    oldCenters = copy(centers);
    while (maxIterations > 0)
        println("iterations left: $maxIterations");
        maxIterations = maxIterations - 1;
        for i = 1:n
            Xexp = repmat(X[:,i],1,k);
            dists = sum(abs(centers - Xexp),1);
            classes[i] = find(min(dists) .== dists)[1];
        end
        for i = 1:k
            inClass = X[:,find(classes .== i)];
            if (isempty(inClass))
                centers[:,i] = nan_;
            else
                centers[:,i] = mean(inClass,2);
            end
        end
        if (all(oldCenters .== centers))
            break;
        end
        oldCenters = copy(centers);
    end
    (centers, classes)
end

