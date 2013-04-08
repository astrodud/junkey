import Base.size, Base.ndims, Base.similar, Base.copy, Base.ref

typealias Indices{T<:Int} Union(T, AbstractVector{T}, Range1{T})

type ArrayView{T,N,A<:AbstractArray,I<:(AbstractVector...)} <: AbstractArray{T,N}
    parent::A
    indexes::I
    dims::Dims
    
    ArrayView(p::A, i::I) = new(p, i, map(length, i))
end

view{T,N}(A::AbstractArray{T,N}, i::NTuple{N,AbstractVector}) =
    ArrayView{T,N,typeof(A),typeof(i)}(A, i)

#change integer indexes into Range1 objects
view(A::AbstractArray, i::Indices...) =
    view(A, ntuple(length(i), j -> isa(i[j],AbstractVector) ? i[j] :
                                                            (i[j]:i[j])))

view(A::ArrayView, i::Indices...) =
    view(A.parent, ntuple(length(i), j -> isa(i[j],AbstractVector) ? A.indexes[j][i[j]] :
                                                                    (A.indexes[j][i[j]]):(A.indexes[j][i[j]])))

size(s::ArrayView) = s.dims
ndims{T,N}(s::ArrayView{T,N}) = N

copy(s::ArrayView) = copy_to(similar(s.parent, size(s)), s)
similar(s::ArrayView, T::Type, dims::Dims) = similar(s.parent, T, dims)

ref(s::ArrayView) = s

ref{T}(s::ArrayView{T}, i::Indices...) =
    ref(s.parent, i) ##ntuple(length(i), j -> s.indexes[j][i[j]]))

ref{T}(s::ArrayView{T,2}, i::Indices, j::Indices) =
    s.parent[s.indexes[1][i], s.indexes[2][j]]

# ref{T}(s::ArrayView{T,1}, i::Int) = s.parent[s.indexes[1][i]]

# ref{T}(s::ArrayView{T,2}, i::Int, j::Int) =
#     s.parent[s.indexes[1][i], s.indexes[2][j]]

# ref(s::ArrayView, is::Int...) = s.parent[map(ref, s.indexes, is)...]

# ref(s::ArrayView, i::Int) = s[ind2sub(size(s), i)...]

# ref{T}(s::ArrayView{T,2}, i::AbstractVector, j::AbstractVector) =
#     s.parent[s.indexes[1][i], s.indexes[2][j]]

# ref{T}(s::ArrayView{T,2}, I::Range1{Int}, J::Range1{Int}) =
#     ref(s.parent, s.indexes[1][I], s.indexes[2][J])

# ref{T}(s::ArrayView{T,2}, i::Int, J::Range1{Int}) =
#     ref(s.parent, s.indexes[1][i], s.indexes[2][J])

# ref{T}(s::ArrayView{T,2}, I::Range1{Int}, j::Int) =
#     ref(s.parent, s.indexes[1][I], s.indexes[2][j])

# ref{T}(s::ArrayView{T,2}, I::Range1{Int}, j::Int) =
#     ref(s.parent, s.indexes[1][I], s.indexes[2][j])

# ref{T}(s::ArrayView{T,2}, I::Range{Int}, J::Range{Int}) =
#     ref(s.parent, s.indexes[1][I], s.indexes[2][J])

# ## added

# ref{T}(s::ArrayView{T,1}, I::Range1{Int}) =
#     ref(s.parent, s.indexes[1][I]) #first_index+(first(I)-1)*s.strides[1]):s.strides[1]:(s.first_index+(last(I)-1)*s.strides[1]))

# ref{T}(s::ArrayView{T,1}, I::Range{Int}) =
#     ref(s.parent, s.indexes[1][I]) #(s.first_index+(first(I)-1)*s.strides[1]):(s.strides[1]*step(I)):(s.first_index+(last(I)-1)*s.strides[1]))

# ref{T,S<:Integer}(s::ArrayView{T,1}, I::AbstractVector{S}) =
#     ref(s.parent, s.indexes[1][I])

# ref{T,S<:Integer}(s::ArrayView{T,1}, I::AbstractVector{S}) =
#     ref(s.parent, s.indexes[1][I])


if false
    x = reshape([1:100],10,10)
    xx = view(x,[1,3,5],1:3)
end
