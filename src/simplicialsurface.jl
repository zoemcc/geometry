using GeometryBasics
import GeometryBasics.faces
using FileIO
using StaticArrays
using SparseArrays

struct HalfEdge
    twin::UInt64
    next::UInt64
end

twin(he::HalfEdge) = he.twin
next(he::HalfEdge) = he.next


struct SimplicialSurface{P <: AbstractPoint}
    vertices::Array{P, 1}
    halfedges::Array{HalfEdge, 1} # not implemented yet
    edges::Array{SVector{2, UInt64}, 1}
    faces::Array{SVector{3, UInt64}, 1}
    edgetoindex::Dict{SVector{2, UInt64}, UInt64}
    facetoindex::Dict{SVector{3, UInt64}, UInt64}
    edgevertexadjacency::SparseMatrixCSC{UInt8, UInt64}
    faceedgeadjacency::SparseMatrixCSC{UInt8, UInt64}
end

vertices(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.vertices
halfedges(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.halfedges
edges(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.edges
faces(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.faces
edgetoindex(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.edgetoindex
facetoindex(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.facetoindex
edgevertexadjacency(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.edgevertexadjacency
faceedgeadjacency(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.faceedgeadjacency

sv(x::Vararg{Integer}) = SVector{length(x), UInt64}(x...)
sv(x::NTuple{N, Integer}) where N = sv(x...)

sortsv(x::Vararg{Integer}) = sort(sv(x...))
sortsv(x::NTuple{N, Integer}) where N = sortsv(x...)


function constructsimplicialsurface(vertices::AbstractArray{P, 1}, tris::AbstractArray{NgonFace{3, T}}) where {P <: AbstractPoint, T} 
    normedtris = map(NgonFace{3, UInt64}, tris)
    numtris = length(normedtris)

    @show P
    faces = Array{SVector{3, UInt64}, 1}(undef, numtris)
    facetoindex = Dict{SVector{3, UInt64}, UInt64}()
    edges = Array{SVector{2, UInt64}, 1}(undef, 0)
    edgetoindex = Dict{SVector{2, UInt64}, UInt64}()
    numedges::UInt64 = zero(UInt64)
    edgevertexadjacencyCOO = Array{Tuple{UInt64, UInt64, UInt8}, 1}(undef, 0)
    faceedgeadjacencyCOO = Array{Tuple{UInt64, UInt64, UInt8}, 1}(undef, 0)

    indiceshelper = ((1, 2), (2, 3), (1, 3))

    # Iterate over all faces, constructing indices and adjacency data as we go
    for i in 1:numtris
        # Construct face data
        triindcart = sortsv(normedtris[i]...)
        faces[i] = triindcart
        facetoindex[triindcart] = i

        # Iterate over face vertices and edges
        for j in 1:3
            # Construct edge data
            ind1 = triindcart[indiceshelper[j][1]]
            ind2 = triindcart[indiceshelper[j][2]]
            edgeindcart = sv(ind1, ind2)
            @show edgeindcart

            # New edge
            if !haskey(edgetoindex, edgeindcart)
                numedges += one(UInt64)
                @show numedges
                push!(edges, edgeindcart)
                edgetoindex[edgeindcart] = numedges

                # Construct edge vertex adjacency matrix data
                push!(edgevertexadjacencyCOO, (numedges, ind1, one(UInt8)))
                push!(edgevertexadjacencyCOO, (numedges, ind2, one(UInt8)))
            end

            # Construct Adjacency matrix data
            edgeindlin = edgetoindex[edgeindcart]
            push!(faceedgeadjacencyCOO, (i, edgeindlin, one(UInt8)))
        end
    end

    @show edgevertexadjacencyCOO
    @show faceedgeadjacencyCOO

    edgevertexadjacency = sparse(unzip(edgevertexadjacencyCOO)...)
    @show edgevertexadjacency
    faceedgeadjacency = sparse(unzip(faceedgeadjacencyCOO)...)
    @show faceedgeadjacency

    #halfedges = Array{HalfEdge, 1}(undef, numedges * 2)
    halfedges = Array{HalfEdge, 1}(undef, 0)
    simplicsurf = SimplicialSurface{P}(vertices, halfedges, edges, faces, edgetoindex, facetoindex, edgevertexadjacency, faceedgeadjacency)

    #return simplicsurf
    return simplicsurf
end

function constructsimplicialsurface(mesh::AbstractMesh)
    return constructsimplicialsurface(coordinates(mesh), faces(mesh))
end

function initfortest()
    mesh = load("meshes/original/tri.ply")
    return mesh
end

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
