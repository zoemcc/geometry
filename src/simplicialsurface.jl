using GeometryBasics
import GeometryBasics.faces
using FileIO
using StaticArrays
using SparseArrays

struct HalfEdge
    twin::UInt32
    next::UInt32
end

twin(he::HalfEdge) = he.twin
next(he::HalfEdge) = he.next


struct SimplicialSurface{P <: AbstractPoint}
    vertices::Array{P, 1}
    #halfedges::Array{HalfEdge, 1} # not implemented yet
    edges::Array{SVector{2, UInt32}, 1}
    faces::Array{SVector{3, UInt32}, 1}
    edgetoindex::Dict{SVector{2, UInt32}, UInt32}
    facetoindex::Dict{SVector{3, UInt32}, UInt32}
    edgevertexadjacency::SparseMatrixCSC{UInt32, UInt32}
    faceedgeadjacency::SparseMatrixCSC{UInt32, UInt32}
end

vertices(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.vertices
#halfedges(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.halfedges
edges(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.edges
faces(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.faces
edgetoindex(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.edgetoindex
facetoindex(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.facetoindex
edgevertexadjacency(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.edgevertexadjacency
faceedgeadjacency(simplicsurf::SimplicialSurface{<: AbstractPoint}) = simplicsurf.faceedgeadjacency

struct SurfaceSubset
    vertices::Array{UInt32, 1}
    edges::Array{UInt32, 1}
    faces::Array{UInt32, 1}
end

vertices(subset::SurfaceSubset) = subset.vertices
edges(subset::SurfaceSubset) = subset.edges
faces(subset::SurfaceSubset) = subset.faces

function buildsubsetvector(subset::SurfaceSubset, surface::SimplicialSurface, subsetfunc)::Array{UInt32, 1}
    cursubset = subsetfunc(subset)
    numsubsettotal = length(subsetfunc(surface))
    subsetvector = zeros(UInt32, numsubsettotal)
    for i in 1:length(cursubset)
        subsetvector[cursubset[i]] = one(UInt32)
    end
    return subsetvector
end

buildvertexvector(subset::SurfaceSubset, surface::SimplicialSurface) = buildsubsetvector(subset, surface, vertices)
buildedgevector(subset::SurfaceSubset, surface::SimplicialSurface) = buildsubsetvector(subset, surface, edges)
buildfacevector(subset::SurfaceSubset, surface::SimplicialSurface) = buildsubsetvector(subset, surface, faces)

sv(x::Vararg{Integer}) = SVector{length(x), UInt32}(x...)
sv(x::NTuple{N, Integer}) where N = sv(x...)

sortsv(x::Vararg{Integer}) = sort(sv(x...))
sortsv(x::NTuple{N, Integer}) where N = sortsv(x...)

function star(subset::SurfaceSubset, surface::SimplicialSurface)::SurfaceSubset
    originalvertices = buildvertexvector(subset, surface)
    originaledges = buildedgevector(subset, surface)
    originalfaces = buildfacevector(subset, surface)
    starvertices = originalvertices
    staredges = subsetnormalize(originaledges + edgevertexadjacency(surface) * starvertices)
    starfaces = subsetnormalize(originalfaces + faceedgeadjacency(surface) * staredges)
    starsubset = SurfaceSubset(map(tononzeroindices, (starvertices, staredges, starfaces))...)
    return starsubset
end

function closure(subset::SurfaceSubset, surface::SimplicialSurface)::SurfaceSubset
    originalfaces = buildfacevector(subset, surface)
    originaledges = buildedgevector(subset, surface)
    originalvertices = buildvertexvector(subset, surface)
    closurefaces = originalfaces
    closureedges = subsetnormalize(originaledges + faceedgeadjacency(surface)' * closurefaces)
    closurevertices = subsetnormalize(originalvertices + edgevertexadjacency(surface)' * closureedges)
    closuresubset = SurfaceSubset(map(tononzeroindices, (closurevertices, closureedges, closurefaces))...)
    return closuresubset
end

function link(subset::SurfaceSubset, surface::SimplicialSurface)::SurfaceSubset
    starofclosure = star(closure(subset, surface), surface)
    closureofstar = closure(star(subset, surface), surface)
    return subsetdifference(closureofstar, starofclosure, surface)
end

function subsetdifference(subset1::SurfaceSubset, subset2::SurfaceSubset, surface::SimplicialSurface)
    vertices1 = buildvertexvector(subset1, surface)
    edges1 = buildedgevector(subset1, surface)
    faces1 = buildfacevector(subset1, surface)

    vertices2 = buildvertexvector(subset2, surface)
    edges2 = buildedgevector(subset2, surface)
    faces2 = buildfacevector(subset2, surface)

    return SurfaceSubset(map(tononzeroindices ∘ subsetdifferencehelper, ((vertices1, vertices2), (edges1, edges2), (faces1, faces2)))...)
end

subsetnormalize(x::Array{UInt32, 1})::Array{UInt32, 1} = map(a->min(a, one(UInt32)), x)

function subsetdifferencehelper(arrs::Tuple{Array{UInt32, 1}, Array{UInt32, 1}})::Array{UInt32, 1}
    (x, y) = arrs
    lambdacompare((a, b)) = UInt32((a == one(UInt32)) && (b == zero(UInt32)))
    return map(lambdacompare, zip(x, y))
end

function tononzeroindices(x::Array{UInt32, 1})::Array{UInt32, 1}
    nonzeroindices = Array{UInt32, 1}(undef, sum(x))
    curindex = 1
    for i in 1:length(x)
        if x[i] == one(UInt32)
            nonzeroindices[curindex] = i
            curindex += 1
        end
    end
    return nonzeroindices
end


function constructsimplicialsurface(vertices::AbstractArray{P, 1}, tris::AbstractArray{NgonFace{3, T}})::SimplicialSurface{P} where {P <: AbstractPoint, T} 
    normedtris = map(NgonFace{3, UInt32}, tris)
    numtris = length(normedtris)

    #@show P
    faces = Array{SVector{3, UInt32}, 1}(undef, numtris)
    facetoindex = Dict{SVector{3, UInt32}, UInt32}()
    edges = Array{SVector{2, UInt32}, 1}(undef, 0)
    edgetoindex = Dict{SVector{2, UInt32}, UInt32}()
    numedges::UInt32 = zero(UInt32)
    edgevertexadjacencyCOO = Array{Tuple{UInt32, UInt32, UInt32}, 1}(undef, 0)
    faceedgeadjacencyCOO = Array{Tuple{UInt32, UInt32, UInt32}, 1}(undef, 0)
    #halfedgenext = Dict{UInt32, UInt32}()

    indiceshelper = ((1, 2), (2, 3), (1, 3))

    # Iterate over all faces, constructing indices and adjacency data as we go
    for i in 1:numtris
        # Construct face data
        triindcart = sortsv(normedtris[i]...)
        faces[i] = triindcart
        facetoindex[triindcart] = i

        # Iterate over edges in this face
        for j in 1:3
            # Construct edge data
            ind1 = triindcart[indiceshelper[j][1]]
            ind2 = triindcart[indiceshelper[j][2]]
            edgeindcart = sv(ind1, ind2) # don't need sortsv since already sorted from above
            #@show edgeindcart

            # New edge
            if !haskey(edgetoindex, edgeindcart)
                numedges += one(UInt32)
                #@show numedges
                push!(edges, edgeindcart)
                edgetoindex[edgeindcart] = numedges

                # Construct edge vertex adjacency matrix data
                push!(edgevertexadjacencyCOO, (numedges, ind1, one(UInt32)))
                push!(edgevertexadjacencyCOO, (numedges, ind2, one(UInt32)))
            end

            # Construct face edge adjacency matrix data
            edgeindlin = edgetoindex[edgeindcart]
            push!(faceedgeadjacencyCOO, (i, edgeindlin, one(UInt32)))
        end
    end

    #@show edgevertexadjacencyCOO
    #@show faceedgeadjacencyCOO

    edgevertexadjacency = sparse(unzip(edgevertexadjacencyCOO)...)
    #@show edgevertexadjacency
    faceedgeadjacency = sparse(unzip(faceedgeadjacencyCOO)...)
    #@show faceedgeadjacency

    #halfedges = Array{HalfEdge, 1}(undef, numedges * 2)
    #halfedges = Array{HalfEdge, 1}(undef, 0)


    #simplicsurf = SimplicialSurface{P}(vertices, halfedges, edges, faces, edgetoindex, facetoindex, edgevertexadjacency, faceedgeadjacency)
    simplicsurf = SimplicialSurface{P}(vertices, edges, faces, edgetoindex, facetoindex, edgevertexadjacency, faceedgeadjacency)

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
