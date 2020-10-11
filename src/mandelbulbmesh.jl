using Makie
using DistMesh
using GeometryBasics

include("geom_example.jl")
include("simplicialsurface.jl")

mandeldist = p->-mandelbulbdistfunc(1)(p)

ori = GeometryBasics.Point{3, Float64}(-2, -2, -2)
wid = GeometryBasics.Point{3, Float64}(4, 4, 4)

result = distmesh(mandeldist, HUniform(), 0.2, DistMeshSetup(); origin=ori, widths=wid)
#result = distmesh(spheredist, HUniform(), 0.1)

#Makie.scatter(result.points)

function boundary_triangles(tetrahedra::Array{SimplexFace{4, Int32}, 1})
    tritonum = Dict{SVector{3, UInt32}, Int}()
    @show tritonum
    tritonum[sv(1, 2, 3)] = 1
    @show tritonum

    tricombinations = [sv(1,2,3), sv(1,2,4), sv(1,3,4), sv(2,3,4)]
    @show tricombinations
    @show tricombinations[1]
    for (i, tet) in enumerate(tetrahedra)
        for combination in tricombinations
            index = sortsv(tet[combination]...)
            if !haskey(tritonum, index)
                tritonum[index] = 1
            else
                tritonum[index] += 1
            end
        end
    end
    @show length(tritonum)

    boundarytris = Array{NgonFace{3, UInt32}, 1}(undef, 0)
    for (key, val) in pairs(tritonum)
        if val == one(UInt32)
            tri = NgonFace{3, UInt32}(key...)
            push!(boundarytris, tri)
        end
    end

    return boundarytris
end

boundarytris = boundary_triangles(result.tetrahedra)

boundarymesh = GeometryBasics.Mesh(result.points, boundarytris)

scene = Makie.Scene()
Makie.mesh!(scene, boundarymesh, color=[norm(v[1:2]) for v in coordinates(boundarymesh)])

Makie.wireframe!(scene, boundarymesh)

