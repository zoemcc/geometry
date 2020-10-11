import Makie
using DistMesh
using GeometryBasics

include("geom_example.jl")
include("simplicialsurface.jl")

mandeldist = p->-mandelbulbdistfunc(1)(p)

ori = GeometryBasics.Point{3, Float64}(-1.2, -1.2, -1.2)
wid = GeometryBasics.Point{3, Float64}(2.4, 2.4, 2.4)
#ori = GeometryBasics.Point{3, Float64}(-1.5, -1.5, -1.5)
#wid = GeometryBasics.Point{3, Float64}(3.0, 3.0, 3.0)

#result = distmesh(mandeldist, HUniform(), 0.04, DistMeshSetup(); origin=ori, widths=wid, maxiters=100)
#result = distmesh(spheredist, HUniform(), 0.1, DistMeshSetup(); origin=ori, widths=wid, maxiters=100)
function curveofshapes(curve::AbstractArray{Tuple{Point3{Float64}, Float64}, 1}, shapefunc::Function)
    function innerdistance(v)
        distances = map(p->shapefunc((v - p[1]) .* (1 / p[2])), curve)
        return minimum(distances)
    end
    return innerdistance
end
boxdistcur(v) = boxdist(v, Point3{Float64}(0.2, 0.2, 0.2))
curve = Array{Tuple{Point3{Float64}, Float64}, 1}(undef, 0)
push!(curve, (Point3{Float64}(-0.4, -0.4, -0.4), 1.0))
push!(curve, (Point3{Float64}(0.4, -0.4, -0.4), 0.7))
push!(curve, (Point3{Float64}(0.4, 0.4, -0.4), 0.4))
push!(curve, (Point3{Float64}(0.4, 0.4, 0.4), 0.3))
multibox = curveofshapes(curve, boxdistcur)
result = distmesh(multibox, HUniform(), 0.07, DistMeshSetup(); origin=ori, widths=wid, maxiters=100)

#Makie.scatter(result.points)

function boundary_triangles(tetrahedra::Array{SimplexFace{4, Int32}, 1})
    tritonum = Dict{SVector{3, UInt32}, Int}()
    @show tritonum
    tritonum[sv(1, 2, 3)] = 1
    @show tritonum

    tricombinations = [sv(1,2,3), sv(1,2,4), sv(1,3,4), sv(2,3,4)]
    @show tricombinations
    @show tricombinations[1]
    alltris = Array{NgonFace{3, UInt32}, 1}(undef, 0)
    for (i, tet) in enumerate(tetrahedra)
        for combination in tricombinations
            index = sortsv(tet[combination]...)
            if !haskey(tritonum, index)
                tritonum[index] = 1
            else
                tritonum[index] += 1
            end
            tri = NgonFace{3, UInt32}(index...)
            push!(alltris, tri)
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

    return boundarytris, alltris
end

boundarytris, alltris = boundary_triangles(result.tetrahedra)

boundarymesh = GeometryBasics.Mesh(result.points, boundarytris)
vertexsets = connectedvertices(faces(boundarytris))
removedmesh = removesmallcomponents(boundarymesh, vertexsets; sizetokeep=2)
allmesh = GeometryBasics.Mesh(result.points, alltris)

scene = Makie.Scene()
Makie.mesh!(scene, boundarymesh, color=[norm(v[1:2]) for v in coordinates(boundarymesh)])

Makie.wireframe!(scene, boundarymesh)
Makie.wireframe!(scene, allmesh)
Makie.wireframe!(scene, removedmesh)
Makie.mesh!(scene, removedmesh, color=[norm(v[1:2]) for v in coordinates(removedmesh)])

save("meshes/generated/reallynicespherebigg.stl", boundarymesh)
mesh = load("meshes/generated/reallynicespherebigg.stl")
scene = Makie.Scene()
Makie.mesh!(scene, mesh, color=[norm(v[1:2]) for v in coordinates(mesh)])
Makie.wireframe!(scene, mesh)
