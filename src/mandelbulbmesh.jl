import Makie
using ColorTypes
using RegionTrees
using AdaptiveDistanceFields
using DistMesh
using Zygote
using GeometryBasics

include("geom_example.jl")
include("simplicialsurface.jl")

testp = Point3{Float64}(1., -1., 0.64)
tau = 2pi

mandeldist = p->-mandelbulbdistfunc(1)(p)

ori = Point3{Float64}(-1.1, -1.1, -1.1)
wid = Point3{Float64}(2.2, 2.2, 2.2)
orisv = SVector(ori)
widsv = SVector(wid)
#ori = GeometryBasics.Point{3, Float64}(-1.5, -1.5, -1.5)
#wid = GeometryBasics.Point{3, Float64}(3.0, 3.0, 3.0)

#result = distmesh(mandeldist, HUniform(), 0.04, DistMeshSetup(); origin=ori, widths=wid, maxiters=100)
#result = distmesh(spheredist, HUniform(), 0.1, DistMeshSetup(); origin=ori, widths=wid, maxiters=100)

outerbox(v) = boxdist(v, Point3{Float64}(0.8, 0.8, 0.8))

boxdistcur(v) = boxdist(v, Point3{Float64}(0.2, 0.2, 0.2))

curve = Array{Tuple{Point3{Float64}, Float64}, 1}(undef, 0)
push!(curve, (Point3{Float64}(-0.4, -0.4, -0.4), 1.0))
push!(curve, (Point3{Float64}(0.4, -0.4, -0.4), 0.7))
push!(curve, (Point3{Float64}(0.4, 0.4, -0.4), 0.45))
push!(curve, (Point3{Float64}(0.4, 0.4, 0.4), 0.3))
multibox = curveofshapes(curve, boxdistcur)

linecurveextrudedist(v) = curvedistgen(linecurve, -10., 10.)(v) - 0.4

linecurveextrude_intersect_box_dist(v) = max(outerbox(v), linecurveextrudedist(v))

#=
function centraldiff(f::Function,p::VT) where VT
    println("Using our central diff")
    @show f, VT
    deps = sqrt(eps(eltype(VT)))
    dx = (f(p.+VT(deps,0,0)) - f(p.-VT(deps,0,0)))
    dy = (f(p.+VT(0,deps,0)) - f(p.-VT(0,deps,0)))
    dz = (f(p.+VT(0,0,deps)) - f(p.-VT(0,0,deps)))
    grad = VT(dx,dy,dz)./(2deps) #normalize?

end
=#

helixcurve(t) = Point3{Number}(0.3*cos(tau * t), 0.3*sin(tau * t), t)
#helixcurve(t) = Point3{Float64}(0., 0., 1) .* t
#helixcurvedist = curvedistgen(helixcurve, -10., 10.)
#disttocurve(t) = normsq(helixcurve(t) - testp)
#Ddisttocurve = D(disttocurve)
helixcurveextrudedist(v) = curvedistgen(helixcurve, -10., 10.)(v) - 0.2
helixcurveextrude_intersect_box_dist(v) = max(outerbox(v), helixcurveextrudedist(v))

# v expensive!!!
#helixadf = AdaptiveDistanceField(helixcurveextrude_intersect_box_dist, orisv, widsv)


weirdcurve(t) = Point3{Number}(0.5*cos(tau * t), 0.5*sin(tau * t), 0.4*cos(2tau*t))
weirdcurveextrudedist(v) = curvedistgen(weirdcurve, -1.1, 1.1)(v) - 0.125
weirdcurveextrude_intersect_box_dist(v) = max(outerbox(v), weirdcurveextrudedist(v))


begin
    curdist = helixadf
    curh0 = 0.015

    fh(v) = 0.05 + (0.5 * curdist(v))
    result = distmesh(curdist, fh, curh0, DistMeshSetup(); origin=ori, widths=wid, maxiters=100)
    resultouterbox = distmesh(outerbox, HUniform(), 0.5, DistMeshSetup(); origin=ori, widths=wid, maxiters=100)
    #resultouterbox = distmesh(outerbox, HUniform(), 5curh0, DistMeshSetup(); origin=ori, widths=wid, maxiters=2)

    #Makie.scatter(result.points)

    boundarytris, alltris = boundary_triangles(result.tetrahedra)
    boundarytrisouterbox, alltrisouterbox = boundary_triangles(resultouterbox.tetrahedra)

    boundarymesh = GeometryBasics.Mesh(result.points, boundarytris)
    boundarymeshouterbox = GeometryBasics.Mesh(resultouterbox.points, boundarytrisouterbox)
    #vertexsets = connectedvertices(faces(boundarytris))
    #removedmesh = removesmallcomponents(boundarymesh, vertexsets; sizetokeep=2)
    #allmesh = GeometryBasics.Mesh(result.points, alltris)

    show_axis = true

    scene = Makie.Scene()
    Makie.mesh!(scene, boundarymesh, color=[norm(v[1:2]) for v in coordinates(boundarymesh)], show_axis=show_axis)

    Makie.wireframe!(scene, boundarymesh, show_axis=show_axis)
    Makie.wireframe!(scene, boundarymeshouterbox, show_axis=show_axis, color=RGB(0.8, 0.3, 0.8))
end

#Makie.wireframe!(scene, allmesh)
#Makie.wireframe!(scene, removedmesh)
#Makie.mesh!(scene, removedmesh, color=[norm(v[1:2]) for v in coordinates(removedmesh)])
aleaf = nothing
outv = nothing
for leaf in allleaves(adf.root)
    aleaf = leaf
    outv = RegionTrees.center(leaf)
    #outv = hcat(collect(RegionTrees.vertices(leaf.boundary))...)
    
end

adfpoints = map(RegionTrees.center, allleaves(helixadf.root))

scene = Makie.Scene()
Makie.scatter!(scene, adfpoints, size=1e-4, color=[norm(v) for v in adfpoints])


savename = "meshes/generated/superhighreshelix.stl"
save(savename, boundarymesh)
mesh = load(savename)
scene = Makie.Scene()
Makie.mesh!(scene, mesh, color=[norm(v[1:2]) for v in coordinates(mesh)])
Makie.wireframe!(scene, mesh)
