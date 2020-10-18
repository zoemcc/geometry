using FractalGeometry

using StaticArrays
using GeometryBasics
using FileIO
import Makie

begin
    testp = Point3(0.6, 0.3, -0.4)
    const tau = 2pi
    const τ = tau
    const TAU = τ
    ori = Point3(-1.1, -1.1, -1.1)
    wid = Point3(2.2, 2.2, 2.2)
    orisv = SVector(ori)
    widsv = SVector(wid)

    outerbox(v) = FractalGeometry.boxdist(v, Point3(0.8, 0.8, 0.8))
    outerboxsdf = FractalGeometry.FunctionSDF(outerbox)
end

begin
    arc_length_sampling_freq = 0.01
    linecurveγ = FractalGeometry.FunctionCurve(FractalGeometry.linecurve, (-1., 1.))
    linecurvesdf = FractalGeometry.curve_to_kdtree_sdf(linecurveγ, 0.125, arc_length_sampling_freq)
    distancep = linecurvesdf(testp)

    linecurveboxsdf = FractalGeometry.IntersectSDF(linecurvesdf, outerboxsdf)
end

begin
    using FractalGeometry
    arc_length_sampling_freq = 0.01
    helixcurve1(t) = Point3(0.3*cos(tau * t), 0.3*sin(tau * t), t)
    helixcurve2(t) = Point3(-0.3*cos(tau * t), -0.3*sin(tau * t), t)
    helixcurveconn1(t) = helixcurve1(0.4) * t + helixcurve2(0.4) * (1 - t)
    helixcurveconn2(t) = helixcurve1(-0.4) * t + helixcurve2(-0.4) * (1 - t)

    curves = [helixcurve1, helixcurve2, helixcurveconn1, helixcurveconn2]
    tspans = [(-1., 1.), (-1., 1.), (0., 1.), (0., 1.)]
    extrusions = [0.125, 0.125, 0.05, 0.05]

    union_sdf = FractalGeometry.curves_to_union_sdf(curves, tspans, extrusions, arc_length_sampling_freq; boundary_sdf=outerboxsdf)

    distancep = union_sdf(testp)

end


begin
    curdist = union_sdf
    curh0 = 0.03

    result = FractalGeometry.distmesh(curdist, FractalGeometry.HUniform(), curh0, FractalGeometry.DistMeshSetup(); origin=ori, widths=wid, maxiters=100)
    resultouterbox = FractalGeometry.distmesh(outerbox, HUniform(), 0.5, DistMeshSetup(); origin=ori, widths=wid, maxiters=100)
    boundarytris, alltris = FractalGeometry.boundary_triangles(result.tetrahedra)
    boundarytrisouterbox, alltrisouterbox = FractalGeometry.boundary_triangles(resultouterbox.tetrahedra)

    boundarymesh = GeometryBasics.Mesh(result.points, boundarytris)
    boundarymeshouterbox = GeometryBasics.Mesh(resultouterbox.points, boundarytrisouterbox)
end
begin

    show_axis = true

    scene = Makie.Scene()
    Makie.mesh!(scene, boundarymesh, color=[norm(v[1:2]) for v in coordinates(boundarymesh)], show_axis=show_axis)

    Makie.wireframe!(scene, boundarymesh, show_axis=show_axis)
    Makie.wireframe!(scene, boundarymeshouterbox, show_axis=show_axis, color=RGB(0.8, 0.3, 0.8))
end

begin 
    savefile = "meshes/generated/doublehelixclipped.stl"
    save(savefile, boundarymesh)

end

begin 
    mesh = load(savefile)
    scene = Makie.Scene()
    Makie.mesh!(scene, mesh, color=[norm(v[1:2]) for v in coordinates(mesh)], show_axis=show_axis)

    Makie.wireframe!(scene, mesh, show_axis=show_axis)

end
