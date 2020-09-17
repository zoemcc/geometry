include("geom_example.jl")
include("simplicialsurface.jl")

function bunnydance()
    bunny = loadbunny()
    p1 = bunnyeyepot(0.)
    bs = constructsimplicialsurface(bunny)
    p2 = solvepoissonproblem(bs, p1)
    scene = Makie.mesh(bunny, color=p2, show_axis=false)
    angle = -pi / 4 - pi / 10
    quat = Makie.Quaternion((0., 0., sin(angle / 2), cos(angle / 2)))
    Makie.rotate!(scene, quat)
    #trans = Point3{Float32}(100., -100., 100.)
    #Makie.translate!(scene, trans...)
    Makie.update_cam!(scene, Makie.Vec3f0(1.721643, 2.5823247, 0.6416913), Makie.Vec3f0(0., 0., -0.1))
    return scene
end

