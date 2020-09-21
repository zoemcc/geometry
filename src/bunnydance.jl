using ColorSchemes
include("geom_example.jl")
include("simplicialsurface.jl")

#= notes for the audiomp4:
starting beats:
0.05s : starts
[0.05s, 0.2s, 0.4s, 0.6s, 0.8s (big), 1.1s, 1.3s, 1.5s, 1.7s]
1.85s : beats end
2.1s-2.3s : RUN 
2.5s : the drop starts
every 0.2s roughly a beat
first drum style: [3.3s, 5.1s, 6.9ss, 8.6s, 10.2s]
second drum snare? [2.7s, 3s, 3.2s, 4.7s, 6.1s, 6.4s, 6.6s, 8.2s, 9.5s, 9.8s, 10s]
big waaamk [2.5s, 6s, 9.6s]
10.4s : the drop ends
=#

function makescene(scene)  
    #mesh = Makie.mesh!(scene, bunny, color=p1, show_axis=false, colormap=ColorSchemes.linear_bmw_5_95_c86_n256)
    #angle = -pi / 4 - pi / 10
    #quat = Makie.Quaternion((0., 0., sin(angle / 2), cos(angle / 2)))
    #Makie.rotate!(scene, quat)
    Makie.update_cam!(scene, 0.8 * Makie.Vec3f0(-1.4153199, 3.1283882, 0.1344321), Makie.Vec3f0(0., 0., -0.1))
    Makie.update!(scene)
    # 1.721643, 2.5823247, 0.6416913
    #Float32[-1.9936941, 3.4563544, 0.2008853]
    # realclose Float32[-1.7143533, 2.6653938, 0.2666605]
    # adjusted -1.4153199, 3.1283882, 0.1344321
    return scene
end

function bunnydance(scene)
    bunny = loadbunny()
    transvec = [Point3{Float32}(0., 0., 0.)]
    bunnyverts = coordinates(bunny)
    bunnyverts .+= transvec
    bs = constructsimplicialsurface(bunny)
    p1 = 3.0 .* bunnyeyepot2(0.)
    p2 = p1
    #p2 = solvepoissonproblem(bs, p1)

    #trans = Point3{Float32}(100., -100., 100.)
    #Makie.translate!(scene, trans...)
    #Makie.update!(scene)
    framerate = 60
    starttime = 0.5
    endtime = 10.5
    time = starttime:1/framerate:endtime
    
    numverts = length(vertices(bs))
    velocity = zeros(Float64, numverts)

    beatfreq = 2.0 * π * (1.0 / 1.6)
    beatfreqsmooth = 2.0 * π * 5.0
    forcingvalue = [cos(beatfreq * (t - 2.7)) for t in time]
    smoothingsteps = [cos(beatfreqsmooth * (t - 2.7)) for t in time]
    for i in 1:length(smoothingsteps)
        if smoothingsteps[i] < 0.0
            smoothingsteps[i] .* 0.1
        end
    end
    j = 0
    Makie.record(scene, "videos/generated/bunnydance.mp4", time; framerate = framerate) do i
        if mod(j, 10)  == 0
            @show minimum(p2)
            @show maximum(p2)
        end
        j += 1
        forcing = 1.0 .* forcingvalue[j].* bunnyeyepot2(0.)
        if i > 2.7
            p2[:], velocity[:] = solvewavestep(bs, p2, velocity, forcing, 100.0, 0.6)
            #@show forcing[1311]
            #@show norm(p2)
            #p2[:] = solveheatstep(bs, p2, 0.02)
            #p2[:], velocity[:] = solvewavestep(bs, p2, velocity, zeros(Float64, numverts), 1000.0, 0.1)
            newverts = solvecurvatureflowstep(bs, 0.00017 * smoothingsteps[j])
            updatemeshandsurface!(bunny, bs, newverts)
        end
        scene.plots[1][:mesh][] = bunny
        scene.plots[1][:color][] = p2
    end
    return scene
end

function bunnyeyepot2(t::Float64)
    eyepot = Array(sparsevec(
        # eyeleft, eyeright, max size, booty, ear left, ear right
        [1311, 1015, 1430, 680, 568, 747],
         [Float64(1 * cos(t)), Float64(1 * cos(t)), Float64(0), Float64(0), Float64(0), Float64(0)]))
    return eyepot
end


function generatesoundcurves(framerate=10)
    time = 0:1/framerate:10.5
    
    startingpulses = [0.1, 0.2, 0.4, 0.6, 0.8, 1.1, 1.3, 1.5, 1.7]
    minstart, maxstart = minimum(startingpulses), maximum(startingpulses)

    run = [2.1]

    mainbeat = collect(2.5:0.2:10.5)
    bigwaamk = [2.5, 6, 9.6]

    drums = [3.3, 5.1, 6.9, 8.6, 10.2]
    snaredrums = [2.7, 3, 3.2, 4.7, 6.1, 6.4, 6.6, 8.2, 9.5, 9.8, 10]

    heatstepcurve = zeros(Float64, length(time))
    for t in time
        @show t
    end


    return collect(time), heatstepcurve
    #return time, startingpulses, run, mainbeat, bigwaamk, drums, snaredrums, heatstepcurve
end

# bunny = loadbunny(); bs = constructsimplicialsurface(bunny); newverts = solvecurvatureflowstep(bs, 0.001); updatemeshandsurface!(bunny, bs, newverts); Makie.mesh(bunny, color=bunnyeyepot(0.), show_axis=false)

