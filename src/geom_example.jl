using Meshing
using FileIO
using MeshIO
using GeometryBasics
import GeometryBasics
import GeometryTypes
using LinearAlgebra
using CUDA
#using LightGraphs
import Makie
using SparseArrays


defaultsamplesregion() = ((60, 60, 60), Rect(Vec(-2,-2,-2), Vec(4,4,4)))
samples, region = defaultsamplesregion()

function gyroid()
    gyroid(v) = cos(v[1])*sin(v[2])+cos(v[2])*sin(v[3])+cos(v[3])*sin(v[1])
    gyroid_shell(v) = max(gyroid(v)-0.4,-gyroid(v)-0.4)


    # generate directly using GeometryBasics API
    # Rect specifies the sampling intervals
    gy_mesh = Mesh(gyroid_shell, Rect(Vec(0,0,0),Vec(pi*4,pi*4,pi*4)),
                    MarchingCubes(), samples=(50,50,50))


    return gy_mesh
end

spheredist(v) = norm(v) - 1.0

function boxdist(v, b)
    q = abs.(v) - b
    return norm(max.(q, 0.0)) + min.(maximum(q), 0.0)
end


function sphere()


    # generate directly using GeometryBasics API
    # Rect specifies the sampling intervals
    sphere_mesh = Mesh(spheredist, Rect(Vec(-1.2,-1.2,-1.2),Vec(2.4, 2.4, 2.4)),
                    MarchingCubes(), samples=(10,10,10))


    return sphere_mesh
end


function rectGBtoGT(rect::GeometryBasics.HyperRectangle{N, T}) where {N, T <: Number}
    return GeometryTypes.HyperRectangle{N, T}(origin(rect), widths(rect))
end

function rectGTtoGB(rect::GeometryTypes.HyperRectangle{N, T}) where {N, T <: Number}
    return GeometryBasics.HyperRectangle{N, T}(rect.origin, rect.widths)
end

function origin(rect::GeometryTypes.HyperRectangle{N, T}) where {N, T <: Number}
    return rect.origin
end

function widths(rect::GeometryTypes.HyperRectangle{N, T}) where {N, T <: Number}
    return rect.widths
end


function spheresdf(samples::NTuple{3, T}, region::Rect) where {T <: Integer}
    regionsamples = sampleregion(samples, region)
    return spheredist.(regionsamples)
end

function sampleregion(samples::NTuple{3, T}, region::Rect) where {T <: Integer}
    regionsamples = zeros(Point3{Float64}, samples[1], samples[2], samples[3])
    regionwidths = widths(region)
    regionorigin = origin(region)
    regionincrements = regionwidths ./ (samples .- 1)
    for i = 1:samples[1], j = 1:samples[2], k = 1:samples[3]
        samplepositions = ((i, j, k) .- 1) .* regionincrements .+ regionorigin
        regionsamples[i, j, k] = Point3{Float64}(samplepositions...)
    end
    return regionsamples
end


function facearray(mesh)
    mesh_faces = faces(mesh)
    m = length(mesh_faces)
    face_arr = Array{NgonFace{3, Int64}, 1}(undef, m)
    for i = 1:m
        face_arr[i] = Int64.(mesh_faces[i])
    end
    return face_arr
end


function singletri()
    verts = Array{Point{3, Float64}, 1}(undef, 3)
    pointtype = Point{3, Float64}
    verts[1] = pointtype([0.0, 0.0, 1.0])
    verts[2] = pointtype([0.0, 1.0, 0.0])
    verts[3] = pointtype([1.0, 0.0, 0.0])
    faces = Array{GLTriangleFace, 1}(undef, 1)
    faces[1] = GLTriangleFace(1, 2, 3)
    geom_mesh = Mesh(verts, faces)
    return geom_mesh
end

function connectedvertices(faces::AbstractArray{NgonFace{3, N}}) where {N <: Integer}
    Ntype = eltype(N)
    connectedvertexsets = Vector{Set{Ntype}}(undef, 0)
    #@show connectedvertexsets
    #@show length(connectedvertexsets)
    for i =  1:length(faces)
    #for i =  1:10
        if i % 100 == 0
            @show i
        end
        #@show i
        curset = Set{Ntype}([Ntype(faces[i][v]) for v = 1:3])
        #@show curset
        #@show typeof(curset)
        for v = 1:3
            for k = length(connectedvertexsets):-1:1
            #@show k
                vert = Ntype(faces[i][v])
                if vert in connectedvertexsets[k]
                    #@show typeof(curset)
                    #@show connectedvertexsets[k]
                    vertexset = popat!(connectedvertexsets, k)
                    #@show vertexset
                    curset = union!(vertexset, curset)
                    #@show typeof(curset)
                end
            end
        end
        #@show typeof(connectedvertexsets)
        push!(connectedvertexsets, curset)
    end
    return connectedvertexsets
end

function removesmallcomponents(mesh, connectedvertexsets; sizetokeep=1)
    sizessorted = sort(length.(connectedvertexsets))
    componentsizetokeep = sizessorted[sizetokeep]
    biggestcomponent = collect(filter(x -> length(x) == componentsizetokeep, connectedvertexsets)[1])
    oldtonew = Dict{UInt32, UInt32}()
    numnewverts = length(biggestcomponent)
    newverts = Array{Point3{Float64}, 1}(undef, numnewverts)
    oldverts = coordinates(mesh)
    #for i = 1:4
    for i = 1:numnewverts
        oldindex = biggestcomponent[i]
        oldtonew[oldindex] = i
        #@show oldverts[oldindex]
        newverts[i] = oldverts[oldindex]
    end
    @show length(newverts)

    oldfaces = faces(mesh)
    newfaces = Array{NgonFace{3, UInt32}, 1}(undef, 0)
    @show newfaces
    for i = 1:length(oldfaces)
    #for i = 1:100
        #@show i
        #@show oldfaces[i]
        newindices = Array{UInt32, 1}(undef, 0)
        for j = 1:3
            oldindex = UInt32(oldfaces[i][j])
            #@show oldindex
            if haskey(oldtonew, oldindex)
                push!(newindices, oldtonew[oldindex])
            end
        end
        #@show newindices
        if length(newindices) == 3
            newface = NgonFace{3, UInt32}(newindices...)
            push!(newfaces, newface)
        end
        #newindices = [oldtonew[] for j = 1:3]
        #@show newindices
        #newface = NgonFace{3, UInt32}([])
        
    end
    #@show newfaces
    newmesh = GeometryBasics.Mesh(newverts, newfaces)
    return newmesh
    #return biggestcomponent

end


function mandelbulbmesh(samples::NTuple{3, T}, region::Rect, numfractaliter::Integer) where {T <: Integer}
    mandelbulbdist = mandelbulbdistfunc(numfractaliter)
    return Mesh(mandelbulbdist, region, samples, MarchingTetrahedra(), pointtype=Point3{Float64})
end



function meshview(inmesh, color)
    return Makie.mesh(inmesh, color=color, show_axis=false)
end

function meshview(inmesh)
    return Makie.mesh(inmesh, color=[norm(v) for v in coordinates(inmesh)], show_axis=false)
end


function normsq(p::AbstractArray{T}) where {T <: Number}
    return sum(p .^ 2)
end

function inversesqrt(num::T) where {T <: Number}
    return 1.0 / sqrt(num)
end

function mandelbulbdistfunc(numfractaliter::Integer) 
    function mandelbulbdist(p::Point{3, T}) where {T <: Real}
        px, py, pz = p[1], p[2], p[3]
        dz = 1.0
        w = Array{T, 1}(undef, 3)
        w .= 0

        stop_iter = 0

        w .+= p
        #@show typeof(w)
        m = normsq(w)
        #trap = vcat(abs.(w), m)

        for i = 1:numfractaliter

            w2 = w .* w
            w4 = w2 .* w2

            x, y, z = w[1], w[2], w[3]
            x2, y2, z2 = w2[1], w2[2], w2[3]
            x4, y4, z4 = w4[1], w4[2], w4[3]

            m2 = m * m
            m4 = m2 * m2

            dz *= 8.0 * sqrt(m4 * m2 * m)
            dz += 1.0

            k3 = x2 + z2
            k2 = inversesqrt(k3 * k3 * k3 * k3 * k3 * k3 * k3)
            k1 = x4 + y4 + z4 - 
                6.0 * y2 * z2 - 
                6.0 * x2 * y2 + 
                2.0 * z2 * x2
            k4 = x2 - y2 + z2

            neww = [
                px + 
                64.0 * x * y * z * (x2 - z2) * k4 * 
                (x4 - 6.0 * x2 * z2 + z4) * k1 * k2,
                py +
                -16.0 * y2 * k3 * k4 * k4 + k1 * k1,
                pz +
                -8.0 * y * k4 * (x4 * x4 - 
                28.0 * x4 * x2 * z2 +
                70.0 * x4 * z4 -
                28.0 * x2 * z2 * z4 +
                z4 * z4) * k1 * k2
                ]
            w[:] = neww[:]

            #trap[:, :, :] = min.(trap, vcat(abs.(w), m)) .* (1.0 .- stop_iter) .+ stop_iter .* trap
            m = normsq(w)
            if m > 256.0
                return -0.0005
            end
            #m[:, :, :] = min.(normsq(w), 256.0)
            #m[map(isnan, m)] .= 256.0
            #dz[map(isnan, dz)] .= 1.0

            #@show i
            #@show w
            #@show m
            #@show stop_iter
            #@show dz
            stop_iter = i

        end

        dist = 0.25 * log(m) * sqrt(m) / dz
        return dist
    end
    return mandelbulbdist
end

function init()
    blob2 = mandelbulbmesh(samples, region, 2)
    connecteds = connectedvertices(faces(blob2))
    #newblob2 = removesmallcomponents(blob2, connecteds)
    return blob2, connecteds #, newblob2
end

function loadbunny()
    bunnymeshpre = load("meshes/original/small_bunny.ply")
    vertspre, tris = coordinates(bunnymeshpre), faces(bunnymeshpre)
    transformmat = [1 0 0; 0 0 1; 0 1 0]
    verts = map(vert->Point3{Float32}(transformmat * vert), vertspre)

    return Mesh(verts, tris)
end

function bunnyeyepot(t::Float64)
    eyepot = Array(sparsevec(
        # eyeleft, eyeright, max size, booty, ear left, ear right
        [1311, 1015, 1430, 680, 568, 747],
         [Float64(3 * cos(t)), Float64(3 * cos(t)), Float64(0), Float64(-2), Float64(0.5), Float64(0.5)]))
    return eyepot
end


#begin
    #p1 = solvepoissonproblem(bs, bunnyeyepot(0.))
    #for i in 1:100
        #p1 = solveheatstep(bs, p1, 0.04)
        #meshview(bunny, p1)
    #end
#end

