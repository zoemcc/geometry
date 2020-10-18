abstract type SignedDistanceField <: Function end

struct SDFKDTree <: SignedDistanceField
    tree::KDTree
    extrusion::Float64
end


function (sdf::SignedDistanceField)(point::AbstractPoint)
    @show "No sdf distance implemented for this type:", typeof(sdf), dump(sdf) 
end

@inline function (sdf::SDFKDTree)(point::AbstractPoint)
    distance = knn(sdf.tree, point, 1)[2][1] - sdf.extrusion
end

function curve_to_kdtree_sdf(γ::FunctionCurve, extrusion::Float64, arc_length_sampling_freq::Float64)
    t_to_s, s_to_t, γs, sspan = arc_length_parameterization(γ.γ, γ.tspan)
    mins, maxs = sspan[1:2]

    γsamples = (mins:arc_length_sampling_freq:maxs)
    γsampled_points = map(γs, γsamples)

    γkdtree = KDTree(γsampled_points)
    SDFKDTree(γkdtree, extrusion)
end

struct UnionSDF{SDF1 <: SignedDistanceField, SDF2 <: SignedDistanceField} <: SignedDistanceField
    sdf1::SDF1
    sdf2::SDF2
end

@inline function (sdf::UnionSDF)(point::AbstractPoint)
    distance = min(sdf.sdf1(point), sdf.sdf2(point))
end

struct IntersectSDF{SDF1 <: SignedDistanceField, SDF2 <: SignedDistanceField} <: SignedDistanceField
    sdf1::SDF1
    sdf2::SDF2
end

@inline function (sdf::IntersectSDF)(point::AbstractPoint)
    distance = max(sdf.sdf1(point), sdf.sdf2(point))
end

struct FunctionSDF <: SignedDistanceField
    sdf::Function
end

@inline function (sdf::FunctionSDF)(point::AbstractPoint)
    distance = sdf.sdf(point)
end

function curves_to_union_sdf(γs::Vector{Function}, tspans::Vector{Tuple{Float64, Float64}}, 
        extrusions::Vector{Float64}, arc_length_parameterization::Float64; boundary_sdf=nothing)
    γcurves = [FunctionCurve(γ, tspan) for (γ, tspan) in zip(γs, tspans)]
    γsdfs = [curve_to_kdtree_sdf(γ, extrusion, arc_length_parameterization) for (γ, extrusion) in zip(γcurves, extrusions)]
    unionsdf = reduce((sdf1, sdf2) -> UnionSDF(sdf1, sdf2), γsdfs)
    if !(boundary_sdf isa Nothing) && boundary_sdf isa SignedDistanceField
        sdf = IntersectSDF(unionsdf, boundary_sdf)
    else
        sdf = unionsdf
    end
end

