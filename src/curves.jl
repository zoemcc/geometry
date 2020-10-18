abstract type AbstractCurve <: Function end

struct FunctionCurve <: AbstractCurve
    γ::Function
    tspan::Tuple{Float64, Float64}
end

function (γ::FunctionCurve)(t::Real)
    γ.γ(t)
end

