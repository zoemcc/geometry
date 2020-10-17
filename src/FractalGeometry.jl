module FractalGeometry

using Meshing
using FileIO
using Roots
using MeshIO
using ForwardDiff
using GeometryBasics
import GeometryTypes
using LinearAlgebra
using CUDA
#using LightGraphs
import Makie
using SparseArrays
import GeometryBasics.faces
using StaticArrays
using NearestNeighbors
using Plots
using ForwardDiff
using DifferentialEquations

include("geom_example.jl")
include("simplicialsurface.jl")

end