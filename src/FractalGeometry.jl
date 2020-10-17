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

include("geom_example.jl")
include("simplicialsurface.jl")

end