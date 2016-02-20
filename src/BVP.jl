VERSION >= v"0.4.0-dev+6521" && __precompile__()

module BVP
include("odesolve.jl")
include("shooting.jl")
include("fem.jl")
include("colloc.jl")
end #module
