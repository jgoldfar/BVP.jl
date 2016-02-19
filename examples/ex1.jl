# Example 1 (MATLAB documentation for bvp4c)
include(joinpath("..", "src","bvp.jl"))
# using BVP
# Solve the example in the MATLAB documentation,
# y''+|y|=0,
# y(0)=0, y(4)=-2
function ex1(n::Int = 100)
  odefun{T<:Real}(x::T, y::Vector{T}) = [y[2], -abs(y[1])]
  bcfun{T<:Real}(ya::Vector{T}, yb::Vector{T}) = [ya[1], yb[1] + convert(T, 2)]
  initx = linspace(0, 4, n)
  # first solution:
  # inity = [10.0, 0.0]
  # second solution:
  inity = [1.0, 0.0]

  solinit = BVP.BVPSol(initx, inity)

  BVP.bvp1!(odefun, bcfun, solinit)
  # or
  BVP.bvp11!(odefun, bcfun, solinit)
end
# Uncomment to plot solution
# using Gadfly
# plot(x = solinit.x, y = solinit.y[1, :], Geom.line)
