# Example 1 (Harmonic eivenvalue problem)
include(joinpath("..", "src","bvp.jl"))
# Solve
# y''+\lambda y=0,
# y(0)=0, y(2*pi)=0, y'(0)=1.0 (analytic solution should be sin(x))

function harmonic(n::Int = 1000)
  odefun{T<:Real}(x::T, y::Vector{T}, param::T) = [y[2], -param * y[1]]
  bcfun{T<:Real}(ya::Vector{T}, yb::Vector{T}, param::T) = [ya[1], ya[2] - convert(T, 1), yb[1]]
  initx = linspace(0, pi, n)
  inity = [2.0, 0.0]
  initparam = 3.0

  solinit = BVP.BVPSol(initx, inity, hasparam = true, initparam = initparam,
                       xreltolerance = 1e-9, fabstolerance = 1e-9)

  BVP.bvp1!(odefun, bcfun, solinit)
  # or
  BVP.bvp11!(odefun, bcfun, solinit)

  # View resulting solution structure
  return solinit
end
sol = harmonic()
dump(solinit)
# Make plot of solution
using Gadfly
plot(layer(x = solinit.x, y = solinit.y[1, :], Geom.line),
     layer(x = solinit.x, y = [sin(x) for x in solinit.x], Geom.line),
     )

# Convergence to the solution is qualitatively ok with bvp1, but finds the "wrong eigenvalue."
# evidently, this requires a different solver.
