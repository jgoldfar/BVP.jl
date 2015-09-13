# typealias GridType{T} Union(Vector{T}, LinSpace{T})
# include(joinpath(dirname(@__FILE__), "..", "src", "fem.jl"))
# using Base.Test
const fem_test_tol = 1e-4
function fem_test1(N::Int = 100)
  r(x) = zero(x)
  q(x) = zero(x)
  const va = 0.0
  const vb = 0.0
  utrue(x) = zero(x)
  const xg, uc = fem_spectral(r, q, va, vb, N, Val{:neumann})
  const ut = map(utrue, xg)
  @test norm(ut - uc) < fem_test_tol
  return 0
end
fem_test1(50)

function fem_test2(N::Int = 100)
  r(x) = zero(x)
  q(x) = zero(x)
  const va = 1.0
  const vb = 1.0
  utrue(x) = x
  const xg, uc = fem_spectral(r, q, va, vb, N, Val{:neumann})
#   return xg, uc
  const ut = map(utrue, xg)
  @test norm(ut - uc) < fem_test_tol
  return 0
end
# fem_test2(20)
# xg, uc = fem_test2(20)

# using Gadfly
# draw(PNG("elem.png", 12cm, 6cm), plot(x=xg, y=uc, Geom.point))
