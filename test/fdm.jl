# include(joinpath(dirname(@__FILE__),"..","src","colloc.jl"))
const fdm_test_tol = 1e-4
function fdm_test1(N::Int = 100)
  p(x) = zero(x)
  r(x) = zero(x)
  q(x) = zero(x)
  const va = 0.0
  const vb = 0.0
  utrue(x) = zero(x)
  const xg, uc = fdm(p, r, q, va, vb, N)
  const ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test1()
function fdm_test2(N::Int = 100)
  p(x) = zero(x)
  r(x) = zero(x)
  q(x) = zero(x)
  const va = 0.0
  const vb = 1.0
  utrue(x) = x
  const xg, uc = fdm(p, r, q, va, vb, N)
  const ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test2()
function fdm_test3(N::Int = 100)
  p(x) = zero(x)
  r(x) = zero(x)
  q(x) = exp(x)
  const va = 1.0
  const vb = e
  utrue(x) = exp(x)
  const xg, uc = fdm(p, r, q, va, vb, N)
  const ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test3()
function fdm_test4(N::Int = 200)
  p(x) = zero(x)
  r(x) = -one(x)
  q(x) = zero(x)
  const va = 0.0
  const vb = sin(1.0)
  utrue(x) = sin(x)
  const xg, uc = fdm(p, r, q, va, vb, N)
  const ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test4()
function fdm_test5(N::Int = 200)
  p(x) = -one(x)
  r(x) = zero(x)
  q(x) = zero(x)
  const va = 1.0
  const vb = 0.0
  utrue(x) = (e*exp(-x) - 1)/(e-1)
  const xg, uc = fdm(p, r, q, va, vb, N)
  const ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test5()
