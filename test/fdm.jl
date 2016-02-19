# include(joinpath(dirname(@__FILE__),"..","src","colloc.jl"))
const fdm_test_tol = 1e-4
function fdm_test1(N::Int = 100)
  p(x) = zero(x)
  r(x) = zero(x)
  q(x) = zero(x)
  va = 0.0
  vb = 0.0
  utrue(x) = zero(x)
  xg, uc = fdm(p, r, q, va, vb, N)
  ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test1()
function fdm_test2(N::Int = 100)
  p(x) = zero(x)
  r(x) = zero(x)
  q(x) = zero(x)
  va = 0.0
  vb = 1.0
  utrue(x) = x
  xg, uc = fdm(p, r, q, va, vb, N)
  ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test2()
function fdm_test3(N::Int = 100)
  p(x) = zero(x)
  r(x) = zero(x)
  q(x) = exp(x)
  va = 1.0
  vb = e
  utrue(x) = exp(x)
  xg, uc = fdm(p, r, q, va, vb, N)
  ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test3()
function fdm_test4(N::Int = 200)
  p(x) = zero(x)
  r(x) = -one(x)
  q(x) = zero(x)
  va = 0.0
  vb = sin(1.0)
  utrue(x) = sin(x)
  xg, uc = fdm(p, r, q, va, vb, N)
  ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test4()
function fdm_test5(N::Int = 200)
  p(x) = -one(x)
  r(x) = zero(x)
  q(x) = zero(x)
  va = 1.0
  vb = 0.0
  utrue(x) = (e*exp(-x) - 1)/(e-1)
  xg, uc = fdm(p, r, q, va, vb, N)
  ut = map(utrue, xg)
  @test (norm(ut - uc) < fdm_test_tol)
  return 0
end
fdm_test5()
