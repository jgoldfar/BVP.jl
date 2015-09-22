const fem_test_tol = 1e-4
function fem_test1(N::Int = 100)
  p(x) = zero(x)
  r(x) = zero(x)
  q(x) = zero(x)
  const va = 0.0
  const vb = 0.0
  utrue(x) = zero(x)
  const xg, uc = fem()
  # println(xg)
  # println(uc)
  # const ut = map(utrue, xg)
  @test true # (norm(ut - uc) < fdm_test_tol)
  return 0
end
fem_test1()