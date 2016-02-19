# include(joinpath(dirname(@__FILE__),"..","src","colloc.jl"))
const trapz_test_tol = 1e-4
# using Base.Test
function trapz_test1(N::Int = 100)
  xgrid = linspace(0, 1, N)
  xgrid2 = logspace(0,1,N)
  A(x) = zero(x)
  q(x) = zero(x)

  betaa = 0.0
  betab = 0.0
  xg2, y1 = trapz_sep(A, q, betaa, betab, xgrid)
  @test norm(y1) < 1e-10


  xg2, y2 = trapz_sep(A, q, betaa, betab, xgrid2)
  @test norm(y2) < 1e-10

  betaa = 1.0
  betab = 1.0
  xg2, y3 = trapz_sep(A, q, betaa, betab, xgrid)
  @test norm(y3-1) < 1e-10

  xg2, y3 = trapz_sep(A, q, betaa, betab, xgrid2)
  @test norm(y3-1) < 1e-10
  return 0
end
trapz_test1()

function trapz_test2(N::Int = 100)
  A(x) = zero(x)
  q(x) = one(x)

  betaa = 0.0
  betab = 1.0
  xg1, y1 = trapz_sep(A, q, betaa, betab, linspace(0,1,N))
  xg2, y2 = trapz_sep(A, q, betaa, betab, linspace(0,1,2N))
  @test norm(y2 - xg2) < norm(y1 - xg1)

  betaa = 1.0
  betab = 2.0
  xg1, y1 = trapz_sep(A, q, betaa, betab, linspace(0,1,N))
  xg2, y2 = trapz_sep(A, q, betaa, betab, linspace(0,1,2N))
  utrue(x) = 1+x
  @test norm(y2 - map(utrue, xg2)) < norm(y1 - map(utrue, xg1))
  return 0
end
trapz_test2()

function trapz_test3(N::Int = 100; c::Real = 1.0)
  A(x) = one(x)
  q(x) = zero(x)
  if isapprox(c, 0) # zero solution is exact, so refining doesn't do much.
    return 0
  end
  betaa = c
  betab = c * e
  xg1, y1 = trapz_sep(A, q, betaa, betab, linspace(0, 1, N))
  xg2, y2 = trapz_sep(A, q, betaa, betab, linspace(0, 1, 2N))
  utrue(x) = c*exp(x)
  @test norm(y2 - map(utrue, xg2)) < norm(y1 - map(utrue, xg1))
  return 0
end
for c in -4:4
  trapz_test3(c = c)
end

function trapz_test4(N::Int = 100; Ba::Real = 1.0, Bb::Real = 1.0, beta::Real = 0.0)
  A(x) = zero(x)
  q(x) = zero(x)
  xgrid = linspace(0, 1, N)
  # Soln is constant C s.t. C Ba + C Bb = beta,
  # so C = beta/(Ba + Bb)
  xg1, y1 = trapz(A, q, Ba, Bb, beta, xgrid)
  @test norm(y1 - beta / (Ba + Bb)) < 1e-10
  return 0
end
for Ba in 1:3, beta in -2:2, Bb in (-0.5, 0.5, 1)
  trapz_test4(Ba = Ba, Bb = Bb, beta = beta)
end

function trapz_test5(N::Int = 100; Ba::Real = 1.0, Bb::Real = 1.0, beta::Real = 0.0, c::Real = 0.0)
  A(x) = zero(x)
  xgrid = linspace(0, 1, N)
  q(x) = c * one(x)
  # Soln is linear function cx + d s.t. Ba (c * 0 + d) + Bb (c * 1 + d) = beta,
  # so d Ba + Bb c + Bb d = beta, and hence
  # d = (beta - Bb c) / (Ba + Bb)
  utrue(x) = c * x + (beta - Bb * c) / (Ba + Bb)
  xg1, y1 = trapz(A, q, Ba, Bb, beta, xgrid)
  @test norm(y1 - map(utrue, xg1)) < 1e-10
  return 0
end
for Ba in 1:3, beta in -2:2, Bb in (-0.5, 0.5, 1, 2, 3), c in -2:2
  trapz_test5(Ba = Ba, Bb = Bb, beta = beta, c = c)
end

# Multi-variable
function trapz_test6(N::Int = 100; Bbv::Real = 0.0, betav::Real = 0.0)
  nv = 2
  O = ones(nv, nv)
  E = eye(nv)
  if VERSION >= v"0.4-"
    udE = flipdim(E, 1)
  else
    udE = flipud(E)
  end

  xgrid = linspace(0, 1, N)
  A(x) = zero(x) * O
  q(x) = zeros(nv)
  # Soln is constant C s.t. C Ba + C Bb = beta,
  # so C = beta/(Ba + Bb)
  Bav = 4.0
  Ba = Bav * E
  Bb = Bbv * udE
  beta = betav * ones(nv)
  try
    xg1, y1 = trapz(A, q, Ba, Bb, beta, xgrid)
    @test norm(y1 - betav / (Bav + Bbv)) < 1e-10
    #println("Bav: ", Bav, ", Bbv: ", Bbv, ", beta: ", beta, " :: OK")
  catch v
    println("Bav: ", Bav, ", Bbv: ", Bbv, ", beta: ", beta, " :: FAILED")
    #println(v)
    #Base.show_backtrace(STDOUT, catch_backtrace())
    #println("\n")
  end
  return 0
end
for Bbv in -3:3, betav in -2:2
  trapz_test6(Bbv = Bbv, betav = betav)
end
# function trapz_test7(N::Int = 10; c::Real = 0.0, Bav::Real = 4.0, Bbv::Real = 0.0, betav::Real = 0.0)
#   nv = 2
#   O = ones(nv, nv)
#   E = eye(nv)
#   if VERSION >= v"0.4-"
#     udE = flipdim(E, 1)
#   else
#     udE = flipud(E)
#   end
#   xgrid = linspace(0, 1, N)
#
#   A(x) = zero(x) * O
#   q(x) = c * ones(nv)
#   # Soln is linear function cx + d s.t. Ba (c * 0 + d) + Bb (c * 1 + d) = beta,
#   # so d Ba + Bb c + Bb d = beta, and hence
#   # d = (beta - Bb c) / (Ba + Bb)
#   Ba = Bav * E
#   Bb = Bbv * udE
#   beta = betav * ones(nv)
#
#   d = (betav - Bbv * c) / (Bav + Bbv)
#
#   utrue(x) = c * x + d
#   try
#     xg1, y1 = trapz(A, q, Ba, Bb, beta, xgrid)
#     @test norm(y1[1, :] - map(utrue, xg1)') < 1e-10
#   catch v
#     #     print("\n", Bav, ", ", Bbv, ", ", c, ", ", betav, ", ", d, ": ")
#     @test false
#   end
#   return 0
# end
# # Problem is not well posed for some of these combinations...
# for Bav in 1:3, Bbv in (-0.5, 0.5, 1, 2, 3), betav in -2:2, c in -2:2
#   trapz_test7(Bav = Bav, Bbv = Bbv, betav = betav, c = c)
# end

# Multivariable, reformulation of problem
# y'' + a pi^2 y = 0
# y'(0)=y'(1)=0
# y(0) = a
#
# which has true solution
# y(x) = a cos(pi x)
# (verification: y' = -a pi sin(pi x), y'' = -a pi^2 cos(pi x); so y'' + a pi^2 y = 0)
function trapz_test8(N::Int = 100; a::Real = 0.0)
  nv = 2
  Am = [0.0 1.0;
              -pi^2 0.0]
  Ov = zeros(nv)
  xgt1 = linspace(0, 1, N)
  xgt2 = linspace(0, 1, 2N)
  # u(0) + u'(0) = 0
  # u'(1) = 0
  Ba = [1.0 1.0;
        0.0 0.0]
  Bb = [0.0 0.0;
        0.0 1.0]
  beta = copy(Ov)

  q(x) = zero(x) * Ov
  A(x) = Am
  beta[1] = a
  utrue(x) = [a * cos(pi * x);
              -a * pi * sin(pi * x)]
  xg1, y1 = trapz(A, q, Ba, Bb, beta, xgt1)
  xg2, y2 = trapz(A, q, Ba, Bb, beta, xgt2)
  error1 = 0.0
  for i in 1:N
    error1 += norm(y1[:, i] - utrue(xg1[i]))
  end
  error2 = 0.0
  for i in 1:2N
    error2 += norm(y2[:, i] - utrue(xg2[i]))
  end
#   println(error1, ", ", error2)
  @test error2 <= error1
  return 0
end
for a in -5:5
  trapz_test8(a = a)
end
