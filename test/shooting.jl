# typealias GridType{T} Union(Vector{T}, LinSpace{T})
# include(joinpath(dirname(@__FILE__), "..", "src", "odesolve.jl"))
# include(joinpath(dirname(@__FILE__), "..", "src", "shooting.jl"))
# using Base.Test
const sht_test_tol = 1e-4
function sht_test1(;N::Int = 100, b::Real = 1)
  A(x) = zero(x)
  q(x) = zero(x)
  xgrid = linspace(0, b, N)
  xg1, y1 = sht_single(A, q,
                       [1.0 1.0; 0.0 0.0],
                       [0.0 0.0; 0.0 1.0],
                       [0.0, 0.0],
                       xgrid)
  @test norm(y1) < 1e-10
  return 0
end
for b in 1:10
  sht_test1(b=b)
end

# Ascher et. a. example 4.2
function sht_test2(;N::Int = 50, lambda::Real = 1, b::Real = 1)
  xgrid1 = linspace(0, b, N)
  xgrid2 = linspace(0, b, 2N)
  A(x) = [0 1.0;
          lambda^2 0.0]
  q(x) = (1 - lambda^2) * [0, exp(x)]

  xg1, y1 = sht_single(A, q,
                       [1.0 0.0;
                        0.0 0.0],
                       [0.0 0.0;
                        1.0 0.0],
                       [1.0, exp(b)],
                       xgrid1)
  xg2, y2 = sht_single(A, q,
                       [1.0 0.0;
                        0.0 0.0],
                       [0.0 0.0;
                        1.0 0.0],
                       [1.0, exp(b)],
                       xgrid2)
  utrue(x) = [exp(x), exp(x)]
  error1 = 0.0
  error2 = 0.0
  for (i, x) in enumerate(xg1)
    error1 += norm(y1[:, i] - utrue(x))
  end
  for (i, x) in enumerate(xg2)
    error2 += norm(y2[:, i] - utrue(x))
  end
  #     println(error1, ", ", error2)
  @test error2 < error1
  return 0
end
for lambda in 0.5:5
  for b in 1:10
    sht_test2(lambda=lambda, b=b)
  end
end

function sht_test3(N::Int = 50; a::Real = 1.0)
  nv = 2
  # u = a cos(pi x)
  # u' = -a pi sin(pi x)
  # u'' = -a pi^2 cos(pi x) = -pi^2 u
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

  xg1, y1 = sht_single(A, q, Ba, Bb, beta, xgt1)
  xg2, y2 = sht_single(A, q, Ba, Bb, beta, xgt2)
  error1 = 0.0
  for i in 1:N
    error1 += norm(y1[:, i] - utrue(xg1[i]))
  end
  error2 = 0.0
  for i in 1:2N
    error2 += norm(y2[:, i] - utrue(xg2[i]))
  end
  #   println(error1, ", ", error2)
  @test error2 < error1
  return 0
end
for a in 1:3
  sht_test3(a = a)
end
