using Sundials
function fem(p::Function, q::Function, f::Function, a2::Real = 1.0, a3::Real = 1.0, b2::Real = 1.0, b3::Real = 1.0, N::Int = 100)
  # Solve u'' + q u' + p u == f
  # a1 u(0) +a2 u'(0) = a3
  # b1 u(1) +b2 u'(1) = b3
  int_q(a::Real, b::Real) = first(quadgk(q, a, b))
  int_p(a::Real, b::Real) = first(quadgk(p, a, b))
  int_f(a::Real, b::Real) = first(quadgk(f, a, b))
  function int_tq(a::Real, b::Real)
    tq(t) = (t-a)*q(t)
    return first(quadgk(tq, a, b))
  end
  function int_tp(a::Real, b::Real)
    tp(t) = (t-a)*p(t)
    return first(quadgk(tp, a, b))
  end
  function int_tf(a::Real, b::Real)
    tf(t) = (t-a)*f(t)
    return first(quadgk(tf, a, b))
  end
  function int_tsqp(a::Real, b::Real)
    tsqp(t) = (t-a)*(t-a)*p(t)
    return first(quadgk(tsqp, a, b))
  end
  const xgrid = linspace(0,1,N)
  const x1 = xgrid[2]
  function fs(y, fy)
    const y1 = y[1]
    xc, xn = x1, x1
    y1m = (a3 - a2*y1)
    fy[1] = y1 * y1 * (int_tq(0, x1) - x1 + int_tsqp(0, x1)) -
            y1 * int_tf(0, x1) - y1m * (int_f(0, x1) + (a3/a2)) +0
            # y1 * y1m * (int_q(0, x1) + 2*int_tp(0, x1)) + 
            # y1m * y1m * (int_p(0, x1) + 1/a2)
    accum = xn * y1
    for i in 2:(N-1)
      yc = y[i]
      xc, xn = xn, xgrid[i + 1]
      hk = xn - xc
      fy[i] = -yc * int_tf(xc, xn) - accum * int_f(xc, xn) +
              yc * yc * (-hk + int_tq(xc, xn) + int_tsqp(xc, xn)) + yc * accum * (2*int_tp(xc, xn) + int_q(xc, xn)) +
              + accum * accum * int_p(xc, xn)
      accum += yc * hk
    end
    const yn = y[N]
    const hnm = xn - xc
    const ynm = (b3 - (b2 - hnm) * yn)
    fy[N] = -yn * int_tf(xc, xn) - ynm * int_f(xc, xn) +
            yn * yn * (-hnm + int_tq(xc, xn) + int_tsqp(xc, xn)) +  
            yn * ynm * (int_q(xc, xn) + 2*int_tp(xc, xn)) +
            ynm * ynm * int_p(xc, xn) + yn * (b3 - b2 * yn)
    return 0
  end
  res = Sundials.kinsol(fs, ones(N))
  function eval_f(t::Real)
    if t >= 1
      accum = 0
      for i in 1:(N - 1)
        accum += res[i] * (xgrid[i+1] - xgrid[i])
      end
      return res[N] * (t-xgrid[N]) + accum
    else
    const j = findfirst(j->(t >= xgrid[j]) && (t < xgrid[j + 1]), 1:(N-1))
    if j < 1
      d = a3 / res[1] - a2
      return res[1] * (t + d)
    else
      accum = 0
      for i in 1:(j - 1)
        accum += res[i] * (xgrid[i+1] - xgrid[i])
      end
      return res[j] * (t-xgrid[j]) + accum
    end
    end
    return zero(eltype(res))
  end
  return xgrid, map(eval_f, xgrid)
end

using Base.Test
const fem_test_tol = 1e-4
function fem_test1(N::Int = 200)
  p(x) = zero(x)
  r(x) = -one(x)
  f(x) = zero(x)
  const a2 = 1.0
  const a3 = 1.0
  const b2 = 1.0
  const b3 = sin(1.0) + cos(1.0)
  utrue(x) = sin(x)
  const xg, uc = fem(p, r, f, a2, a3, b2, b3, N)
  const ut = map(utrue, xg)
  println(norm(ut - uc))
  # @test (norm(ut - uc) < fem_test_tol)
  return 0
end
fem_test1()

# Previous eval_f
# function eval_f(t::Real)
#   if t >= 1
#     println("Extrapolating")
#   else
#   const j = findfirst(j->(t >= xgrid[j]) && (t < xgrid[j + 1]), 1:(N-1))
#   if j < 1
#     println("In first interval.")
#     d = a3 - a2 * res[1]
#   else
#     println("In interval ", j)
#   end
#   end
#   return 0
# end
# eval_f(-xgrid[2])
# eval_f(xgrid[1])
# eval_f(xgrid[1] + (xgrid[2] - xgrid[1]) / 2)
# eval_f(xgrid[2])
# eval_f(xgrid[2] + (xgrid[3]-xgrid[2])/2)
# eval_f(xgrid[end - 1])
# eval_f(xgrid[end])
# eval_f(xgrid[end] + (xgrid[end]-xgrid[end-1])/2)