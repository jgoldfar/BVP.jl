function fem_spectral(r::Function, q::Function, va::Real, vb::Real, N::Int, ::Type{Val{:neumann}}, napprox::Int = N)
  # Solver for y'' - r y == q, y'(0)=va, y'(1)=vb
  T = typeof(va)
  xgrid = linspace(0, 1, N)
  A = Array(T, napprox, napprox)
  b = Array(T, napprox)
  for i in 1:napprox
    hi(x::Real) = cos(i * pi * x) + sin(i * pi * x)
    hip(x::Real) = i *pi * (cos(i * pi * x) - sin(i * pi * x))
    for j in 1:napprox
      hj(x::Real) = cos(j * pi * x) + sin(j * pi * x)
      hjp(x::Real) = j *pi * (cos(j * pi * x) - sin(j * pi * x))
      t = (x) -> hip(x) * hjp(x) + r(x) * hi(x) * hj(x)# + p(x) * hip(x) * hj(x) + r(x) * hi(x) * hj(x)
      A[j, i] = first(quadgk(t, 0, 1, abstol=1e-5))
    end
    t(x::Real) = q(x) * hi(x)
    b[i] = vb * hi(1) - va - first(quadgk(t, 0, 1, abstol=1e-5))
  end
#   println(A)
#   println(b)
  beta_i = \(A, b)
  rec(x::Real) = sinser_reconstruct(x, beta_i)
  return xgrid, map(rec, xgrid)
end
function sinser_reconstruct{T<:Real}(x::Real, beta_i::Vector{T})
  N = length(beta_i)
  vo = 0.0
  for i in 1:(N-1)
    vo += beta_i[i] * sin(i*pi*x)
  end
  return vo
end
function fem_galerk(p::Function, r::Function, q::Function, va::Real, vb::Real, N::Int, ::Type{Val{:neumann}})
  # Solver for y'' - p y' - r y == q, y'(0)=va, y'(1)=vb
  T = typeof(va)
  xgrid = linspace(0, 1, N)
  A = Array(T, N-1, N-1)
  b = Array(T, N-1)
  dx = xgrid[2]
  xni = 0.0
  for i in 1:(N-1)
    xci, xni = xni, xgrid[i]
    hi(x) = elem_1d_pwlin(x, xci-dx, xni+dx)
    hip(x) = elem_1d_pwlin_der(x, xci-dx, xni+dx)

    xnj = 0.0
    for j in 1:(N-1)
      xcj, xnj = xnj, xgrid[j]
      hj(x) = elem_1d_pwlin(x, xcj-dx, xnj+dx)
      hjp(x) = elem_1d_pwlin_der(x, xcj-dx, xnj+dx)
      A[j, i] = first(quadgk(x->(hip(x) * hjp(x) + p(x) * hip(x) * hj(x) + r(x) * hi(x) * hj(x)), 0, 1))
    end
    b[i] = vb * hi(1) - va * hi(0) - first(quadgk(x->q(x) * hi(x), 0, 1))
  end
  #     println(A)
  beta_i = \(A, b)
  rec(x) = elem_1d_pwlin_reconstruct(x, xgrid, beta_i)
  return xgrid, map(rec, xgrid)
end
function elem_1d_pwlin_reconstruct{T<:Real}(x::Real, xgrid::AbstractVector, beta_i::Vector{T})
  N = length(xgrid)
  dx = xgrid[2]
  xni = 0.0
  vo = 0.0
  for i in 1:(N-1)
    xci, xni = xni, xgrid[i]
    vo += beta_i[i] * elem_1d_pwlin(x, xci-dx, xni+dx)
  end
  return vo
end
function elem_1d_pwlin(x::Real, a::Real, b::Real)
  # Tent function with maximum 1 at a+b/2
  if x <= a || x >= b || isapprox(a, b)
    return zero(x)
  else
    return 1-2abs(x-(a+b)/2)/(b-a)
  end
end
# using Gadfly
# xg = linspace(0,1)
# yv = map(x->elem_1d_pwlin(x, 0.25, 0.75), xg)
# draw(PNG("elem.png", 12cm, 6cm), plot(x=xg, y=yv, Geom.point))
function elem_1d_pwlin_der(x::Real, a::Real, b::Real)
  # Ptwise a.e. derivative of tent function with maximum 1 at a+b/2
  if x <= a || x >= b || isapprox(a, b)
    return zero(x)
  else
    return 1-2sign(x-(a+b)/2)/(b-a)
  end
end
