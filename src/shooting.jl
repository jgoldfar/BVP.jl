function sht_single{T<:Real}(A::Function, q::Function, Ba::Matrix{T}, Bb::Matrix{T}, beta::Vector{T}, xgrid::GridType{T})
  const nbeta, nv = size(Ba)
  const iv_test = zeros(nv)
  const Yb = Array(T, nv, nv)
  odefun1(x, y) = A(x) * y
  for i in 1:nv
    iv_test[i] = 1.0
    tmp = fwdeuler_final(odefun1, xgrid, iv_test)
    for j in 1:nv
      Yb[j, i] = tmp[j]
    end
    iv_test[i] = 0.0
  end
  const Q = Ba + Bb * Yb

  iv_test[1] = 1.0 # Avoid catching (only) zero solution.
  odefun2(x, y) = A(x) * y + q(x)
  const vb = fwdeuler_final(odefun2, xgrid, iv_test)

  const betahat = beta - Ba * iv_test - Bb * vb
  const ivt = iv_test + vec(\(Q, betahat))
  return xgrid, fwdeuler(odefun2, xgrid, ivt)
end

#TODO: Consider AD for calculation of ODE jacobian and BC Jacobian
function sht_single{T<:Real}(odefun::Function, odejacfun::Function, bcfun::Function, bcjacfun::Function, nv::Int, nc::Int, xgrid::GridType{T}, nlsolver::Val{:newton_internal})
  # Solve the problem
  # y' = odefun(x, y)
  # bcfun(y(a), y(b)) = 0 (nc*1)
  # where y is a nv-dimensional vector.
  # The inputs assumed to be of the form
  # odefun :: R x R^(nv) -> R^(nv)
  # odejacfun :: R x R^(nv) -> R^((1+nv) * nv)
  # bcfun :: R^(nv) x R^(nv) -> R^(nc)
  # bcjacfun :: R^(nv) x R^(nv) -> R^(2nv * nc)

  const iv_curr = zeros(nv)
  const iv_prev = copy(iv_curr)
  error("Unimplemented.")
end

# Default solver is internal
sht_single{T<:Real}(odefun::Function, odejacfun::Function, bcfun::Function, bcjacfun::Function, nv::Int, nc::Int, xgrid::GridType{T}) = sht_single(odefun, odejacfun, bcfun, bcjacfun, nv, nc, xgrid, Val{:newton_internal})
