# Implementation of FDM and collocation methods for BVP
function fdm(p::Function, r::Function, q::Function, va::Real, vb::Real, N::Int)
  # Solver for y'' - p y' - r y == q, y(0)=va, y(1)=vb
  const T = typeof(va)
  #TODO: Replace strict version comparison
  if VERSION == v"0.4.0-rc1"
    const xgrid = linspace(0, 1, N + 1)
  else
    const xgrid = linrange(0, 1, N + 1)
  end
  const h = step(xgrid)
  const subdiag, maindiag, supdiag = Array(T, N), Array(T, N + 1), Array(T, N)
  const A = Tridiagonal(subdiag, maindiag, supdiag)
  const rhs = Array(T, N + 1)

  maindiag[1] = one(T)
  supdiag[1] = zero(T)
  rhs[1] = va
  for i in 2:N
    xi = xgrid[i]
    vtmp = (h / 2) * p(xi)
    subdiag[i - 1] = 1 + vtmp
    maindiag[i] = -2 - h * h * r(xi)
    supdiag[i] = 1 - vtmp

    rhs[i] = h * h * q(xi)
  end
  subdiag[N] = zero(T)
  maindiag[N + 1] = one(T)
  rhs[N + 1] = vb

  return xgrid, A \ rhs
end

include("blockmat.jl")

function trapz{T<:Real}(A::Function, q::Function, Ba::Matrix{T}, Bb::Matrix{T}, beta::Vector{T}, xgrid::GridType{T})
  const n, n2 = size(Ba)
  @assert (n == n2) && (size(Bb) == (n, n2))
  const nx = length(xgrid)
  xn = xgrid[1]
  qn = q(xn)
  An = A(xn)
  const B = BlockMatrix{typeof(Ba)}(2nx)
  const E = eye(n)
  const rhs = Array(T, n * nx)
  # Set bc
  cind = 1
  B.rows[cind] = nx; B.cols[cind] = 1; B.mats[cind] = Ba
  cind += 1
  B.rows[cind] = nx; B.cols[cind] = nx; B.mats[cind] = Bb
  cind += 1
  for i in 1:(nx - 1)
    xc, xn = xn, xgrid[i + 1]
    qc, qn = qn, q(xn)
    Ac, An = An, A(xn)
    hi = xn - xc
    hir = 1 / hi
    B.rows[cind] = i; B.cols[cind] = i
    B.mats[cind] = -E * hir - Ac / 2
    cind += 1
    B.rows[cind] = i; B.cols[cind] = i + 1
    B.mats[cind] = E * hir - An / 2
    cind += 1
    rhs[((i - 1) * n + 1):(i * n)] = (qn + qc) / 2
  end
  rhs[((nx - 1) * n) + 1:(n * nx)] = beta
  # Debugging only:
  #     S = sparse(B)
  #     println([full(S) rhs])
  const yv = sparse(B) \ rhs
  return xgrid, reshape(yv, n, nx)
end
function trapz{T<:Real}(A::Function, q::Function, Ba::Real, Bb::Real, beta::Real, xgrid::GridType{T})
  const nx = length(xgrid)
  xn = xgrid[1]
  qn = q(xn)
  An = A(xn)
  const To = typeof(qn)
  const M = zeros(T, nx, nx)
  const rhs = Array(To, nx)

  for i in 1:(nx - 1)
    xc, xn = xn, xgrid[i + 1]
    qc, qn = qn, q(xn)
    Ac, An = An, A(xn)
    hi = xn - xc
    hir = 1 / hi
    M[i, i + 1] = hir - Ac / 2
    M[i, i] = -hir - An / 2
    rhs[i] = (qn + qc) / 2
  end
  M[nx, 1] = Ba
  M[nx, nx] = Bb
  rhs[nx] = beta
  return xgrid, M \ rhs
end
function trapz_sep{T<:Real}(A::Function, q::Function, Ba::Matrix{T}, Bb::Matrix{T}, betaa::Vector{T}, betab::Vector{T}, xgrid::GridType{T})
  const n, n2 = size(Ba)
  @assert (n == n2) && (size(Bb) == (n, n2))
  const nx = length(xgrid)
  xn = xgrid[1]
  qn = q(xn)
  An = A(xn)
  #   const B = BlockMatrix{typeof(Ba)}(2nx)
  const B = BlockMatrix{typeof(Ba)}()
  const E = eye(n)
  const rhs = Array(T, n * nx)
  # Set bc
  #   cind = 1
  #   B.rows[cind] = 1; B.cols[cind] = 1; B.mats[cind] = Ba
  push!(B.rows, 1); push!(B.cols, 1); push!(B.mats, Ba)
  #   cind += 1
  #   B.rows[cind] = nx; B.cols[cind] = nx; B.mats[cind] = Bb
  push!(B.rows, nx); push!(B.cols, nx); push!(B.mats, Bb)
  #   cind += 1
  rhs[1:n] = betaa

  for i in 2:(nx - 1)
    xc, xn = xn, xgrid[i + 1]
    qc, qn = qn, q(xn)
    Ac, An = An, A(xn)
    hi = xn - xc
    hir = 1 / hi
    #     B.rows[cind] = i; B.cols[cind] = i
    #     B.mats[cind] = -E * hir - Ac / 2
    push!(B.rows, i); push!(B.cols, i)
    push!(B.mats, -E * hir - Ac / 2)
    #     cind += 1
    push!(B.rows, i); push!(B.cols, i + 1)
    push!(B.mats, E * hir - An / 2)
    #     cind += 1
    rhs[((i - 1) * n + 1):(i * n)] = (qn + qc) / 2
  end
  rhs[((nx - 1) * n) + 1:(n * nx)] = betab
  # Debugging only:
  #   S = sparse(B)
  #   println("\n", [full(S) rhs])
  #   println("\n", full(S))
  #   const yv = sparse(S) \ rhs

  const yv = sparse(B) \ rhs
  return xgrid, reshape(yv, n, nx)
end

function trapz_sep{T<:Real}(A::Function, q::Function, betaa::Real, betab::Real, xgrid::GridType{T})
  const nx = length(xgrid)
  xn = xgrid[1]
  qn = q(xn)
  An = A(xn)
  const To = typeof(qn)
  const subdiag, maindiag, supdiag = zeros(To, nx - 1), Array(To, nx), Array(To, nx - 1)
  const M = Tridiagonal(subdiag, maindiag, supdiag)
  const rhs = Array(To, nx)
  maindiag[1] = one(To)
  rhs[1] = betaa
  supdiag[1] = zero(To)
  for i in 2:(nx - 1)
    xc, xn = xn, xgrid[i + 1]
    qc, qn = qn, q(xn)
    Ac, An = An, A(xn)
    hi = xn - xc
    hir = 1 / hi
    maindiag[i] = -hir - Ac / 2
    supdiag[i] = hir - An / 2
    rhs[i] = (qn + qc) / 2
  end
  rhs[nx] = betab
  maindiag[nx] = one(To)
  return xgrid, M \ rhs
end
