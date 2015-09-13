## Begin forward Euler
function fwdeuler{T<:Real}(odefun::Function,
                           xgrid::GridType{T},
                           iv::Vector{T},
                           param::Real)
  const nv = length(iv)
  const nxgrid = length(xgrid)
  const yvout = Array(T, nv, nxgrid)
  for i in 1:nv
    yvout[i, 1] = iv[i]
  end
  xcurr = xgrid[1]
  const yprev = copy(iv)
  for i in 2:nxgrid
    # Save x values/indexing between iterations
    xprev, xcurr = xcurr, xgrid[i]
    dx = xcurr - xprev

    # Euler step
    ycurr = yprev + dx * odefun(xcurr, yprev, param)

    # Set into output array
    for j in 1:nv
      yvout[j, i] = ycurr[j]
      yprev[j] = ycurr[j]
    end
  end
  return yvout
end
function fwdeuler{T<:Real}(odefun::Function,
                           xgrid::GridType{T},
                           iv::Vector{T})
  const nv = length(iv)
  const nxgrid = length(xgrid)
  const yvout = Array(T, nv, nxgrid)
  for i in 1:nv
    yvout[i, 1] = iv[i]
  end
  xcurr = xgrid[1]
  const yprev = copy(iv)
  for i in 2:nxgrid
    # Save x values/indexing between iterations
    xprev, xcurr = xcurr, xgrid[i]
    dx = xcurr - xprev

    # Euler step
    ycurr = yprev + dx * odefun(xcurr, yprev)

    # Set into output array
    for j in 1:nv
      yvout[j, i] = ycurr[j]
      yprev[j] = ycurr[j]
    end
  end
  return yvout
end

# Implementations saving only the value at the final moment
function fwdeuler_final{T<:Real}(odefun::Function,
                                 xgrid::GridType{T},
                                 iv::Vector{T},
                                 param::Real)
  const nxgrid = length(xgrid)
  const nv = length(iv)
  xcurr = xgrid[1]
  const ycurr = copy(iv)
  for i in 2:nxgrid
    # Save x values/indexing between iterations
    xprev, xcurr = xcurr, xgrid[i]
    dx = xcurr - xprev

    # Euler step
    tmp = odefun(xcurr, ycurr, param)
    for j in 1:nv
      ycurr[j] += dx * tmp[j]
    end
  end
  return ycurr
end
function fwdeuler_final{T<:Real}(odefun::Function,
                                 xgrid::GridType{T},
                                 iv::Vector{T})
  const nxgrid = length(xgrid)
  const nv = length(iv)
  xcurr = xgrid[1]
  const ycurr = copy(iv)
  for i in 2:nxgrid
    # Save x values/indexing between iterations
    xprev, xcurr = xcurr, xgrid[i]
    dx = xcurr - xprev

    # Euler step
    tmp = odefun(xcurr, ycurr)
    for j in 1:nv
      ycurr[j] += dx * tmp[j]
    end
  end
  return ycurr
end
## End forward Euler method
