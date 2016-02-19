## Begin RK4
function rk4{T<:Real}(odefun::Function,
                      xgrid::AbstractVector,
                      iv::Vector{T},
                      param::Real)
  nv = length(iv)
  nxgrid = length(xgrid)
  yvout = Array(T, nv, nxgrid)
  for i in 1:nv
    yvout[i, 1] = iv[i]
  end
  xcurr = xgrid[1]
  yprev = copy(iv)
  for i in 2:nxgrid
    # Save x values/indexing between iterations
    xprev, xcurr = xcurr, xgrid[i]
    dx = xcurr - xprev
    dxh = dx / 2
    y1 = odefun(xprev, yprev)
    y2 = odefun(xprev + dxh, yprev + y1 * dxh)
    y3 = odefun(xprev + dxh, yprev + y2 * dxh)
    y4 = odefun(xprev + dx, yprev + y3 * dx)
    # RK4 step
    ycurr = yprev + dx/6 * (y1 + y2 + y3 + y4)

    # Set into output array
    for j in 1:nv
      yvout[j, i] = ycurr[j]
      yprev[j] = ycurr[j]
    end
  end
  return yvout
end
function rk4{T<:Real}(odefun::Function,
                      xgrid::AbstractVector,
                      iv::Vector{T})
  nv = length(iv)
  nxgrid = length(xgrid)
  yvout = Array(T, nv, nxgrid)
  for i in 1:nv
    yvout[i, 1] = iv[i]
  end
  xcurr = xgrid[1]
  yprev = iv
  for i in 2:nxgrid
    # Save x values/indexing between iterations
    xprev, xcurr = xcurr, xgrid[i]
    dx = xcurr - xprev
    dxh = dx / 2
    y1 = odefun(xprev, yprev)
    y2 = odefun(xprev + dxh, yprev + y1 * dxh)
    y3 = odefun(xprev + dxh, yprev + y2 * dxh)
    y4 = odefun(xprev + dx, yprev + y3 * dx)
    # RK4 step
    ycurr = yprev + dx/6 * (y1 + y2 + y3 + y4)

    # Set into output array
    for j in 1:nv
      yvout[j, i] = ycurr[j]
      yprev[j] = ycurr[j]
    end
  end
  return yvout
end

# Implementations saving only the value at the final moment
function rk4_final{T<:Real}(odefun::Function,
                            xgrid::AbstractVector,
                            iv::Vector{T},
                            param::Real)
  nxgrid = length(xgrid)
  nv = length(iv)
  xcurr = xgrid[1]
  ycurr = copy(iv)
  for i in 2:nxgrid
    # Save x values/indexing between iterations
    xprev, xcurr = xcurr, xgrid[i]
    dx = xcurr - xprev

    dxh = dx / 2
    y1 = odefun(xprev, ycurr)
    y2 = odefun(xprev + dxh, ycurr + y1 * dxh)
    y3 = odefun(xprev + dxh, ycurr + y2 * dxh)
    y4 = odefun(xprev + dx, ycurr + y3 * dx)
    # RK4 step
    for j in 1:nv
      ycurr[j] += dx/6 * (y1[j] + y2[j] + y3[j] + y4[j])
    end
  end
  return ycurr
end
function rk4_final{T<:Real}(odefun::Function,
                            xgrid::AbstractVector,
                            iv::Vector{T})
  nxgrid = length(xgrid)
  nv = length(iv)
  xcurr = xgrid[1]
  ycurr = copy(iv)
  for i in 2:nxgrid
    # Save x values/indexing between iterations
    xprev, xcurr = xcurr, xgrid[i]
    dx = xcurr - xprev

    dxh = dx / 2
    y1 = odefun(xprev, ycurr)
    y2 = odefun(xprev + dxh, ycurr + y1 * dxh)
    y3 = odefun(xprev + dxh, ycurr + y2 * dxh)
    y4 = odefun(xprev + dx, ycurr + y3 * dx)
    # RK4 step
    for j in 1:nv
      ycurr[j] += dx/6 * (y1[j] + y2[j] + y3[j] + y4[j])
    end
  end
  return ycurr
end
## End RK4 method
