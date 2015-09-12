# odesolve.jl: Solvers I've written myself
## Begin ODE solvers
include("fwdeuler.jl")
#TODO: Implement Runge-Kutta ODE Solver method
function rk{T<:Real}(odefun::Function,
                     xgrid::GridType{T},
                     iv::Vector{T},
                     order::Int,
                     param::Real = NaN)
#   const yvout = zeros(T, length(iv), length(xgrid))
  error("rk method is unimplemented.")
  return yvout
end


## Sundials wrapper
# With parameter
# function sundials_wrap{T<:Real}(odefun::Function,
#                                 xgrid::GridType{T},
#                                 iv::Vector{T},
#                                 param::T)
#   function tfunc{T1<:Real}(t::T1, y::Vector{T}, ydot::Vector{T})
#     ydot[:] = odefun(t, y, param)
#   end
#   res = Sundials.cvode(tfunc, iv, xgrid)
#   return res
# end

# # Without parameter
# function sundials_wrap{T<:Real}(odefun::Function,
#                                 xgrid::GridType{T},
#                                 iv::Vector{T})
#   function tfunc{T1<:Real}(t::T1, y::Vector{T}, ydot::Vector{T})
#     ydot[:] = odefun(t, y)
#   end
#   res = Sundials.cvode(tfunc, iv, xgrid)
#   return res
# end
