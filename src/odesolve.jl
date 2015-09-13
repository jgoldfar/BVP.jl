# odesolve.jl: Solvers I've written myself
## Begin ODE solvers
include("fwdeuler.jl")
include("rkode.jl")
include("odewrap.jl")


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
