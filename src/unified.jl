type BVPsol{T, T1, T2}
  node::Int
  npar::Int
  leftbc::Int
  npts::Int
  flag::Int
  mxnsub::Int
  xmsh::Vector{T}
  ymsh::Matrix{T1}
  parameters::Vector{T}
  iwork::Vector{Int}
  work::Vector{T2}
  solverparam::Dict{Symbol, Int}
  BVPsol{T, T1, T2} = new{T, T1, T2}()
end

# BVPinit!{}(node::Int, leftbc::Int, xint::Tuple{Real, Real}, y::Real, parameters_init::Vector{T}, mxnsub::Int, [solverparam::Dict])
# BVPinit!{}(node::Int, leftbc::Int, xint::Tuple{Real, Real}, y::Function, parameters_init::Vector{T}, mxnsub::Int)
# BVPinit!{}(node::Int, leftbc::Int, xmsh::Vector{T}, y::Real, parameters_init::Vector{T}, mxnsub::Int)
# BVPinit!{}(node::Int, leftbc::Int, xmsh::Vector{T}, y::Function, parameters_init::Vector{T}, mxnsub::Int)
# BVPinit!{}(node::Int, leftbc::Int, xmsh::Vector{T}, y::Real, mxnsub::Int)
# BVPinit!{}(node::Int, leftbc::Int, xmsh::Vector{T}, y::Function, mxnsub::Int)
# solve!(::BVPsol, ::Type{Val{:method}}, f::Function, bc::Function, df::Function, dbc::Function, do_trace::Bool)
# BVPeval{}(::BVPsol, x::Real, [yo::Vector{T}])
# BVPparam{}(::BVPsol) -> length(BVPsol.parameters) > 0 ? BVPsol.parameters : (warn("No parameters defined."); BVPsol.parameters)
