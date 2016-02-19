type BlockMatrix{T<:AbstractMatrix}
  mats::Vector{T}
  rows::Vector{Int}
  cols::Vector{Int}
  BlockMatrix() = new(Array(T, 0), Array(Int, 0), Array(Int, 0))
  BlockMatrix(n::Int) = new(Array(T, n), Array(Int, n), Array(Int, n))
end
#
# The two implementations of Block to sparse
# are currently interchangeable for small block
# sizes, but the first outperforms the second by
# ~an order of magnitude at block size 100.
function Base.sparse{T}(B::BlockMatrix{T})
  # Assume all blocks are nxn. TODO: Remove assumption.
  M = B.mats[1]
  n, n2 = size(M)
  vT = eltype(M)
  nm = length(B.mats)
  Iv, Jv = Array(Int, nm*n*n), Array(Int, nm*n*n)
  Ev = Array(vT, nm*n*n)
  for i in 1:nm
    Mi = B.mats[i]
    r = (B.rows[i] - 1) * n + 1
    c = (B.cols[i] - 1) * n + 1
    for j in 1:n
      for k in 1:n
        ind = sub2ind((nm, n, n), i, j, k)
        Iv[ind] = j - 1 + r
        Jv[ind] = k - 1 + c
        Ev[ind] = Mi[j, k]
      end
    end
  end
  sparse(Iv, Jv, Ev)
end
# function Base.sparse{T}(B::BlockMatrix{T})
#   # Assume all blocks are nxn. TODO: Remove assumption.
#   M = B.mats[1]
#   (n, n2) = size(M)
#   vT = eltype(M)
#   Iv, Jv = Array(Int, 0), Array(Int, 0)
#   Ev = Array(vT, 0)
#   for (i, M) in enumerate(B.mats)
#     r = (B.rows[i] - 1) * n + 1
#     c = (B.cols[i] - 1) * n + 1
#     for i in 1:n
#       for j in 1:n
#         push!(Iv, i - 1 + r)
#         push!(Jv, j - 1 + c)
#         push!(Ev, M[i, j])
#       end
#     end
#   end
#   sparse(Iv, Jv, Ev)
# end
function Base.show{T}(io::IO, B::BlockMatrix{T})
  println(io, "BlockMatrix{", T, "}: [", max(B.rows...), "x", max(B.cols...), "] blocks.")
end
dump(B::BlockMatrix) = show(STDOUT, B)
# Planned interface:
# push!(B::BlockMatrix{T}, (Block::T, row::Int, col::Int))
# pop!(B::BlockMatrix{T}) -> Block::T
# eltype(B::BlockMatrix{T})
# size(B::BlockMatrix{T}, [n::Int])
# length(B::BlockMatrix{T})
# copy(B::BlockMatrix{T})
# show(B::BlockMatrix{T})
# (maybe?) getindex(B::BlockMatrix{T}, i_1,...)
# nnz(B::BlockMatrix{T})
# sparse(B::BlockMatrix{T})
# full(B::BlockMatrix{T})
# eachindex(B::BlockMatrix{T})
