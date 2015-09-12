# Recommended: Run with  `find . -name '*.jl' | peat 'julia --color test/runtests.jl'`

if VERSION < v"0.4-"
  typealias LinSpace{T} Vector{T}
end
using Compat
typealias GridType{T} Union(Vector{T}, LinSpace{T})
# include(joinpath(dirname(@__FILE__), "..", "src", "bvp.jl"))
include(joinpath(dirname(@__FILE__), "..", "src", "odesolve.jl"))
include(joinpath(dirname(@__FILE__), "..", "src", "colloc.jl"))
include(joinpath(dirname(@__FILE__), "..", "src", "shooting.jl"))
include(joinpath(dirname(@__FILE__), "..", "src", "fem.jl"))

using Base.Test
test_handler(r::Test.Success) = print_with_color(:green, ".")
const rpt_fails = false
test_handler(r::Test.Failure) = (print_with_color(:red, "."); rpt_fails && println(STDERR, "\n", r))
test_handler(r::Test.Error) = (print_with_color(:red, "x"); println(STDERR, "\n", r))
const testfiles = (
  "mlcompare.jl",
  "fdm.jl",
  "trapz.jl",
  "shooting.jl",
  "fem.jl"
  )

tic()
Test.with_handler(test_handler) do
  for testfile in testfiles
    try
      println(" Running test $testfile.\n")
      include(joinpath(dirname(@__FILE__),testfile))
      println("\n")
    catch v
      println(" Error $v\n in $testfile\n")
      Base.show_backtrace(STDERR, catch_backtrace())
      print("\n")
    end
  end
end
toc()
