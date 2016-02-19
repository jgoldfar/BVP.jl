# Recommended: Run with  `find . -name '*.jl' | peat 'julia --color test/runtests.jl'`

using Compat
const pkgbasedir = joinpath(dirname(@__FILE__), "..")
# include(joinpath(dirname(@__FILE__), "..", "src", "bvp.jl"))
include(joinpath(pkgbasedir, "src", "odesolve.jl"))
include(joinpath(pkgbasedir, "src", "colloc.jl"))
include(joinpath(pkgbasedir, "src", "shooting.jl"))
include(joinpath(pkgbasedir, "src", "fem.jl"))

if VERSION >= v"0.5-"
  using Base.Test
else
  using BaseTestNext
  const Test = BaseTestNext
end

tic()
@testset "All" begin
  @testset "MATLAB Compare" begin
  include("mlcompare.jl")
  @test true
  end
  @testset "FDM" begin
  include("fdm.jl")
  end
  @testset "Trapz" begin
  include("trapz.jl")
  end
  @testset "Shooting" begin
  include("shooting.jl")
  end
  @testset "FEM" begin
  include("fem.jl")
  end
end
toc()
