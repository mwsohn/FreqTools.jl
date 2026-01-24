using Test, RDatasets, FreqTools, FreqTables

lbw = dataset("COUNT", "lbw")

@testset "One-Way Frequency Tables" begin

    # simple table
    o = tab(lbw, :Low)
    @test isapprox(o.omat[:,1], [130, 59, 189])
    @test isapprox(o.omat[:,2], [68.78, 31.22, 100.00], atol = 0.01)
    @test isapprox(o.omat[:,3], [68.78, 100.00, 100.00], atol = 0.01)

    # summarize
    o = tab(lbw, :Smoke, summarize = :BWt)
    @test isapprox(o.omat[:,1], [115, 74, 189 ])
    @test isapprox(o.omat[:, 2], [3054.96, 2772.30, 2944.29], atol=0.01)
    @test isapprox(o.omat[:, 3], [752.41, 659.81, 729.02], atol=0.01)

    # sort + summarize
    o = tab(lbw, :Race, summarize=:BWt, sort = true)
    @test isapprox(o.omat[:, 1], [115, 74, 189])
    @test isapprox(o.omat[:, 2], [3054.96, 2772.30, 2944.29], atol=0.01)
    @test isapprox(o.omat[:, 3], [752.41, 659.81, 729.02], atol=0.01)


end