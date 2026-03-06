using Test, RDatasets, FreqTools, FreqTables

lbw = dataset("COUNT", "lbw")
bfi = dataset("psych", "bfi")

@testset "One-Way Frequency Tables" begin

    # simple table
    o = tab(lbw, :Low)
    @test isapprox(o.omat[:, 1], [130, 59, 189])
    @test isapprox(o.omat[:, 2], [68.78, 31.22, 100.00], atol=0.01)
    @test isapprox(o.omat[:, 3], [68.78, 100.00, 100.00], atol=0.01)

    # summarize
    o = tab(lbw, :Smoke, summarize=:BWt)
    @test isapprox(o.omat[:, 1], [115, 74, 189])
    @test isapprox(o.omat[:, 2], [3054.96, 2772.30, 2944.29], atol=0.01)
    @test isapprox(o.omat[:, 3], [752.41, 659.81, 729.02], atol=0.01)

    # sort
    o = tab(lbw, :Race, sort=true)
    @test isapprox(o.omat[:, 1], [96, 67, 26, 189])
    @test isapprox(o.omat[:, 2], [50.79, 35.45, 13.76, 100.00], atol=0.01)
    @test isapprox(o.omat[:, 3], [50.79, 86.24, 100.00, 100.00], atol=0.01)

    # sort + summarize
    o = tab(lbw, :Race, sort=true, summarize=:BWt)
    @test isapprox(o.omat[:, 1], [96, 67, 26, 189])
    @test isapprox(o.omat[:, 2], [3103.01, 2804.01, 2719.69, 2944.29], atol=0.01)
    @test isapprox(o.omat[:, 3], [727.87, 721.30, 638.68, 729.02], atol=0.01)

    # skipmissing + sort
    o = tab(bfi, :Education, skipmissing=false, sort=true)
    @test isapprox(o.omat[:, 1], [1249, 418, 394, 292, 224, 223, 2800])
    @test isapprox(o.omat[:, 2], [44.61, 14.93, 14.07, 10.43, 8.00, 7.96, 100.00], atol=0.01)
    @test isapprox(o.omat[:, 3], [44.61, 59.54, 73.61, 84.04, 92.04, 100.00, 100.00], atol=0.01)

    # skipmissing + sort + summarize
    o = tab(bfi, :Education, summarize=:Age, skipmissing=false, sort=true)
    @test isapprox(o.omat[:, 1], [1249, 418, 394, 292, 224, 223, 2800])
    @test isapprox(o.omat[:, 2], [27.23, 35.30, 32.98, 31.51, 25.13, 17.95, 28.78], atol=0.01)
    @test isapprox(o.omat[:, 3], [9.45, 10.98, 10.33, 12.25, 10.40, 8.52, 11.13], atol=0.01)
end

@testset "Two-Way Frequency Tables" begin

    # simple table (no percentages)
    o = tab(lbw, :Race, :Low, pct=nothing)
    @test isapprox(map(x -> x[1], o.omat), [73 23 96; 15 11 26; 42 25 67; 130 59 189])

    # simple table (row percentages only)
    o = tab(lbw, :Race, :Low, pct=:r)
    @test isapprox(map(x -> x[1], o.omat), [73 23 96; 15 11 26; 42 25 67; 130 59 189])
    @test isapprox(map(x -> x[2], o.omat), [76.04 23.96 100.0; 57.69 42.31 100.0; 62.69 37.31 100.0; 68.78 31.22 100.0], atol=0.01)

    # simple table (column percentages only)
    o = tab(lbw, :Race, :Low, pct=:c)
    @test isapprox(map(x -> x[1], o.omat), [73 23 96; 15 11 26; 42 25 67; 130 59 189])
    @test isapprox(map(x -> x[2], o.omat), [56.15 38.98 50.79; 11.54 18.64 13.76; 32.31 42.37 35.45; 100.0 100.0 100.0], atol=0.01)

    # simple table (cell percentages only)
    o = tab(lbw, :Race, :Low, pct=:e)
    @test isapprox(map(x -> x[1], o.omat), [73 23 96; 15 11 26; 42 25 67; 130 59 189])
    @test isapprox(map(x -> x[2], o.omat), [38.62 12.17 50.79; 7.94 5.82 13.76; 22.22 13.23 35.45; 68.78 31.22 100.0], atol=0.01)

    # simple table (pct = :rc)
    o = tab(lbw, :Race, :Low, pct=:rc)
    @test isapprox(map(x -> x[1], o.omat), [73 23 96; 15 11 26; 42 25 67; 130 59 189])
    @test isapprox(map(x -> x[2], o.omat), [76.04 23.96 100.0; 57.69 42.31 100.0; 62.69 37.31 100.0; 68.78 31.22 100.0], atol=0.01)
    @test isapprox(map(x -> x[3], o.omat), [56.15 38.98 50.79; 11.54 18.64 13.76; 32.31 42.37 35.45; 100.0 100.0 100.0], atol=0.01)

    # simple table (pct = :rce)
    o = tab(lbw, :Race, :Low, pct=:rce)
    @test isapprox(map(x -> x[1], o.omat), [73 23 96; 15 11 26; 42 25 67; 130 59 189])
    @test isapprox(map(x -> x[2], o.omat), [76.04 23.96 100.0; 57.69 42.31 100.0; 62.69 37.31 100.0; 68.78 31.22 100.0], atol=0.01)
    @test isapprox(map(x -> x[3], o.omat), [56.15 38.98 50.79; 11.54 18.64 13.76; 32.31 42.37 35.45; 100.0 100.0 100.0], atol=0.01)
    @test isapprox(map(x -> x[4], o.omat), [38.62 12.17 50.79; 7.94 5.82 13.76; 22.22 13.23 35.45; 68.78 31.22 100.0], atol=0.01)

    # skipmissing & pct = :rce
    o = tab(bfi, :Gender, :Education, skipmissing=false, pct=:rce)
    @test isapprox(map(x -> x[1], o.omat), [93 103 356 134 152 81 919; 131 189 893 260 266 142 1881; 224 292 1249 394 418 223 2800])
    @test isapprox(map(x -> x[2], o.omat), [10.12 11.21 38.74 14.58 16.54 8.81 100.0; 6.96 10.05 47.47 13.82 14.14 7.55 100.0; 8.0 10.43 44.61 14.07 14.93 7.96 100.0], atol=0.1)
    @test isapprox(map(x -> x[3], o.omat), [41.52 35.27 28.5 34.01 36.36 36.32 32.82; 58.48 64.73 71.5 65.99 63.64 63.68 67.18; 100.0 100.0 100.0 100.0 100.0 100.0 100.0], atol=0.1)
    @test isapprox(map(x -> x[4], o.omat), [3.32 3.68 12.71 4.79 5.43 2.89 32.82; 4.68 6.75 31.89 9.29 9.5 5.07 67.18; 8.0 10.43 44.61 14.07 14.93 7.96 100.0], atol=0.1)

    # summarize
    o = tab(lbw, :Race, :Low, summarize=:BWt)
    @test isapprox(map(x -> x[1], o.omat), [3407.32 2137.17 3103.01; 3151.8 2130.45 2719.69; 3255.5 2045.52 2804.01; 3328.78 2097.08 2944.29], atol=0.1)
    @test isapprox(map(x -> x[2], o.omat), [522.61 333.75 727.87; 370.52 406.93 638.68; 405.89 439.67 721.3; 478.11 390.88 729.02], atol=0.1)
    @test isapprox(map(x -> x[3], o.omat), [73 23 96; 15 11 26; 42 25 67; 130 59 189])

    # skipmissing + summarize
    o = tab(bfi, :Gender, :Education, summarize=:Age, skipmissing=false)
    @test isapprox(map(x -> x[1], o.omat), [25.15 31.47 25.39 33.17 33.95 18.88 28.02; 25.12 31.54 27.96 32.88 36.08 17.42 29.15; 25.13 31.51 27.23 32.98 35.3 17.95 28.78], atol=0.1)
    @test isapprox(map(x -> x[2], o.omat), [9.82 12.27 8.13 10.67 12.0 9.37 11.03; 10.83 12.27 9.83 10.17 10.29 7.98 11.16; 10.4 12.25 9.45 10.33 10.98 8.52 11.13], atol=0.1)
    @test isapprox(map(x -> x[3], o.omat), [93.0 103.0 356.0 134.0 152.0 81.0 919.0; 131.0 189.0 893.0 260.0 266.0 142.0 1881.0; 224.0 292.0 1249.0 394.0 418.0 223.0 2800.0], atol=0.1)

    # immediate table to conduct chisq test - p-value is for interactive use only
    o = tab([73 23; 15 11; 42 25])
    @test isapprox(map(x -> x[1], o.omat), [73 23 96; 15 11 26; 42 25 67; 130 59 189])
end


@testset "Three-Way Frequency" begin
    # simple table (no percentages)
    o = tab(lbw, :Race, :Low, :Smoke, pct=nothing)
    @test isapprox(map(x -> x[1], o.omat[1]), [40 4 44; 11 5 16; 35 20 55; 86 29 115])
    @test isapprox(map(x -> x[1], o.omat[2]), [33 19 52; 4 6 10; 7 5 12; 44 30 74])

    # simple table (row percentages only)
    o = tab(lbw, :Race, :Low, :Smoke, pct=:r)
    @test isapprox(map(x -> x[2], o.omat[1]), [90.91 9.09 100.0; 68.75 31.25 100.0; 63.64 36.36 100.0; 74.78 25.22 100.0], atol=0.1)
    @test isapprox(map(x -> x[2], o.omat[2]), [63.46 36.54 100.0; 40.0 60.0 100.0; 58.33 41.67 100.0; 59.46 40.54 100.0], atol=0.1)

    # simple table (column percentages only)
    o = tab(lbw, :Race, :Low, :Smoke, pct=:c)
    @test isapprox(map(x -> x[2], o.omat[1]), [46.51 13.79 38.26; 12.79 17.24 13.91; 40.70 68.97 47.83; 100.0 100.0 100.0], atol=0.01)
    @test isapprox(map(x -> x[2], o.omat[2]), [75.00 63.33 70.27; 9.09 20.00 13.51; 15.91 16.67 16.22; 100.0 100.0 100.0], atol=0.01)

    # # simple table (cell percentages only)
    o = tab(lbw, :Race, :Low, :Smoke, pct=:e)
    @test isapprox(map(x -> x[2], o.omat[1]), [34.78 3.48 38.26; 9.57 4.35 13.91; 30.43 17.39 47.83; 74.78 25.22 100.0], atol=0.1)
    @test isapprox(map(x -> x[2], o.omat[2]), [44.59 25.68 70.27; 5.41 8.11 13.51; 9.46 6.76 16.22; 59.46 40.54 100.0], atol=0.1)

    # simple table (pct = :rc)
    o = tab(lbw, :Race, :Low, :Smoke, pct=:rc)
    @test isapprox(map(x -> x[2], o.omat[1]), [90.91 9.09 100.0; 68.75 31.25 100.0; 63.64 36.36 100.0; 74.78 25.22 100.0], atol=0.1)
    @test isapprox(map(x -> x[2], o.omat[2]), [63.46 36.54 100.0; 40.0 60.0 100.0; 58.33 41.67 100.0; 59.46 40.54 100.0], atol=0.1)
    @test isapprox(map(x -> x[3], o.omat[1]), [46.51 13.79 38.26; 12.79 17.24 13.91; 40.70 68.97 47.83; 100.0 100.0 100.0], atol=0.1)
    @test isapprox(map(x -> x[3], o.omat[2]), [75.00 63.33 70.27; 9.09 20.00 13.51; 15.91 16.67 16.22; 100.0 100.0 100.0], atol=0.1)

    # simple table (pct = :rce)
    o = tab(lbw, :Race, :Low, :Smoke, pct=:rce)
    @test isapprox(map(x -> x[2], o.omat[1]), [90.91 9.09 100.0; 68.75 31.25 100.0; 63.64 36.36 100.0; 74.78 25.22 100.0], atol=0.1)
    @test isapprox(map(x -> x[2], o.omat[2]), [63.46 36.54 100.0; 40.0 60.0 100.0; 58.33 41.67 100.0; 59.46 40.54 100.0], atol=0.1)
    @test isapprox(map(x -> x[3], o.omat[1]), [46.51 13.79 38.26; 12.79 17.24 13.91; 40.70 68.97 47.83; 100.0 100.0 100.0], atol=0.1)
    @test isapprox(map(x -> x[3], o.omat[2]), [75.00 63.33 70.27; 9.09 20.00 13.51; 15.91 16.67 16.22; 100.0 100.0 100.0], atol=0.1)
    @test isapprox(map(x -> x[4], o.omat[1]), [34.78 3.48 38.26; 9.57 4.35 13.91; 30.43 17.39 47.83; 74.78 25.22 100.0], atol=0.1)
    @test isapprox(map(x -> x[4], o.omat[2]), [44.59 25.68 70.27; 5.41 8.11 13.51; 9.46 6.76 16.22; 59.46 40.54 100.0], atol=0.1)

    # skipmissing & pct = :rce
    # missing values are artrarily created
    lbw2 = copy(lbw)
    lbw2.Race = convert(Vector{Union{Missing,Int8}}, lbw2.Race)
    lbw2.Low = convert(Vector{Union{Missing,Int8}}, lbw2.Low)
    lbw2.Smoke = convert(Vector{Union{Missing,Int8}}, lbw2.Smoke)
    lbw2[[10, 11, 12], :Race] .= missing
    lbw2[[15, 21, 32], :Low] .= missing
    lbw2[[1, 5, 12], :Smoke] .= missing
    o = tab(lbw2, :Race, :Low, :Smoke, skipmissing=false, pct=:rce)
    # 1 - counts
    @test isapprox(map(x -> x[1], o.omat[1]), [40 4 0 44; 10 5 0 15; 32 20 1 53; 1 0 0 1; 83 29 1 113])
    @test isapprox(map(x -> x[1], o.omat[2]), [29 19 2 50; 4 6 0 10; 7 5 0 12; 1 0 0 1; 41 30 2 73])
    @test isapprox(map(x -> x[1], o.omat[3]), [1 0 0 1; 1 0 0 1; 0 0 0 0; 1 0 0 1; 3 0 0 3])
    # 2 - row percentages
    @test isapprox(map(x -> x[2], o.omat[1]), [90.91 9.09 0.00 100.0; 66.67 33.33 0.00 100.0; 60.38 37.74 1.89 100.0; 100.00 0.00 0.00 100.0; 73.45 25.66 0.88 100.0], atol=0.1)
    @test isapprox(map(x -> x[2], o.omat[2]), [58.00 38.00 4.00 100.0; 40.00 60.00 0.00 100.0; 58.33 41.67 0.00 100.0; 100.0 0.00 0.00 100.0; 56.16 41.10 2.74 100.00], atol=0.1)
    # @test isapprox(map(x -> x[2], o.omat[3]), [33.33 NaN NaN 33.33; 33.33 NaN NaN 33.33; 0.0 NaN NaN 0.0; 33.33 NaN NaN 33.33; 100.0 NaN NaN 100.0], atol=0.1)
    # 3 - column percentages
    @test isapprox(map(x -> x[3], o.omat[1]), [48.19 13.79 0.00 38.94; 12.05 17.24 0.00 13.27; 38.55 68.97 100.00 46.90; 1.20 0.00 0.00 0.88; 100.0 100.0 100.0 100.0], atol=0.1)
    @test isapprox(map(x -> x[3], o.omat[2]), [70.73 63.33 100.0 68.49; 9.76 20.00 0.00 13.70; 17.07 16.67 0.00 16.44; 2.44 0.00 0.00 1.37; 100.0 100.0 100.0 100.0], atol=0.1)
    # @test isapprox(map(x -> x[3], o.omat[3]), [[33.33 NaN NaN 33.33; 33.33 NaN NaN 33.33; 0.0 NaN NaN 0.0; 33.33 NaN NaN 33.33; 100.0 NaN NaN 100.0]], atol=0.1)
    @test isapprox(map(x -> x[4], o.omat[1]), [35.40 3.54 0.00 38.94; 8.85 4.42 0.00 13.27; 28.32 17.70 0.88 46.90; 0.88 0.00 0.00 0.88; 73.45 25.66 0.88 100.00], atol=0.1)
    @test isapprox(map(x -> x[4], o.omat[2]), [39.73 26.03 2.74 68.49; 5.48 8.22 0.00 13.70; 9.59 6.85 0.00 16.44; 1.37 0.00 0.00 1.37; 56.16 41.10 2.74 100.00], atol=0.1)

    # summarize
    o = tab(lbw, :Race, :Low, :Smoke, summarize=:BWt)
    # 1 - mean
    @test isapprox(map(x -> x[1], o.omat[1]), [3582.73 1889.0 3428.75; 3240 2149.8 2876.6; 3292.5313 2056.6 2817.17; 3427.68 2049.55 3067.63], atol=0.1)
    @test isapprox(map(x -> x[1], o.omat[2]), [], atol=0.1)
    @test isapprox(map(x -> x[1], o.omat[3]), [], atol=0.1)
    @test isapprox(map(x -> x[2], o.omat), [522.61 333.75 727.87; 370.52 406.93 638.68; 405.89 439.67 721.3; 478.11 390.88 729.02], atol=0.1)
    @test isapprox(map(x -> x[3], o.omat), [73 23 96; 15 11 26; 42 25 67; 130 59 189])

    # # skipmissing + summarize
    # o = tab(bfi, :Gender, :Education, summarize=:Age, skipmissing=false)
    # @test isapprox(map(x -> x[1], o.omat), [25.15 31.47 25.39 33.17 33.95 18.88 28.02; 25.12 31.54 27.96 32.88 36.08 17.42 29.15; 25.13 31.51 27.23 32.98 35.3 17.95 28.78], atol=0.1)
    # @test isapprox(map(x -> x[2], o.omat), [9.82 12.27 8.13 10.67 12.0 9.37 11.03; 10.83 12.27 9.83 10.17 10.29 7.98 11.16; 10.4 12.25 9.45 10.33 10.98 8.52 11.13], atol=0.1)
    # @test isapprox(map(x -> x[3], o.omat), [93.0 103.0 356.0 134.0 152.0 81.0 919.0; 131.0 189.0 893.0 260.0 266.0 142.0 1881.0; 224.0 292.0 1249.0 394.0 418.0 223.0 2800.0], atol=0.1)



end
