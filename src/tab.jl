"""
    tab(df::DataFrame,vars::Symbol...; sort = true, skipmissing=true)
    tab(na::NamedArray)
    tab(m::Matrix)

Produces an one-way or two-way frequency table from a DataFrame, other Tables objects, or a NamedArray returned from
`freqtable` function. `tab` is mainly a wrapper for the excellent `FreqTables` package.

Use `skipmissing = false` to obtain frequencies that include `missing` values.

For an one-way frequency table, the table can be
ordered by the frequency by specifying `sort = true` as an option.

For an one-way or two-way table, summary values (mean, standard deviation, and count) 
of a continuous variable can be requested by specifying the variable name
as an option `summarize = :var`.

For a two-way table, percentages can be requested by specifying
a `pct = :rce` option in any combination 
where `r` indicates "row" percents, `c` "column" percents, and `e` "cell" percents
(default is `:rce`). They can be use in any combination or any order.

If input is a matrix of counts, a Pearson chi-square test will be performed. 
"""
function tab(na::NamedArray; skipmissing=true, pct=:rce, sort = false, digits=2)

    len = length(na.dimnames)
    if len == 1
        _tab1(na; skipmissing=skipmissing, sort = sort, digits=digits)
    elseif len == 2
        _tab2(na; skipmissing=skipmissing, pct=pct)
    elseif len == 3
        _tab3(na; skipmissing=skipmissing, pct=pct)
    end
    throw(ArgumentError("Crosstabs for more than 3 variables are not currently supported."))
end
function tab(indf, var::Union{Symbol,String}; skipmissing=true, sort=false, summarize=nothing, digits = 2)
    s = Tables.schema(indf)
    if in(Symbol(var), s.names) == false
        throw(ArgumentError("$var is not found in the input table."))
        return nothing
    end
    if summarize != nothing
        return _tab1summarize(Tables.getcolumn(indf, var), Tables.getcolumn(indf, summarize), skipmissing=skipmissing, digits = digits, sort = sort, varname=string(var))
    end
    _tab1(freqtable(indf, var, skipmissing=skipmissing), sort=sort, digits=digits)
end
function tab(ivar::AbstractVector; skipmissing=true, sort=false, digits = 2)
    _tab1(freqtable(ivar, skipmissing=skipmissing); sort=sort, digits = digits)
end
function tab(indf,var1::Union{Symbol,String}, var2::Union{Symbol,String}; pct=:rce, maxrows=-1, maxcols=20, skipmissing=true, summarize=nothing, digits=2)

    s = Tables.schema(indf)
    if in(Symbol(var1), s.names) == false
        throw(ArgumentError("$var1 is not found in the input DataFrame."))
        return nothing
    end
    if in(Symbol(var2), s.names) == false
        throw(ArgumentError("$var2 is not found in the input DataFrame."))
        return nothing
    end
    if summarize == nothing
        return _tab2(freqtable(indf, var1, var2, skipmissing=skipmissing); pct=pct, maxrows=maxrows, maxcols=maxcols)
    end

    _tab2summarize(Tables.getcolumn(indf,var1), 
        Tables.getcolumn(indf,var2), 
        Tables.getcolumn(indf,summarize); 
        maxrows=-1, 
        maxcols=20,
        varnames=string(var1, " ╲ ", var2),
        digits=digits)
end
function tab(indf, var1::Union{Symbol,String}, var2::Union{Symbol,String}, var3::Union{Symbol,String};
    pct=:rce, maxrows=-1, maxcols=20, skipmissing=true, summarize=nothing, digits=2)

    s = Tables.schema(indf)
    for v in (var1, var2, var3)
        if in(Symbol(v), s.names) == false
            throw(ArgumentError("$v is not found in the input DataFrame."))
            return nothing
        end
    end

    m = []
    vomat = []
    if summarize == nothing
        na = freqtable(indf, var1, var2, var3, skipmissing=skipmissing)

        # stratify the var3 (na.dimnames[3])
        n3 = size(na, 3)
        vals = na.dicts[3].keys
        thirdnm = string(na.dimnames[3])

        for i in 1:n3
            m = FreqTools._tab2(na[:, :, i], pct=pct, digits=digits, tests=true, maxrows=maxrows, maxcols=maxcols)
            push!(vomat, m.omat)
        end
        return TAB3OUT(vomat, m.rownames, m.colnames, m.varnames, thirdnm, vals, maxrows, maxcols, digits, tests)
    end

    if skipmissing
        vals = sort(unique(collect(skipmissing(Tables.getcolumn(indf, var3)))))
    else
        vals = sort(unique(Tables.getcolumn(indf, var3)))
    end
    n3 = length(vals)
    thirdnm = string(var3)

    for v in vals
        subdf = filter(x -> x[var3] == v, indf)
        m = _tab2summarize(Tables.getcolumn(subdf, var1),
            Tables.getcolumn(subdf, var2),
            Tables.getcolumn(subdf, summarize); maxrows=maxrows, maxcols=maxcols,
            digits=digits, varnames=string(var1, " ╲ ", var2))
        push!(vomat, m.omat)
    end
    return TAB3OUT(vomat, m.rownames, m.colnames, m.varnames, thirdnm, vals, maxrows, maxcols, digits, tests)
end
function tab(a::AbstractArray; pct=:rce, digits=2)
    if length(size(a)) == 2 && all(x -> x >= 2, a)
        _tab2(NamedArray(a), pct=pct, digits=digits)
    else
        throw(ArgumentError("Input array must be 2x2 and have at least two levels on each dimension."))
    end
end

function _tab1(na::NamedArray; sort=false, digits=2)

    if sort
        s = sortperm(na, rev=true)
        arry = na.array[s]
        rownames = vcat(string.(names(na)[1][s]), "Total")
    else
        arry = na.array
        rownames = vcat(string.(names(na)[1]), "Total")
    end

    # counts - the last row has the total
    counts = vcat(arry, sum(na, dims=1))

    # percents
    percents = 100 .* counts ./ counts[end]

    # cumulative percents
    cumpct = 100 .* vcat(cumsum(arry, dims=1), counts[end]) ./ counts[end]

    # combine them
    omat = hcat(counts, percents, cumpct)

    return TAB1OUT(omat, rownames, string(dimnames(na)[1]), digits)
end

function _tab1summarize(var, sumvar; skipmissing=false, digits=2, sort = false, varname=nothing)

    df = DataFrame(t=var, s=sumvar)
    if skipmissing
        df = dropmissing(df)
    else
        df = df[completecases(df[:, [:s]]), :]
    end

    len = length(unique(df[:, :t]))
    omat = Matrix{Any}(missing, len + 1, 3)
    rownames = []
    @inbounds for (i, subdf) in enumerate(groupby(df, :t, sort=true))
        push!(rownames, subdf[1, :t])
        if size(subdf, 1) == 0 || subdf == nothing
            omat[i, 1:3] .= (0, NaN, NaN)
        else
            omat[i, 1:3] .= (size(subdf, 1), mean(subdf[:, :s]), std(subdf[:, :s]))
        end
    end
    omat[len+1, 1:3] .= (size(df, 1), mean(df[:, :s]), std(df[:, :s]))

    push!(rownames, "Total")
    if sort == true
        p = sortperm(omat[:, 1], rev=true)
        push!(p, popfirst!(p))
        omat = omat[p, :]
        rownames = rownames[p]
    end

    return TAB1OUT2(omat, string.(rownames), ["N", "Mean", "SD"], varname == nothing ? "" : varname, digits)
end

function _tab2(na::NamedArray; pct=:rce, digits = 2, tests=true, maxrows=-1, maxcols=20)

    # counts
    counts = na.array
    counts = vcat(counts, sum(counts, dims=1)) # column sum
    counts = hcat(counts, sum(counts, dims=2)) # row sum
    (nrow, ncol) = size(counts)

    # row and column labels and "Total"
    rownames = vcat(string.(names(na)[1]), "Total")
    colnames = vcat(string.(names(na)[2]), "Total")

    # compute percentages
    if pct == nothing
        pctstr = ""
    else
        pctstr = string(pct)
    end
    d = Matrix{Any}(missing, nrow, ncol)
    for i in 1:nrow
        for j in 1:ncol
            vv = []
            push!(vv, counts[i, j])
            for c in pctstr
                if c == 'r' # row percents
                    push!(vv, 100 * counts[i, j] / counts[i, ncol])
                elseif c == 'c' # column percents
                    push!(vv, 100 * counts[i, j] / counts[nrow, j])
                elseif c == 'e' # cell percents
                    push!(vv, 100 * counts[i, j] / counts[nrow, ncol])
                end
            end
            d[i, j] = vv
        end
    end

    return TAB2OUT(d, rownames, colnames, string(join(dimnames(na), " ╲ ")), maxrows, maxcols, digits, tests)
end

function idxdict(vv1,vv2)
    dd = Dict()
    for (i, v1) in enumerate(vv1)
        for (j, v2) in enumerate(vv2)
            dd[(v1, v2)] = (i, j)
        end
    end
    return dd
end

function _tab2summarize(var1, var2, sumvar; maxrows=-1, maxcols=20, skipmissing=nothing, varnames = nothing, digits=2)

    if skipmissing == true
        df = dropmissing(DataFrame(t1=var1, t2=var2, tsumvar=sumvar))
    else
        df = DataFrame(t1=var1, t2=var2, tsumvar=sumvar)[.!ismissing.(sumvar),:]
    end

    vv1 = sort(unique(df[:, :t1]))
    vv2 = sort(unique(df[:, :t2]))
    idx = idxdict(vv1,vv2)
    nrows = length(vv1)
    ncols = length(vv2)

    # allocate memory for output matrix
    omat = Matrix{Union{Missing,Any}}(missing, nrows + 1, ncols + 1)

    # let's implement this using groupby
    for subdf in groupby(df, [:t1, :t2], sort=true)
        (i,j) = idx[(subdf[1,:t1], subdf[1,:t2])]
        omat[i, j] = (mean(subdf[:, :tsumvar]), std(subdf[:, :tsumvar]), size(subdf, 1))
    end
    for (i,subdf) in enumerate(groupby(df, :t1, sort=true))
        omat[i, ncols+1] = (mean(subdf[:, :tsumvar]), std(subdf[:, :tsumvar]), size(subdf, 1))
    end
    for (j,subdf) in enumerate(groupby(df, :t2, sort=true))
        omat[nrows+1, j] = (mean(subdf[:, :tsumvar]), std(subdf[:, :tsumvar]), size(subdf, 1))
    end
    omat[nrows+1, ncols+1] = (mean(df[:, :tsumvar]), std(df[:, :tsumvar]), size(df, 1))

    # value labels and "Total"
    rownames = string.(vcat(vv1, "Total"))

    # colunm names
    colnames = string.(vcat(vv2, "Total"))

    return TAB2OUT(omat, rownames, colnames, varnames == nothing ? "" : varnames, maxrows, maxcols, digits, false)
end
