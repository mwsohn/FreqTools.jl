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
function tab(na::NamedArray; skipmissing=true, pct=:rce, digits=2)

    len = length(na.dimnames)
    if len == 1
        _tab1(na; skipmissing=skipmissing, digits=digits)
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
        return _tab1summarize(Tables.getcolumn(indf, var), Tables.getcolumn(indf, summarize), skipmissing=skipmissing, digits = digits, varname=string(var))
    end
    _tab1(freqtable(indf, var, skipmissing=skipmissing); sort=sort, digits=digits)
end
function tab(ivar::AbstractVector; skipmissing=true, sort=false)
    _tab1(freqtable(ivar, skipmissing=skipmissing); sort=sort)
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

    if summarize == nothing
        na = freqtable(indf, var1, var2, var3, skipmissing=skipmissing)
        # stratify the var3 (na.dimnames[3])
        n3 = size(na, 3)
        vals = na.dicts[3].keys

        for i in 1:n3
            println("\n\n", na.dimnames[3], " = ", vals[i], "\n")

            _tab2(na[:, :, i]; pct=pct, maxrows=maxrows, maxcols=maxcols)
        end
    else
        n3 = sort(unique(indf[!, var3]))
        for v in n3
            println("\n\n", var3, " = ", v, "\n")

            subdf = filter(x -> x[var3] == v, indf)
            _tab2summarize(Tables.getcolumn(subdf,var1), 
                Tables.getcolumn(subdf,var2), 
                Tables.getcolumn(subdf,summarize); maxrows=-1, maxcols=20, 
                digits=digits, varnames=string(var1, " ╲ ", var2))
        end
    end
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
        arry = na[s]
    else
        arry = na.array
    end

    # value labels and "Total"
    if sort
        rownames = vcat(string.(names(na)[1][s]), "Total")
    else
        rownames = vcat(string.(names(na)[1]), "Total")
    end

    # counts - the last row has the total
    counts = vcat(arry, sum(na, dims=1))

    # percents
    percents = 100 .* counts ./ counts[end]

    # cumulative percents
    cumpct = 100 .* vcat(cumsum(arry, dims=1), counts[end]) ./ counts[end]

    omat = hcat(counts, percents, cumpct)

    return TAB1OUT(omat, rownames, string(dimnames(na)[1]), digits)

    # # do not output rows with zeros
    # z = findall(x -> x != 0, na.array)
    # arry = na.array[z]

    # if sort
    #     s = sortperm(arry, rev=true)
    #     arry = arry[s]
    # end

    # # value labels and "Total"
    # if sort
    #     rownames = vcat(names(na)[1][z][s], "Total")
    # else
    #     rownames = vcat(names(na)[1][z], "Total")
    # end

    # # counts - the last row has the total
    # counts = vcat(arry, sum(na, dims=1))

    # # percents
    # percents = 100 .* counts ./ counts[end]

    # # cumulative percents
    # cumpct = 100 .* vcat(cumsum(arry, dims=1), counts[end]) ./ counts[end]

    # ar = hcat(counts, percents, cumpct)

    # fmt = Printf.Format("%.$(digits)f")
    # PrettyTables.pretty_table(ar,
    #     row_labels = rownames,
    #     row_label_column_title=string(na.dimnames[1]),
    #     header=["Counts", " Percent", "Cum Pct"],
    #     formatters= (v,i,j) -> j in (2,3) ? Printf.format(fmt,v) : @sprintf("%.0d", v), 
    #     crop=:none,
    #     hlines=[0, 1, length(rownames), length(rownames) + 1],
    #     vlines=[1])
end
struct TAB1OUT
    omat::Matrix
    rownames::Vector{String}
    varname::String
    digits::Int8
end

function Base.show(io::IO, m::TAB1OUT)
    fmt = Printf.Format("%.$(m.digits)f")
    PrettyTables.pretty_table(io, m.omat,
        row_labels=m.rownames,
        row_label_column_title=string(m.varname),
        header=["Counts", " Percent", "Cum Pct"],
        formatters=(v, i, j) -> j in (2, 3) ? Printf.format(fmt, v) : @sprintf("%.0d", v),
        crop=:none,
        hlines=[0, 1, length(m.rownames), length(m.rownames) + 1],
        vlines=[1])
end

function _tab1summarize(var, sumvar; skipmissing=false, digits=2, order = false, varname=nothing)

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
    if order == true
        p = sortperm(omat[:, 1], rev=true)
        push!(p, popfirst!(p))
        omat = omat[p, :]
        rownames = rownames[p]
    end

    return TAB1OUT2(omat, string.(rownames), ["N", "Mean", "SD"], varname == nothing ? "" : varname, digits)
    # df = DataFrame(t=var, s=sumvar)
    # if skipmissing
    #     ba = completecases(df)
    # else
    #     ba = completecases(df[:, [:s]])
    # end
    # df = df[ba, :]

    # by = sort(unique(df[:, :t]))
    # len = size(by, 1)

    # omat = Matrix{Float64}(undef, len + 1, 3)
    # @inbounds for subdf in groupby(df, :t, sort=true)
    #     i = findfirst(x -> x == subdf[1, :t], by)
    #     if size(subdf, 1) == 0 || subdf == nothing
    #         omat[i, 1:3] .= (0, NaN, NaN)
    #     else
    #         omat[i, 1:3] .= (size(subdf, 1), mean(subdf[:, :s]), std(subdf[:, :s]))
    #     end
    # end
    # omat[len+1, 1:3] .= (size(df, 1), mean(df[:, :s]), std(df[:, :s]))

    # fmt = Printf.Format("%.$(digits)f")
    # pretty_table(omat,
    #     row_labels=string.(vcat(by, "Total")),
    #     row_label_column_title=varname == nothing ? "" : varname,
    #     formatters=(v, _, j) -> isnan(v) ? "." : (j % 3 != 1 ? Printf.format(fmt, v) : @sprintf("%.0f", v)),
    #     header=["N", "Mean", "SD"],
    #     crop=:none,
    #     hlines=vcat([0, 1], len + 1, len + 2),
    #     vlines=[1])
end
struct TAB1OUT2
    omat::Matrix
    rownames::Vector{String}
    colnames::Vector{String}
    varname::String
    digits::Int8
end
function Base.show(io::IO, m::TAB1OUT2)
    fmt = Printf.Format("%.$(m.digits)f")
    len = size(m.omat, 1)

    pretty_table(io,
        m.omat,
        row_labels=m.rownames,
        row_label_column_title=m.varname,
        formatters=(v, _, j) -> isnan(v) ? "." : (j % 3 != 1 ? Printf.format(fmt, v) : string(v)),
        header=m.colnames,
        crop=:none,
        hlines=vcat([0, 1], len, len + 1),
        vlines=[1])
end

function _tab2(na::NamedArray; maxrows=-1, maxcols=20, pct=:rce, digits = 2, tests=true)

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

struct TAB2OUT
    omat::Matrix
    rownames::Vector{String}
    colnames::Vector{String}
    varnames::String
    maxrows::Int64
    maxcols::Int64
    digits::Int8
    tests::Bool
end

function Base.show(io::IO, m::TAB2OUT)

    pretty_table(io,
        _tab2matstr(m),
        linebreaks=true,
        row_labels=m.rownames,
        row_label_column_title=m.varnames,
        header=m.colnames,
        crop=:none,
        max_num_of_rows=m.maxrows,
        max_num_of_columns=m.maxcols,
        hlines=vcat([0, 1], [x + 1 for x in 1:length(m.rownames)]),
        vlines=[1])

    if m.tests == true
        (i, j) = size(m.omat)
        testarray = map(x -> x[1], m.omat)[1:i-1, 1:j-1]
        if size(testarray, 1) > 1 && size(testarray, 2) > 1
            c = ChisqTest(testarray)
            pval = pvalue(c)
            println(io, "Pearson chi-square = ", @sprintf("%.4f", c.stat), " (", c.df, "), p ",
                pval < 0.0001 ? "< 0.0001" : string("= ", round(pval, sigdigits=6)))
        end

        if size(testarray) == (2, 2) && all(x -> x > 0, testarray) # 2x2 array
            println(io, "Fisher's exact test = ", @sprintf("%.4f",
                pvalue(HypothesisTests.FisherExactTest((testarray')...))))
        end
    end
end
function _tab2matstr(m)
    # covert omat to string values
    fmt = Printf.Format("%.$(m.digits)f")
    o = copy(m.omat)
    (nrows,ncols) = size(o)

    # non-missing cell
    e = o[findfirst(x -> !ismissing(x), o)]
    etypeint = typeof(e[1]) <: Integer
    elen = length(e)
    for i in 1:nrows
        for j in 1:ncols
            t = o[i, j]
            if ismissing(t)
                if etypeint
                    o[i, j] = string("0",repeat("\n.",elen-1))
                else
                    o[i, j] = string(repeat(".\n",elen-1),"0")
                end
            else
                tt = [ ismissing(x) ? "." : (typeof(x) <: Integer ? string(x) : Printf.format(fmt,x)) for x in t]
                o[i,j] = join(tt,"\n")
            end
        end
    end

    return o
end