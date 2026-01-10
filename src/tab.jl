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
        varnames = string(var1," \\ ", var2),
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
                digits=digits, varnames = string(var1," \\ ",var2))
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

    # do not output rows with zeros
    z = findall(x -> x != 0, na.array)
    arry = na.array[z]

    if sort
        s = sortperm(arry, rev=true)
        arry = arry[s]
    end

    # value labels and "Total"
    if sort
        rownames = vcat(names(na)[1][z][s], "Total")
    else
        rownames = vcat(names(na)[1][z], "Total")
    end

    # counts - the last row has the total
    counts = vcat(arry, sum(na, dims=1))

    # percents
    percents = 100 .* counts ./ counts[end]

    # cumulative percents
    cumpct = 100 .* vcat(cumsum(arry, dims=1), counts[end]) ./ counts[end]

    ar = hcat(counts, percents, cumpct)

    fmt = Printf.Format("%.$(digits)f")

    PrettyTables.pretty_table(ar,
        row_labels = rownames,
        row_label_column_title=string(na.dimnames[1]),
        header=["Counts", " Percent", "Cum Pct"],
        formatters= (v,i,j) -> j in (2,3) ? Printf.format(fmt,v) : @sprintf("%.0d", v), 
        crop=:none,
        hlines=[0, 1, length(rownames), length(rownames) + 1],
        vlines=[1])
end

function _tab2(na::NamedArray; maxrows=-1, maxcols=20, pct=:rce, digits = 2)

    # counts
    counts = na.array
    counts = vcat(counts, sum(counts, dims=1)) # column sum
    counts = hcat(counts, sum(counts, dims=2)) # row sum

    # drop rows or columns whose margin totals == 0
    rz = findall(x -> x != 0, counts[:, end]) # find all columns with non-zero totals
    cz = findall(x -> x != 0, counts[end, :]) # find all rows with non-zero totals
    counts = counts[rz, cz]
    (nrow, ncol) = size(counts)

    # row and column labels and "Total"
    rownames = vcat(names(na)[1], "Total")[rz]
    colnames = vcat(names(na)[2], "Total")[cz]

    # compute percentages
    if pct == nothing
        pctstr = ""
    else
        pctstr = string(pct)
    end
    cnt = length(pctstr) + 1
    d = Matrix{Any}(undef, nrow * cnt, ncol)
    for i in 1:nrow
        r = i + (i - 1) * (cnt - 1)
        for j in 1:ncol

            # counts
            d[r, j] = counts[i, j]

            # percentages
            for (k, v) in enumerate(pctstr)
                if v == 'r' # row percents
                    d[r+k, j] = 100 * counts[i, j] / counts[i, ncol]
                elseif v == 'c' # column percents
                    d[r+k, j] = 100 * counts[i, j] / counts[nrow, j]
                elseif v == 'e' # cell percents
                    d[r+k, j] = 100 * counts[i, j] / counts[nrow, ncol]
                end
            end
        end
    end

    # add blank cells
    rownames2 = vcat([vcat(x, fill(" ", cnt - 1)) for x in rownames]...)

    fmt = Printf.Format("%.$(digits)f")

    pretty_table(d,
        row_labels=rownames2,
        row_label_column_title=string(na.dimnames[1], " \\ ", na.dimnames[2]),
        header=colnames,
        crop=:none,
        formatters=(v, i, _) -> (cnt == 1 || i % cnt == 1) ? @sprintf("%.0f", v) : Printf.format(fmt,v),
        max_num_of_rows=maxrows,
        max_num_of_columns=maxcols,
        hlines=vcat([0, 1], cnt == 1 ? [nrow,nrow+1] : [x * cnt + 1 for x in 1:(nrow+1)]),
        vlines=[1])

    testarray = na.array[rz[1:end-1], cz[1:end-1]]
    if size(testarray, 1) > 1 && size(testarray, 2) > 1
        c = ChisqTest(testarray)
        pval = pvalue(c)
        println("Pearson chi-square = ", @sprintf("%.4f", c.stat), " (", c.df, "), p ",
            pval < 0.0001 ? "< 0.0001" : string("= ", round(pval, sigdigits=6)))
    end

    if size(testarray) == (2, 2) && all(x -> x > 0, testarray) # 2x2 array
        println("Fisher's exact test = ", @sprintf("%.4f",
            pvalue(HypothesisTests.FisherExactTest((testarray')...))))
    end
end

function _tab1summarize(var, sumvar; skipmissing=false, digits = 2, varname = nothing)

    if skipmissing == false
        m = .!(ismissing.(var) .| ismissing.(sumvar))
    else
        m = .!ismissing.(sumvar)
    end
    var2 = var[m]
    sumvar2 = sumvar[m]

    by = sort(unique(var2))
    len = size(by, 1)

    omat = Matrix{Float64}(undef, len + 1, 3)
    @inbounds for (i, v) in enumerate(by)
        if ismissing(v)
            inc = [ismissing(x) ? true : false for x in var2]
        else
            inc = [x == v ? true : false for x in var2]
        end
        omat[i, 1:3] .= (sum(inc), mean(sumvar2[inc]), std(sumvar2[inc]))
    end
    omat[len+1, 1:3] .= (length(var2), mean(sumvar2), std(sumvar2))

    fmt = Printf.Format("%.$(digits)f")
    pretty_table(omat,
        row_labels=string.(vcat(by, "Total")),
        row_label_column_title=varname == nothing ? "" : varname,
        formatters=(v, _, j) -> isnan(v) ? "." : (j % 3 != 1 ? Printf.format(fmt, v) : @sprintf("%.0f", v)),
        header=["N", "Mean", "StDev"],
        crop=:none,
        hlines=vcat([0, 1], len + 1, len + 2),
        vlines=[1])
end

function _tab2summarize(var1, var2, sumvar; maxrows=-1, maxcols=20, skipmissing=nothing, varnames = nothing, digits=2)

    mm = hcat(var1, var2, sumvar)

    if skipmissing == true
        m = vec(sum(ismissing.(mm), dims=2) .== false)
    else
        m = vec(ismissing.(mm[:, 3]) .== false)
    end
    mm2 = mm[m, :]

    vv1 = sort(unique(var1[m]))
    vv2 = sort(unique(var2[m]))
    nrows = length(vv1)
    ncols = length(vv2)

    # allocate memory for output matrix
    omat = Matrix{Union{Missing,Any}}(missing, nrows + 1, ncols + 1)
    for (i, v1) in enumerate(vv1)
        # row margin
        if ismissing(v1)
            subm1 = @view mm2[ismissing.(mm2[:, 1]), :]
        else
            subm1 = @view mm2[[!ismissing(x) && x == v1 for x in mm2[:, 1]], :]
        end
        omat[i, ncols+1] = tuple(mean(subm1[:,3]), std(subm1[:, 3]), size(subm1, 1))

        for (j, v2) in enumerate(vv2)
            if ismissing(v1) && !ismissing(v2)
                subm2 = @view subm1[[ismissing(x[1]) && x[2] == v2 for x in eachrow(subm1[:, 1:2])], :]
            elseif !ismissing(v1) && ismissing(v2)
                subm2 = @view subm1[[ismissing(x[2]) && x[1] == v1 for x in eachrow(subm1[:, 1:2])], :]
            else
                subm2 = @view subm1[[!ismissing(x[1]) && x == (v1, v2) for x in tuple.(subm1[:, 1], subm1[:, 2])], :]
            end
            omat[i, j] = tuple(mean(subm2[:, 3]), std(subm2[:, 3]), size(subm2, 1))
            if i == 1
                if ismissing(v2)
                    subm2 = @view mm2[ismissing.(mm2[:, 2]), :]
                else
                    subm2 = @view mm2[[!ismissing(x) && x == v2 for x in mm2[:, 2]], :]
                end
                omat[nrows+1, j] = tuple(mean(subm2[:, 3]), std(subm2[:, 3]), size(subm2, 1))
            end
        end
    end
    omat[nrows+1, ncols+1] = tuple(mean(mm2[:, 3]), std(mm2[:, 3]), size(mm2, 1))

    # value labels and "Total"
    rownames = string.(vcat(vv1, "Total"))
    # rownames = vcat([[x, " ", " "] for x in rownames]...)

    # colunm names
    colnames = string.(vcat(vv2, "Total"))

    # covert omat to string values
    fmt = Printf.Format("%.$(digits)f")
    for i in 1:nrows+1
        for j in 1:ncols+1
            t = omat[i, j]
            omat[i, j] = string(isnan(t[1]) ? "." : Printf.format(fmt, t[1]), "\n",
                isnan(t[2]) ? "." : Printf.format(fmt, t[2]), "\n",
                t[3])
        end
    end

    # output
    pretty_table(omat,
        linebreaks=true,
        row_labels=rownames,
        row_label_column_title=varnames == nothing ? "" : varnames,
        header=colnames,
        crop=:none,
        max_num_of_rows=maxrows,
        max_num_of_columns=maxcols,
        hlines=vcat([0, 1], [x + 1 for x in 1:(nrows+1)]),
        vlines=[1])
end

