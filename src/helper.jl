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
        _tab2matstr(m.omat, m.digits),
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
function _tab2matstr(m, digits)

    # covert omat to string values
    fmt = Printf.Format("%.$(digits)f")
    o = copy(m)
    (nrows, ncols) = size(o)

    # non-missing cell
    e = o[findfirst(x -> !ismissing(x), o)]
    etypeint = typeof(e[1]) <: Integer
    elen = length(e)
    for i in 1:nrows
        for j in 1:ncols
            t = o[i, j]
            if ismissing(t)
                if etypeint
                    o[i, j] = string("0", repeat("\n.", elen - 1))
                else
                    o[i, j] = string(repeat(".\n", elen - 1), "0")
                end
            else
                tt = [ismissing(x) || isnan(x) ? "." : (typeof(x) <: Integer ? string(x) : Printf.format(fmt, x)) for x in t]
                o[i, j] = join(tt, "\n")
            end
        end
    end

    return o
end

struct TAB3OUT
    omat::Vector
    rownames::Vector{String}
    colnames::Vector{String}
    varnames::String
    thirdname::String
    thirdval::Vector
    maxrows::Int64
    maxcols::Int64
    digits::Int8
    tests::Bool
end

function Base.show(io::IO, m::TAB3OUT)

    fmt = Printf.Format("%.$(m.digits)f")
    for i in 1:length(m.omat)
        println(io, "\n\n", m.thirdname, " = ", m.thirdval[i], "\n")
        pretty_table(io,
            _tab2matstr(m.omat[i], m.digits),
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
            (i, j) = size(m.omat[i])
            testarray = map(x -> x[1], m.omat[i])[1:i-1, 1:j-1]
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
end