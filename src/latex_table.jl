function name_classified_row(row)
    (spike_motif=row[1],  n1=row[2], n2=row[3], t1=row[4], t2=row[5], connection_motif=row[6])
end

sign_str(x) = if x > 0
    "+"
elseif x < 0
    "-"
else
    "0"
end
abs_rel_sym(x1, x2) = if abs(x1) < abs(x2)
    "<"
elseif abs(x1) > abs(x2)
    ">"
else
    "="
end
function prettify_row(row)
    conditions = filter(x -> x != "", [
        if sign_str(row.n1) == sign_str(row.n2) && row.n1 != 0
            "\$\\abs{n_1} " * abs_rel_sym(row.n1, row.n2) * " \\abs{n_2}\$"
        else
            ""
        end, if sign_str(row.t1) == sign_str(row.t2) && row.t1 != 0
            "\$\\abs{t_1} " * abs_rel_sym(row.t1, row.t2) * " \\abs{t_2}\$"
        else
            ""
        end])
    conditions_str = if length(conditions) > 1
        "\\adjustbox{raise=-0.25\\height}{\\shortstack{" * join(conditions, "\\\\") * "}}"
    elseif length(conditions) == 0
        ""
    else
        only(conditions)
    end
    (
        spike_motif = row.spike_motif,
        n1 = sign_str(row.n1),
        n2 = sign_str(row.n2),
        t1 = sign_str(row.t1),
        t2 = sign_str(row.t2),
        connection_motif = roman_encode(row.connection_motif),
        conditions = conditions_str
    )
end

using LaTeXStrings

function latexify(df::DataFrame)
    """
Lag Motif & \$n_1\$ & \$n_2\$ & \$t_1\$ & \$t_2\$ & Constraints & Information Flow Pattern & Configuration \\\\
    \\hline
    """ *
    join(latexify.(eachrow(df)), "\\\\\n")
end

function graphics_scaling_string(n_same_neurons, n_same_times)
    if n_same_neurons == 3 && n_same_times == 3
        return "height=1.5cm"
    else
        return "scale=1.0"
    end
end

function latexify(row::DataFrameRow)

    row_contents = join(getproperty.(Ref(row), ["spike_motif", "n1", "n2", "t1", "t2", "conditions", "connection_motif"]), "&")
    row_contents *= "& \\adjustbox{padding=0em 0.5em 0em 0.5em, raise=-0.33\\height}{\\includegraphics{motif_figs/$(row.spike_motif).png}}"
    return row_contents
end