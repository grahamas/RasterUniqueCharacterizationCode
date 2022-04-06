
function _to_plusminus(tup, i_nz)
    if tup[i_nz] == false
        [(tup[begin:i_nz-1]..., 0, tup[i_nz+1:end]...)]
    else
        [(tup[begin:i_nz-1]..., 1, tup[i_nz+1:end]...),
         (tup[begin:i_nz-1]..., -1, tup[i_nz+1:end]...)]
    end
end

function expand_plusminus(nz_tup)
    # Expand nonzero indicator stuple into plus-minuses
    l_pm_tups = [nz_tup]
    for i_nz = 2:5
        new_tups = [_to_plusminus(tup, i_nz) for tup in l_pm_tups]
        l_pm_tups = reduce(vcat, new_tups)
    end
    return if length(l_pm_tups) > 1
         [("$(tup[begin]).$(i-1)", tup[begin+1:end]...) for (i, tup) in enumerate(l_pm_tups)]
    else
        [("$(tup[begin])", tup[begin+1:end]...) for (i, tup) in enumerate(l_pm_tups)]
    end
end

function _to_ordered(tup)
    ns_ordered_tups = if tup[2] == tup[4] == 1
        [
            ("$(tup[begin])", 1, tup[3], 1, tup[5]),
            ("$(tup[begin])", 2, tup[3], 1, tup[5]),
            ("$(tup[begin])", 1, tup[3], 2, tup[5]),
        ]
    elseif tup[2] == tup[4] == -1
        [
            ("$(tup[begin])", -1, tup[3], -1, tup[5]),
            ("$(tup[begin])", -2, tup[3], -1, tup[5]),
            ("$(tup[begin])", -1, tup[3], -2, tup[5]),
        ]
    else
        [tup]
    end
    ns_ordered_tups = map(ns_ordered_tups) do tup
        if tup[3] == tup[5] == 1
            [
                ("$(tup[begin])", tup[2], 1, tup[4], 1),
                ("$(tup[begin])", tup[2], 2, tup[4], 1),
                ("$(tup[begin])", tup[2], 1, tup[4], 2)
            ]
        elseif tup[3] == tup[5] == -1
            [
                ("$(tup[begin])", tup[2], -1, tup[4], -1),
                ("$(tup[begin])", tup[2], -2, tup[4], -1),
                ("$(tup[begin])", tup[2], -1, tup[4], -2)
            ]
        else
            tup
        end
    end 
    if !(eltype(ns_ordered_tups) <: Tuple)
        ns_ordered_tups = reduce(vcat, ns_ordered_tups)
    end

    if length(ns_ordered_tups) != 1
        [("$(tup[begin]).$(i-1)", tup[begin+1:end]...) for (i, tup) in enumerate(ns_ordered_tups)]
   else
       ns_ordered_tups
   end
end

function expand_ordered(pm_tups)
    ordered_tups = reduce(vcat, map(_to_ordered, pm_tups))
    return ordered_tups
end

function generate_classes_iterate_signs(nz_tup)
    pm_tups = expand_plusminus(nz_tup)
    ordered_tups = expand_ordered(pm_tups)
    return ordered_tups
end

function generate_all_lag_motif_classes()
    nonzero_indicators = [(n1, n2, t1, t2) for 
        t2 ∈ [false, true],
        t1 ∈ [false, true],
        n2 ∈ [false, true],
        n1 ∈ [false, true]
    ]
    labeled_nonzero_indicators = [(i-1, tup...) for (i,tup) in enumerate(nonzero_indicators)]
    reduce(vcat, map(generate_classes_iterate_signs, labeled_nonzero_indicators))
end