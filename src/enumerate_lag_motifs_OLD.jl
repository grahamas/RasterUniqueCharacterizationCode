
function _to_plusminus_OLD(tup, i_nz)
    if tup[i_nz] == false
        [(tup[begin:i_nz-1]..., 0, tup[i_nz+1:end]...)]
    else
        [(tup[begin:i_nz-1]..., 1, tup[i_nz+1:end]...),
         (tup[begin:i_nz-1]..., -1, tup[i_nz+1:end]...)]
    end
end

function expand_plusminus_OLD(nz_tup)
    # Expand nonzero indicator stuple into plus-minuses
    l_pm_tups = [nz_tup]
    for i_nz = 2:5
        new_tups = [_to_plusminus_OLD(tup, i_nz) for tup in l_pm_tups]
        l_pm_tups = reduce(vcat, new_tups)
    end
    return if length(l_pm_tups) > 1
         [("$(tup[begin]).$i", tup[begin+1:end]...) for (i, tup) in enumerate(l_pm_tups)]
    else
        [("$(tup[begin])", tup[begin+1:end]...) for (i, tup) in enumerate(l_pm_tups)]
    end
end

function _to_ordered_OLD(tup)
    ns_ordered_tups = if tup[2] == tup[3] == 1
        [
            ("$(tup[begin])", 2, 1, tup[end-1:end]...),
            ("$(tup[begin])", 1, 2, tup[end-1:end]...),
            ("$(tup[begin])", 1, 1, tup[end-1:end]...),
        ]
    elseif tup[2] == tup[3] == -1
        [
            ("$(tup[begin])", -2, -1, tup[end-1:end]...),
            ("$(tup[begin])", -1, -2, tup[end-1:end]...),
            ("$(tup[begin])", -1, -1, tup[end-1:end]...),
        ]
    else
        [tup]
    end
    ns_ordered_tups = map(ns_ordered_tups) do tup
        if tup[4] == tup[5] == 1
            [
                ("$(tup[begin])", tup[begin+1:begin+2]..., 2, 1),
                ("$(tup[begin])", tup[begin+1:begin+2]..., 1, 2),
                ("$(tup[begin])", tup[begin+1:begin+2]..., 1, 1),
            ]
        elseif tup[4] == tup[5] == -1
            [
                ("$(tup[begin])", tup[begin+1:begin+2]..., -2, -1),
                ("$(tup[begin])", tup[begin+1:begin+2]..., -1, -2),
                ("$(tup[begin])", tup[begin+1:begin+2]..., -1, -1),
            ]
        else
            tup
        end
    end 
    if !(eltype(ns_ordered_tups) <: Tuple)
        ns_ordered_tups = reduce(vcat, ns_ordered_tups)
    end

    if length(ns_ordered_tups) != 1
        [("$(tup[begin]).$i", tup[begin+1:end]...) for (i, tup) in enumerate(ns_ordered_tups)]
   else
       ns_ordered_tups
   end
end

function expand_ordered_OLD(pm_tups)
    ordered_tups = reduce(vcat, map(_to_ordered_OLD, pm_tups))
    return ordered_tups
end

function generate_classes_iterate_signs_OLD(nz_tup)
    pm_tups = expand_plusminus_OLD(nz_tup)
    ordered_tups = expand_ordered_OLD(pm_tups)
    return ordered_tups
end

function generate_all_lag_motif_classes_OLD()
    nonzero_indicators = [(n1, n2, t1, t2) for 
        t2 ∈ [false, true],
        t1 ∈ [false, true],
        n2 ∈ [false, true],
        n1 ∈ [false, true]
    ]
    labeled_nonzero_indicators = [(i, tup...) for (i,tup) in enumerate(nonzero_indicators)]
    reduce(vcat, map(generate_classes_iterate_signs_OLD, labeled_nonzero_indicators))
end

function info_flow_classify_lag_motif_classes(motif_classes::Vector{<:Tuple})
    classified_motif_classes = map(motif_classes) do class_row
        (class_row..., info_flow_classify_lag_motif_class(class_row[begin+1:end]...))
    end
end

