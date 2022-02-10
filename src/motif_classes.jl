
include("roman_encode.jl")

function count_distinct(args...)
    return length(unique(args))
end

function _1_neuron_motif_classification(n1, n2, t1, t2)
    # All neurons are the same
    # Assume t1 <= t2
    n_distinct_times = count_distinct(0, t1, t2)

    if n_distinct_times == 1
        # All neurons and times are the same
        return 1
    elseif n_distinct_times == 2
        # All neurons are the same, two times
        return 2
    elseif n_distinct_times == 3
        # All neurons are the same, all times distinct
        return 3
    else
        error("Invalid number of distinct times: $n_distinct_times")
    end
end

function _2_neuron_motif_classification(n1, n2, t1, t2)
    # Two neurons are the same (two distinct neurons involved)
    # Assume t1 <= t2
    n_distinct_times = count_distinct(0, t1, t2)

    if n_distinct_times == 1
        # Two neurons are the same, all times are the same
        return 4
    elseif n_distinct_times == 2
        if (0 == t1 && 0 == n1) || (t1 == t2 && n1 == n2) || (0 == t2 && 0 == n2)
            # Base and first nodes are same; or first and second
            return 6
        elseif (t1 == 0) || (t2 < 0)
            # Synchrony is first
            return 7
        elseif (t2 == 0) || (t1 > 0)
            # Synchrony is second
            return 8
        else
            error("Shouldn't be here")
        end
    elseif n_distinct_times == 3
        if (n1 == 0)
            if (0 < t1)
                return 10
            elseif (t1 < 0)
                if t2 < 0
                    return 11
                elseif t2 > 0
                    return 10
                else
                    error("Shouldn't be here")
                end
            else
                error("Shouldn't be here")
            end
        elseif (n2 == 0)
            if (t2 < 0)
                return 9
            elseif (t2 > 0)
                if t1 > 0 # in between
                    return 11
                elseif t1 < 0
                    return 9
                else
                    error("Shouldn't be here")
                end
            else
                error("Shouldn't be here")
            end
        elseif (n1 == n2)
            if (0 < t1)
                return 9
            elseif (t1 < 0 < t2)
                return 11
            elseif (t2 < 0)
                return 10
            else
                error("Shouldn't be here")
            end
        else
            error("Shouldn't get here.")
        end
    else
        error("Invalid number of same times: $n_distinct_times")
    end
end

function _3_neuron_motif_classification(n1, n2, t1, t2)
    # All neurons are distinct
    # Assume t1 <= t2
    n_distinct_times = count_distinct(0, t1, t2)

    if n_distinct_times == 1
        # All neurons are distinct, times are the same
        return 5
    elseif n_distinct_times == 2
        if t1 == 0
            return 13
        elseif t2 == 0
            return 12
        elseif t1 == t2
            if t1 > 0
                return 12
            else
                return 13
            end
        else
            error("Shouldn't get here.")
        end 
    elseif n_distinct_times == 3
        return 14 
    else
        error("Invalid number of same times: $n_distinct_times")
    end
end

function info_flow_classify_lag_motif_class(n1, n2, t1, t2)
    n_distinct_neurons = count_distinct(0, n1, n2)

    n1, n2, t1, t2 = if t1 < t2
        (n1, n2, t1, t2)
    else
        (n2, n1, t2, t1)
    end
    # Assume below that t1 <= t2

    if n_distinct_neurons == 1
        # All neurons are the same
        return _1_neuron_motif_classification(n1, n2, t1, t2)
    elseif n_distinct_neurons == 2
        return _2_neuron_motif_classification(n1, n2, t1, t2)
    elseif n_distinct_neurons == 3
        return _3_neuron_motif_classification(n1, n2, t1, t2)
    else
        error("Invalid number of same neurons: $n_distinct_neurons")
    end
end

function _to_plusminus(tup, i_nz)
    if tup[i_nz] == false
        [(tup[begin:i_nz-1]..., 0, tup[i_nz+1:end]...)]
    else
        [(tup[begin:i_nz-1]..., 1, tup[i_nz+1:end]...),
         (tup[begin:i_nz-1]..., -1, tup[i_nz+1:end]...)]
    end
end

function expand_plusminus(nz_tup)
    # Expand nonzero indicator tuple into plus-minuses
    l_pm_tups = [nz_tup]
    for i_nz = 2:5
        new_tups = [_to_plusminus(tup, i_nz) for tup in l_pm_tups]
        l_pm_tups = reduce(vcat, new_tups)
    end
    return if length(l_pm_tups) > 1
         [("$(tup[begin]).$i", tup[begin+1:end]...) for (i, tup) in enumerate(l_pm_tups)]
    else
        [("$(tup[begin])", tup[begin+1:end]...) for (i, tup) in enumerate(l_pm_tups)]
    end
end

function _to_ordered(tup)
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
    labeled_nonzero_indicators = [(i, tup...) for (i,tup) in enumerate(nonzero_indicators)]
    reduce(vcat, map(generate_classes_iterate_signs, labeled_nonzero_indicators))
end

function info_flow_classify_lag_motif_classes(motif_classes::Vector{<:Tuple})
    classified_motif_classes = map(motif_classes) do class_row
        (class_row..., info_flow_classify_lag_motif_class(class_row[begin+1:end]...))
    end
end

