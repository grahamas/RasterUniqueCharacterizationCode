
include(srcdir("main.jl"))

let N_neurons = 150, N_times = 150, window_size = 15, transform=log10;

# pattern_strings_by_property = Dict(
#     "Level of Activity" => "I, II, IV",
#     "Local Dynamics" => "II, III, VII, VIII, IX, X, XI",
#     "Synchrony" => "V, VII, VIII, XII, XIII",
#     "Network Propagation" => "VI, VII, VIII, IX, X, XI, XII, XIII, XIV",
#     "Convergence" => "VII, IX, XIII, XIV",
#     "Divergence" => "VIII, XI, XII, XIV",
#     "Feedback" => "X",
#     "Feedforward" => "XIV"
# )
# pattern_list_by_property = Dict()
# for (property, str_patterns) ∈ pairs(pattern_strings_by_property)
#     pattern_list_by_property[property] = strip.(split(str_patterns, ","))
# end
# property_list_by_pattern = Dict()
# # Invert mapping
# for (property, patterns_list) ∈ pairs(pattern_list_by_property)
#     for this_pattern ∈ patterns_list
#         this_property_list = get!(property_list_by_pattern, this_pattern, String[])
#         push!(this_property_list, property)
#     end
# end

# property_list_string_by_pattern = Dict(
#     "I" => "Level of Activity",
#     "II" => "Level of Activity, Local Dynamics",
#     "III" => "Local Dynamics",
#     "IV" => "Level of Activity, Synchrony",
#     "V" => "Synchrony",
#     "VI" => "Level of Activity, Feedforward",
#     "VII" => "Local Dynamics, Synchrony, Feedforward",
#     "VIII" => "Local Dynamics, Synchrony, Feedforward",
#     "IX" => "Local Dynamics, Feedforward",
#     "X" => "Local Dynamics, Feedforward",
#     "XI" => "Local Dynamics, Feedforward, Feedback",
#     "XII" => "Synchrony, Feedforward, Divergence  ",
#     "XIII" => "Synchrony, Feedforward, Convergence  ",
#     "XIV" => "Feedforward, Convergence, Divergence"
# )

property_list_string_by_pattern = Dict(
    "I" => "Firing",
    "II" => "Firing, Local Dynamics",
    "III" => "Firing, Local Dynamics",
    "IV" => "Firing, Synchrony",
    "V" => "Firing, Synchrony",
    "VI" => "Firing, Synapse",
    "VII" => "Firing, Local Dynamics, Synchrony, Synapse",
    "VIII" => "Firing, Local Dynamics, Synchrony, Synapse",
    "IX" => "Firing, Local Dynamics, Synapse",
    "X" => "Firing, Local Dynamics, Synapse",
    "XI" => "Firing, Local Dynamics, Synapse, Feedback",
    "XII" => "Firing, Synchrony, Synapse, Divergence  ",
    "XIII" => "Firing, Synchrony, Synapse, Convergence  ",
    "XIV" => "Firing, Synapse, Feedforward, Convergence, Divergence"
)
property_list_by_pattern = Dict()
for (pattern, property_list_string) ∈ pairs(property_list_string_by_pattern)
    property_list_by_pattern[pattern] = strip.(split(property_list_string, ","))
end
# Invert mapping
patterns_by_property = Dict()
for (pattern, property_list) ∈ pairs(property_list_by_pattern)
    for this_property ∈ property_list
        this_pattern_list = get!(patterns_by_property, this_property, String[])
        push!(this_pattern_list, pattern)
    end
end

characteristic_patterns_string_by_property = Dict(
    "Local Instantaneous Activity" => "I",
    "Feedback" => "XI",
    "Local Dynamics" => "II, III",
    "Synchrony" => "IV, V",
    "Divergence" => "XII, XIV",
    "Convergence" => "XIII, XIV",
    "Spike Propagation" => "VI",
    "Feedforward" => "XIV"
)
# Turn the strings into lists
characteristic_patterns_by_property = Dict(
    key => strip.(split(val, ",")) for (key, val) in pairs(characteristic_patterns_string_by_property)
)
property_list_by_characteristic_pattern = Dict()
# Invert mapping
for (property, patterns_list) ∈ pairs(characteristic_patterns_by_property)
    for this_pattern ∈ patterns_list
        this_property_list = get!(property_list_by_characteristic_pattern, this_pattern, String[])
        push!(this_property_list, property)
    end
end


neuron_window_size = time_window_size = window_size
raster = BitMatrix(ones(N_neurons, N_times))

network_class_contributions = triple_correlation_network_classifications(raster, neuron_window_size, time_window_size)

contributions_by_characteristic_property = Dict{String,Float64}([(property, 0.) for property in keys(characteristic_patterns_by_property)])
total_contribution = 0
for (i_pattern, contribution) ∈ enumerate(network_class_contributions)
    total_contribution += contribution
    pattern_roman_numeral = roman_encode(i_pattern)
    property_list = get(property_list_by_characteristic_pattern, pattern_roman_numeral, String[])
    for property ∈ property_list
        contributions_by_characteristic_property[property] += contribution
    end
end

for property ∈ keys(contributions_by_characteristic_property)
    contributions_by_characteristic_property[property] /= total_contribution
end

display(contributions_by_characteristic_property)

contributions_by_property = Dict{String,Float64}([(property, 0.) for property in keys(patterns_by_property)])
total_contribution = 0
for (i_pattern, contribution) ∈ enumerate(network_class_contributions)
    total_contribution += contribution
    pattern_roman_numeral = roman_encode(i_pattern)
    property_list = get(property_list_by_pattern, pattern_roman_numeral, String[])
    for property ∈ property_list
        contributions_by_property[property] += contribution
    end
end

for property ∈ keys(contributions_by_property)
    contributions_by_property[property] /= total_contribution
end

display(contributions_by_property)

end