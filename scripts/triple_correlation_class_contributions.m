function [class_contribution,class_count,contribution,nonzero_contributions] = triple_correlation_class_contributions(arr, neuron_window, time_window)
    class_contribution = zeros(1,14);
    nonzero_contributions = zeros(1,14);
    class_count = zeros(1,14);
    
    neuron_lag_range = - floor((neuron_window / 2)):floor((neuron_window / 2));       
    time_lag_range = -floor((time_window / 2)):floor((time_window / 2));

    [N_neurons, N_times] = size(arr);

    for n1 = -floor((neuron_window / 2)):floor((neuron_window / 2)) ;
    for n2 = -floor((neuron_window / 2)):floor((neuron_window / 2)) ; 
    for t1 = -floor((time_window / 2)):floor((time_window / 2));
    for t2 = -floor((time_window / 2)):floor((time_window / 2));
        
        class = network_motif_classification(n1, n2, t1, t2);
        class_count(1,class) = class_count(1,class) + 1;
        contribution = 0; 
        
                for r = (1-min(neuron_lag_range)):(N_neurons-max(neuron_lag_range))
                for c = (1-min(time_lag_range)):(N_times-max(time_lag_range))

                    contribution=contribution + 1;
                    addition = arr(r,c)*arr(r+n1,c+t1)*arr(r+n2,c+t2);
                    
                    if addition ~=0 
                        nonzero_contributions(1,class) =  nonzero_contributions(1,class) + 1;
                    end
                end
                end
                
        class_contribution(1,class) = contribution+class_contribution(1,class);
        
    end
    end
    end
    end
    
end
