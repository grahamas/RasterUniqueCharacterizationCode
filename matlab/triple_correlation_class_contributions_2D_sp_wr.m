function [class_contribution,class_count,contribution]= triple_correlation_class_contributions_2D_periodic(arr, neuron_window, time_window)
    class_contribution = zeros(1,14);
    class_count = zeros(1,14);
    
    max_time_lag = ceil((time_window / 2));
    [N_neurons, N_times] = size(arr);

    for n1 = -floor((neuron_window / 2)):ceil((neuron_window / 2)) 
    
        spatial_wrap1 = circshift(arr,[-n1,0]);

    for t1 = -floor((time_window / 2)):ceil((time_window / 2))

    for n2 = -floor((neuron_window / 2)):ceil((neuron_window / 2)) 
    
        spatial_wrap2 = circshift(arr,[-n2,0]);

    for t2 = -floor((time_window / 2)):ceil((time_window / 2))
        
        class = network_motif_classification(n1, n2, t1, t2);
        class_count(1,class) = class_count(1,class) + 1;
        contribution = 0;
        
                for r = 1:N_neurons
                for c = max_time_lag+1:N_times-max_time_lag

                    contribution = contribution + arr(r,c)*spatial_wrap1(r,c+t1)*spatial_wrap2(r,c+t2);

                end
                end
                
        class_contribution(1,class) = contribution+class_contribution(1,class);
        
    end
    end
    end
    end
    
end
