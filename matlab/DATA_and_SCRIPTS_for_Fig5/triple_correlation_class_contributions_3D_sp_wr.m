function [class_contribution] = triple_correlation_class_contributions_3D_sp_wr(arr, x_window,y_window, time_window)
    
    class_contribution = zeros(1,14);
    class_count = zeros(1,14);
   
    max_time_lag = ceil(time_window / 2);

    [N_x, N_y, N_times] = size(arr);

    for x1 = -floor((x_window / 2)):ceil((x_window / 2)) 
    for y1 = -floor((y_window / 2)):ceil((y_window / 2)) 
         
        spatial_wrap1 = circshift(arr,[-x1,-y1,0]);

    for x2 = -floor((x_window / 2)):ceil((x_window / 2))         
    for y2 = -floor((y_window / 2)):ceil((y_window / 2))  
    
        spatial_wrap2 = circshift(arr,[-x2,-y2,0]);

    for t1 = -floor((time_window / 2)):ceil((time_window / 2))
    for t2 = -floor((time_window / 2)):ceil((time_window / 2))
        
        class = network_motif_classification_3D(x1, x2, y1, y2, t1, t2);
        class_count(1,class) = class_count(1,class) + 1;
        contribution = 0; 
        
            for xr = 1:N_x
            for yr = 1:N_y
            for c = max_time_lag+1:N_times-max_time_lag

                contribution=contribution + arr(xr,yr,c)*spatial_wrap1(xr,yr,c+t1)*spatial_wrap2(xr,yr,c+t2);

            end
            end
            end
                
        class_contribution(1,class) = contribution+class_contribution(1,class);
        
    end
    end
    end
    end
    end
    end
    
end


