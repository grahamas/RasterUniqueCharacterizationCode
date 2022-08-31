function [expectation,unscaled_expectation,p] = theoretical_spiking_expectation_3D(raster, x_window, y_window, time_window)
    
    motif_order = [
            1
            2 
            3 
            2
            3
            2
            3
            3
            3
            3
            3
            3
            3
            3];    

    p = mean(raster,'all');
    scale = numel(raster);

    t_pm = (time_window);
    n_pm = (x_window+1)*(y_window+1) - 1;

    t_m = floor(t_pm / 2);
    t_p = ceil(t_pm / 2);

    triplet_count = [
        1  
        3*t_pm  
        t_pm*(t_pm-1)  
        3*n_pm  
        n_pm*(n_pm-1)  
        3*n_pm*t_pm  
        4*n_pm*t_p + 2*n_pm*t_m  
        4*n_pm*t_m + 2*n_pm*t_p  
        n_pm*t_p*(t_p-1) + 2*n_pm*t_m*t_p + n_pm*t_m*(t_m-1) 
        n_pm*(t_p)*(t_p-1) + n_pm*(t_m)*(t_m-1) + 2*n_pm*t_m*t_p 
        n_pm*t_m*(t_m-1) + 2*n_pm*t_p*t_m + n_pm*t_p*(t_p-1) 
        n_pm*(n_pm-1)*t_p + 2*n_pm*(n_pm-1)*t_m
        n_pm*(n_pm-1)*t_m + 2*n_pm*(n_pm-1)*t_p
        n_pm*(n_pm-1)*t_pm*(t_pm-1)
    ];


    expectation = numel(raster) .* (p .^ motif_order) .* triplet_count ./ scale;
    expectation = expectation';
    unscaled_expectation = numel(raster) .* (p .^ motif_order) .* triplet_count;
    unscaled_expectation = unscaled_expectation';
end
