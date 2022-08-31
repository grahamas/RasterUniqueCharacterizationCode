function [noise_raster] = generate_whole_raster_activity_matched_noise_3D(raster)

    [x,y,t] = size(raster);
    
    noise_raster = zeros(x,y,t);
    
    num_spikes_raster= sum(sum(sum(raster)));
    
    noise_raster(randperm(numel(noise_raster),num_spikes_raster)) = 1;

end