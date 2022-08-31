function [noise_raster] = generate_noise_raster(raster)
    [n,t] = size(raster);
    
    noise_raster = zeros(n,t);
    num_spikes_raster= sum(sum(sum(raster)));
    noise_raster(randperm(numel(noise_raster),num_spikes_raster)) = 1;
end