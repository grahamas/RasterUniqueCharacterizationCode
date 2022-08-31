
% clearvars
% close all
% clc

function rat_cultures_standalone_loop_process_0714 (datafilename)
    
    % read the data file and set the raster variable
    raster = load(datafilename);
    raster = raster.raster_3D;

    % get file details, directory path, filename and extension - here it
    % shoud be .mat
    [fpath, fname, fext] = fileparts(datafilename);
    %let us create log file with same name, it will be text file
    logfilename = append(fname, '_log.txt');
    fidlog = fopen(logfilename, 'at+');

    %load parameters from a file

    var_params = readmatrix("parameter_file.txt");

%     x_window = 9;
%     y_window = 9;
%     time_window = 100;
%     n_noise_implementations = 2;
%     num_epochs = 0 means full, 1 means only one loop
%   max_time_lag (mtl_value) = 0 means calculate, 50 means use this value

    x_window = var_params(1);
    y_window = var_params(2);
    time_window = var_params(3);
    n_noise_implementations = var_params(4);
    num_epochs = var_params(5);
    el = var_params(6);

    fprintf(fidlog, 'Read paramters %s\n', datetime);
    fprintf(fidlog, '%s\n', mat2str(var_params));

    max_time_lag = ceil(time_window/2);

    if el == 0
        epoch_length = 30*500;
    else
        epoch_length = el;
    end

    size_snippet = 2*max_time_lag + epoch_length + 1;
    total_length = size(raster,3);
    
    n_epochs = floor(total_length / epoch_length )

    if (num_epochs == 1)
        n_epochs = 1;
    end
    
    actual_contribution=cell(n_epochs,1);
    actual_noise = cell(n_epochs,1);
    an_ratio = cell(n_epochs,1);

    init = 1;
    post_end = 0;

    fprintf(fidlog, 'Start %d n_epochs %s\n', n_epochs, datetime);
    for ii = init:n_epochs
        tic
        disp(ii)
        fprintf(fidlog, '%d Epoch %s\n', ii, datetime);

        if (post_end + 1 + size_snippet) < size(raster,3)

            slice_start = post_end+1;
    
            pre_snippet_raster = raster(:,:,slice_start:slice_start+max_time_lag-1);
            pre_end = slice_start+max_time_lag-1;
    
            snippet_raster = raster(:,:, pre_end+1: pre_end+1+epoch_length);
            snip_end = pre_end+1+epoch_length;
    
            post_snippet_raster = raster(:,:, snip_end+1 : snip_end+1+max_time_lag-1);
            post_end = snip_end+1+max_time_lag-1;
    
            temp_snippet = cat(3,pre_snippet_raster, snippet_raster, post_snippet_raster);
            
            size(temp_snippet)

            tic
            fprintf(fidlog, 'Starting triple correlation %s\n', datetime);
            [class_contribution]...
              = triple_correlation_class_contributions_3D_sp_wr(temp_snippet, x_window, ...
                            y_window, time_window);
            toc

            actual_contribution{ii} = class_contribution;

            % shuffle within each snippet separately
            tic
            fprintf(fidlog, 'Starting noise correlation %s\n', datetime);
            for jj = 1:n_noise_implementations
                clear noise_raster
                [pre_noise_raster{jj,1}] = generate_whole_raster_activity_matched_noise_3D(pre_snippet_raster);
                [curr_noise_raster{jj,1}] = generate_whole_raster_activity_matched_noise_3D(snippet_raster);
                [post_noise_raster{jj,1}] = generate_whole_raster_activity_matched_noise_3D(post_snippet_raster);
                
                noise_raster{jj,1} = cat(3,pre_noise_raster{jj,1},curr_noise_raster{jj,1},post_noise_raster{jj,1});      
                
                [noise_contribution{jj,1}]...
                    = triple_correlation_class_contributions_3D_sp_wr(noise_raster{jj,1}, x_window, ...
                            y_window, time_window);
            
            end

            toc
        end
        fprintf(fidlog, 'Completed %d Epoch %s \n', ii, datetime);
    end
    resultfile = append('RESULTS_', fname);
    save(resultfile);
    fprintf(fidlog, 'Result saved %s \n', datetime);
    fclose(fidlog);
end
