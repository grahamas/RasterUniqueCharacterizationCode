% master program starts data specific standalon programs that are specific
% to the data file. it uses standalone_process replicated for each data
% file. then run the program in the background. when no file is selected
% the master program will end. however, standalone_process will run in the
% background.

function rat_master_loop_program_0714(start_num, end_num)

    trial_num = str2double(start_num);
    end_trial_num = str2double(end_num);

    active_programs = 0;
    % get the start trail number as input
    %trial_num = input("what is the START trial_num ?");
    %end_trial_num = input("what is the END trial_num ?");
    if isempty(trial_num)
        return
    end

    % open machine specifc log file for master loop
    mlog_name = [getenv('COMPUTERNAME') '_master_loop.log'];
    mlog_FID = fopen(mlog_name, 'a');
    fprintf(mlog_FID, '<<< MASTER LOOP START %s \n', datetime);

    while 1
        
        % exit if end trial num reached
        if trial_num > end_trial_num
            fprintf(mlog_FID, '>>> MASTER LOOP END %s \n', datetime);
            fclose(mlog_FID);
            return
        end

        %check if any program completed. if it for the first time, then there
        %will be none with RESULTS found. So run 8 programs and wait
    
        filelist = {dir('RESULTS_*.mat').name};
        numoffiles = length(filelist);
          
        if isequal (numoffiles, 0) && (active_programs < 8)
            % run 8 concurrent processes
            for i=1:8
                r_status = run_the_program(trial_num, mlog_FID);
                if isequal(r_status, 0)
                    fprintf(mlog_FID, '>>> MASTER LOOP END %s \n', datetime);
                    fclose(mlog_FID);
                    return
                end
                active_programs = active_programs+1;
                trial_num = trial_num+1;
            end
        else
    
            % wait for one of the programs to complete and run the next data
            % file
            % it is possible no program could have been completed. so
            % numoffiles could be 0 but active_programs could be 9
            if numoffiles > 0
                % for each completed RESULTS file, run another program
                for fn=1:numoffiles
                    % move the RESULTS file to COMPLETED file name
                    resultfile = filelist{fn};
                    completedfile = strrep(resultfile, 'RESULTS','COMPLETED');
                    movefile(resultfile, completedfile);

                    r_status = run_the_program(trial_num, mlog_FID);
                    if isequal(r_status, 0)
                        fprintf(mlog_FID, '>>> MASTER LOOP END %s \n', datetime);
                        fclose(mlog_FID);
                        return
                    end
                    trial_num = trial_num+1;
                end
            end
        end
        % sleep for 30 minutes (30*60 = 1800 seconds)and try to check for completed
        % program
        pause(30);
    end
end

%% function to run the program given trail_num

function [r_status] = run_the_program(trial_num, mlog_FID)
  
    r_status = 0;
    datafilename = ['*trial' num2str(trial_num) '.mat'];
    datalist = {dir(datafilename).name};
    if isempty(datalist)==1
        return;
    end

    r_status = 1;
    datafile = datalist{1};
%     fullfilename = fullfile(inputfpath, inputfname);
%     fprintf('selected file %s \n', fullfilename);

    %let us split full name path into path, filename and extension.
    [fpath, fname, fext] = fileparts(datafile);

    %let us create a copy of standalone process that takes this data file
    %make sure standalone_process.exe is the same folder as data file
    standaloneprocinstance = ['rat_cultures_standalone_loop_process_0714.exe'];
    %fprintf('proc file %s \n', standaloneprocinstance);

    datainstanceproc = append(fname,'.exe');
    %fprintf('new proc file %s \n', datainstanceproc);

    %create a copy of standalone process to match data file
    cpcmd = append('copy ', standaloneprocinstance, ' ', datainstanceproc);
    disp(cpcmd);
    %create a copy
    status = system(cpcmd);
    disp(status);

    %now run the datainstance proc in the background with data file as
    %param
    execmd = append(datainstanceproc, ' ', datafile, ' &');
    fprintf(mlog_FID, 'At %s :: Start %s \n', datetime, execmd);

    disp(execmd);
    disp('\n');
    sysstatus = system(execmd);
end

%% 
