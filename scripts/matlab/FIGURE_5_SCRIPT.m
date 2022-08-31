%%
clearvars
close all
clc

for well_num = 1:42
    get_wells = strcat('*COMPLETED_Rat_MEA2Pharma_DIV22_raster_3D*');

    filelist = {dir(get_wells).name};
    numoffiles = length(filelist);
    for i=1:numoffiles
        load(filelist{i});
        if size(actual_contribution,1)>=1
            % for valid file, extract trial number from the file name
            teststr = split(filelist{i},'3D_');
            teststr2 = split(teststr{2}, '_trial');

            if teststr2{1,1} == 'A2'
                well_num = 1;
            elseif teststr2{1,1} == 'A3'
                well_num = 2;
            elseif teststr2{1,1} == 'A4'
                well_num = 3;            
            elseif teststr2{1,1} == 'A5'
                well_num = 4;
            elseif teststr2{1,1} == 'A6'
                well_num = 5;     
            elseif teststr2{1,1} == 'A7'
                well_num = 6;
            elseif teststr2{1,1} == 'A8'
                well_num = 7;  
            elseif teststr2{1,1} == 'B2'
                well_num = 8;
            elseif teststr2{1,1} == 'B3'
                well_num = 9;            
            elseif teststr2{1,1} == 'B4'
                well_num = 10;
            elseif teststr2{1,1} == 'B5'
                well_num = 11;     
            elseif teststr2{1,1} == 'B6'
                well_num = 12;
            elseif teststr2{1,1} == 'B7'
                well_num = 13;               
            elseif teststr2{1,1} == 'B8'
                well_num = 14;  
            elseif teststr2{1,1} == 'C2'
                well_num = 15;
            elseif teststr2{1,1} == 'C3'
                well_num = 16;            
            elseif teststr2{1,1} == 'C4'
                well_num = 17;
            elseif teststr2{1,1} == 'C5'
                well_num = 18;     
            elseif teststr2{1,1} == 'C6'
                well_num = 19;
            elseif teststr2{1,1} == 'C7'
                well_num = 20;               
            elseif teststr2{1,1} == 'C8'
                well_num = 21; 
            elseif teststr2{1,1} == 'D2'
                well_num = 22;
            elseif teststr2{1,1} == 'D3'
                well_num = 23;            
            elseif teststr2{1,1} == 'D4'
                well_num = 24;
            elseif teststr2{1,1} == 'D5'
                well_num = 25;     
            elseif teststr2{1,1} == 'D6'
                well_num = 26;
            elseif teststr2{1,1} == 'D7'
                well_num = 27;               
            elseif teststr2{1,1} == 'D8'
                well_num = 28;   

            elseif teststr2{1,1} == 'E2'
                well_num = 29;
            elseif teststr2{1,1} == 'E3'
                well_num = 30;            
            elseif teststr2{1,1} == 'E4'
                well_num = 31;
            elseif teststr2{1,1} == 'E5'
                well_num = 32;     
            elseif teststr2{1,1} == 'E6'
                well_num = 33;
            elseif teststr2{1,1} == 'E7'
                well_num = 34;               
            elseif teststr2{1,1} == 'E8'
                well_num = 35;  

            elseif teststr2{1,1} == 'F2'
                well_num = 36;
            elseif teststr2{1,1} == 'F3'
                well_num = 37;            
            elseif teststr2{1,1} == 'F4'
                well_num = 38;
            elseif teststr2{1,1} == 'F5'
                well_num = 39;     
            elseif teststr2{1,1} == 'F6'
                well_num = 40;
            elseif teststr2{1,1} == 'F7'
                well_num = 41;               
            elseif teststr2{1,1} == 'F8'
                well_num = 42;                
            end

            % copy the variable
            actual_pharma(well_num,:) = cell2mat(actual_contribution(1))./(numel(snippet_raster));
            raw_actual_pharma(well_num,:) = cell2mat(actual_contribution(1));

            x_window = var_params(1,1);
            y_window = var_params(1,2);
            time_window = var_params(1,3);

            [expectation(well_num,:),unscaled_expectation(well_num,:),p] = theoretical_spiking_expectation_3D(snippet_raster, x_window, y_window, time_window);
    
            act_contribution(well_num,:) = cell2mat(actual_contribution);
            scaled_actual_contribution(well_num,:) = act_contribution(well_num,:)./(numel(snippet_raster));            
            [super_noise_pharma(well_num,:)] = expectation_conditioned_on_constituent_parts_3D(scaled_actual_contribution(well_num,:), snippet_raster,x_window, y_window, time_window);

            AT_pharma(well_num,:) = scaled_actual_contribution(well_num,:)./((super_noise_pharma(well_num,:)));
            
        end
    end

end


for well_num = 1:42
    get_wells = strcat('*COMPLETED_Rat_MEA2Baseline_DIV22_raster_3D*');

    filelist = {dir(get_wells).name};
    numoffiles = length(filelist);
    for i=1:numoffiles
        load(filelist{i});
        if size(actual_contribution,1)>=1
            % for valid file, extract trial number from the file name
            teststr = split(filelist{i},'3D_');
            teststr2 = split(teststr{2}, '_trial');

            if teststr2{1,1} == 'A2'
                well_num = 1;
            elseif teststr2{1,1} == 'A3'
                well_num = 2;
            elseif teststr2{1,1} == 'A4'
                well_num = 3;            
            elseif teststr2{1,1} == 'A5'
                well_num = 4;
            elseif teststr2{1,1} == 'A6'
                well_num = 5;     
            elseif teststr2{1,1} == 'A7'
                well_num = 6;
            elseif teststr2{1,1} == 'A8'
                well_num = 7;  
            elseif teststr2{1,1} == 'B2'
                well_num = 8;
            elseif teststr2{1,1} == 'B3'
                well_num = 9;            
            elseif teststr2{1,1} == 'B4'
                well_num = 10;
            elseif teststr2{1,1} == 'B5'
                well_num = 11;     
            elseif teststr2{1,1} == 'B6'
                well_num = 12;
            elseif teststr2{1,1} == 'B7'
                well_num = 13;               
            elseif teststr2{1,1} == 'B8'
                well_num = 14;  
            elseif teststr2{1,1} == 'C2'
                well_num = 15;
            elseif teststr2{1,1} == 'C3'
                well_num = 16;            
            elseif teststr2{1,1} == 'C4'
                well_num = 17;
            elseif teststr2{1,1} == 'C5'
                well_num = 18;     
            elseif teststr2{1,1} == 'C6'
                well_num = 19;
            elseif teststr2{1,1} == 'C7'
                well_num = 20;               
            elseif teststr2{1,1} == 'C8'
                well_num = 21; 
            elseif teststr2{1,1} == 'D2'
                well_num = 22;
            elseif teststr2{1,1} == 'D3'
                well_num = 23;            
            elseif teststr2{1,1} == 'D4'
                well_num = 24;
            elseif teststr2{1,1} == 'D5'
                well_num = 25;     
            elseif teststr2{1,1} == 'D6'
                well_num = 26;
            elseif teststr2{1,1} == 'D7'
                well_num = 27;               
            elseif teststr2{1,1} == 'D8'
                well_num = 28;   

            elseif teststr2{1,1} == 'E2'
                well_num = 29;
            elseif teststr2{1,1} == 'E3'
                well_num = 30;            
            elseif teststr2{1,1} == 'E4'
                well_num = 31;
            elseif teststr2{1,1} == 'E5'
                well_num = 32;     
            elseif teststr2{1,1} == 'E6'
                well_num = 33;
            elseif teststr2{1,1} == 'E7'
                well_num = 34;               
            elseif teststr2{1,1} == 'E8'
                well_num = 35;  

            elseif teststr2{1,1} == 'F2'
                well_num = 36;
            elseif teststr2{1,1} == 'F3'
                well_num = 37;            
            elseif teststr2{1,1} == 'F4'
                well_num = 38;
            elseif teststr2{1,1} == 'F5'
                well_num = 39;     
            elseif teststr2{1,1} == 'F6'
                well_num = 40;
            elseif teststr2{1,1} == 'F7'
                well_num = 41;               
            elseif teststr2{1,1} == 'F8'
                well_num = 42;                
            end

            % copy the variable
            actual_baseline(well_num,:) = cell2mat(actual_contribution(1))./(numel(snippet_raster));
            raw_actual_baseline(well_num,:) = cell2mat(actual_contribution(1)); 

            x_window = var_params(1,1);
            y_window = var_params(1,2);
            time_window = var_params(1,3);

            [expectation(well_num,:),unscaled_expectation(well_num,:),p] = theoretical_spiking_expectation_3D(snippet_raster, x_window, y_window, time_window);
    
            act_contribution = cell2mat(actual_contribution);
            scaled_actual_contribution = act_contribution./(numel(snippet_raster));            
            [super_noise(well_num,:)] = expectation_conditioned_on_constituent_parts_3D(scaled_actual_contribution, snippet_raster,x_window, y_window, time_window);

            AT_baseline(well_num,:) = scaled_actual_contribution./((super_noise(well_num,:)));

        end
    end

end



condition1 = [2,9,16,23,30,36,37]; %control
condition2 = [1,3,10,17,24,31,38]; %kainic acid
condition3 = [4,8,11,18,25,32,39]; %GABA
condition4 = [5,12,15,19,26,33,40]; %CNQX
condition5 = [6,13,20,22,27,34,41]; %DAP5
condition6 = [7,14,21,28,29,35,42]; %gabazine

y = [AT_pharma(condition1,:)./AT_baseline(condition1,:)]-1;
yy = [AT_pharma(condition2,:)./AT_baseline(condition2,:)]-1;
yyy = [AT_pharma(condition3,:)./AT_baseline(condition3,:)]-1;
yyyy = [AT_pharma(condition4,:)./AT_baseline(condition4,:)]-1;
yyyyy = [AT_pharma(condition5,:)./AT_baseline(condition5,:)]-1;
yyyyyy = [AT_pharma(condition6,:)./AT_baseline(condition6,:)]-1;

figure; hold on; title('Figure 5B')
x = {y,yyy,yyyyyy} ;
yline(0,'k')
motiflabels = {'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'};
boxplotGroup(x,'Colors','krb','interGroupSpace',7,'PrimaryLabels',{'','',''},'SecondaryLabels',motiflabels)
ylabel('$\frac{A/T Pharma}{A/T Baseline} - 1$','Interpreter','latex')
ylim([-1.5 3.5])

figure; hold on; title('Figure 5D')
x = {y,yyyyy} ;
yline(0,'k')
motiflabels = {'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'};
boxplotGroup(x,'Colors','kg','interGroupSpace',7,'PrimaryLabels',{'',''},'SecondaryLabels',motiflabels)
ylabel('$\frac{A/T Pharma}{A/T Baseline} - 1$','Interpreter','latex')

figure; hold on; title('Figure 5C')
x = {y,yyyy} ;
yline(0,'k')
motiflabels = {'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'};
boxplotGroup(x,'Colors','km','interGroupSpace',7,'PrimaryLabels',{'',''},'SecondaryLabels',motiflabels)
ylabel('$\frac{A/T Pharma}{A/T Baseline} - 1$','Interpreter','latex')
ylim([-1.5 2.5])

actual_control = raw_actual_pharma(condition1,:);
actual_kainic = raw_actual_pharma(condition2,:);
actual_GABA = raw_actual_pharma(condition3,:);
actual_CNQX = raw_actual_pharma(condition4,:);
actual_DAP5 = raw_actual_pharma(condition5,:);
actual_gabazine = raw_actual_pharma(condition6,:);

actual_control_baseline = raw_actual_baseline(condition1,:);
actual_kainic_baseline = raw_actual_baseline(condition2,:);
actual_GABA_baseline = raw_actual_baseline(condition3,:);
actual_CNQX_baseline = raw_actual_baseline(condition4,:);
actual_DAP5_baseline = raw_actual_baseline(condition5,:);
actual_gabazine_baseline = raw_actual_baseline(condition6,:);

var1 = (actual_control(:,1)./actual_control_baseline(:,1));
var2 = (actual_kainic(:,1)./actual_kainic_baseline(:,1));
var3 = (actual_GABA(:,1)./actual_GABA_baseline(:,1));
var4 = (actual_CNQX(:,1)./actual_CNQX_baseline(:,1));
var5 = (actual_DAP5(:,1)./actual_DAP5_baseline(:,1));
var6 = (actual_gabazine(:,1)./actual_gabazine_baseline(:,1));

figure; hold on; title('Figure 5A')
x = {var1, var3, var4, var5, var6} ;
yline(1,'k')
motiflabels = {'0'};
boxplotGroup(x,'Colors','krmgb','interGroupSpace',7,'PrimaryLabels',{'','','','',''},'SecondaryLabels',motiflabels)
ylabel('Spike count (% of baseline)')