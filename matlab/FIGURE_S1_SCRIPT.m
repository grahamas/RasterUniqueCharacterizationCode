%% make rasters
clearvars
close all
clc

set(0,'defaultAxesFontSize',30)

%Baseline
load Rat_MEA2Baseline_DIV22_raster_3D_B7_trial13.mat %use for paper

raster_2D(1,:) = raster_3D(1,1,:);
raster_2D(2,:) = raster_3D(1,2,:);
raster_2D(3,:) = raster_3D(1,3,:);
raster_2D(4,:) = raster_3D(1,4,:);
raster_2D(5,:) = raster_3D(2,1,:);
raster_2D(6,:) = raster_3D(2,2,:);
raster_2D(7,:) = raster_3D(2,3,:);
raster_2D(8,:) = raster_3D(2,4,:);
raster_2D(9,:) = raster_3D(3,1,:);
raster_2D(10,:) = raster_3D(3,2,:);
raster_2D(11,:) = raster_3D(3,3,:);
raster_2D(12,:) = raster_3D(3,4,:);
raster_2D(13,:) = raster_3D(4,1,:);
raster_2D(14,:) = raster_3D(4,2,:);
raster_2D(15,:) = raster_3D(4,3,:);
raster_2D(16,:) = raster_3D(4,4,:);

fontsize=24;
linewidth = 10;
markersize =1.5;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster_2D,1);
    for tt=1:30051; 
        if raster_2D(k,tt)~=0; 
            plot(tt,raster_2D(k,tt)-0.1*k,'ko','MarkerFaceColor','k','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Bin # (500 Hz)')
ylabel('Channel #')
ylim([-10 10])
xlim([0 30051])
axis off
set(gcf,'position',[366    20   619   898])

%CONTROL
load Rat_MEA2Pharma_DIV22_raster_3D_A3_trial2.mat %control, use for paper

raster_2D(1,:) = raster_3D(1,1,:);
raster_2D(2,:) = raster_3D(1,2,:);
raster_2D(3,:) = raster_3D(1,3,:);
raster_2D(4,:) = raster_3D(1,4,:);
raster_2D(5,:) = raster_3D(2,1,:);
raster_2D(6,:) = raster_3D(2,2,:);
raster_2D(7,:) = raster_3D(2,3,:);
raster_2D(8,:) = raster_3D(2,4,:);
raster_2D(9,:) = raster_3D(3,1,:);
raster_2D(10,:) = raster_3D(3,2,:);
raster_2D(11,:) = raster_3D(3,3,:);
raster_2D(12,:) = raster_3D(3,4,:);
raster_2D(13,:) = raster_3D(4,1,:);
raster_2D(14,:) = raster_3D(4,2,:);
raster_2D(15,:) = raster_3D(4,3,:);
raster_2D(16,:) = raster_3D(4,4,:);

fontsize=24;
linewidth = 10;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster_2D,1);
    for tt=1:30051; 
        if raster_2D(k,tt)~=0; 
            plot(tt,raster_2D(k,tt)-0.1*k,'ko','MarkerFaceColor','k','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Bin # (500 Hz)')
ylabel('Channel #')
ylim([-10 10])
xlim([0 30051])
axis off
set(gcf,'position',[366    20   619   898])

%GABA
load Rat_MEA2Pharma_DIV22_raster_3D_B2_trial8.mat %use for paper
raster_2D(1,:) = raster_3D(1,1,:);
raster_2D(2,:) = raster_3D(1,2,:);
raster_2D(3,:) = raster_3D(1,3,:);
raster_2D(4,:) = raster_3D(1,4,:);
raster_2D(5,:) = raster_3D(2,1,:);
raster_2D(6,:) = raster_3D(2,2,:);
raster_2D(7,:) = raster_3D(2,3,:);
raster_2D(8,:) = raster_3D(2,4,:);
raster_2D(9,:) = raster_3D(3,1,:);
raster_2D(10,:) = raster_3D(3,2,:);
raster_2D(11,:) = raster_3D(3,3,:);
raster_2D(12,:) = raster_3D(3,4,:);
raster_2D(13,:) = raster_3D(4,1,:);
raster_2D(14,:) = raster_3D(4,2,:);
raster_2D(15,:) = raster_3D(4,3,:);
raster_2D(16,:) = raster_3D(4,4,:);

fontsize=24;
linewidth = 10;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster_2D,1);
    for tt=1:30051; 
        if raster_2D(k,tt)~=0; 
            plot(tt,raster_2D(k,tt)-0.1*k,'ko','MarkerFaceColor','k','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Bin # (500 Hz)')
ylabel('Channel #')
ylim([-10 10])
xlim([0 30051])
axis off
set(gcf,'position',[366    20   619   898])

%CNQX
load Rat_MEA2Pharma_DIV22_raster_3D_E6_trial33.mat %CNQX, use for paper
raster_2D(1,:) = raster_3D(1,1,:);
raster_2D(2,:) = raster_3D(1,2,:);
raster_2D(3,:) = raster_3D(1,3,:);
raster_2D(4,:) = raster_3D(1,4,:);
raster_2D(5,:) = raster_3D(2,1,:);
raster_2D(6,:) = raster_3D(2,2,:);
raster_2D(7,:) = raster_3D(2,3,:);
raster_2D(8,:) = raster_3D(2,4,:);
raster_2D(9,:) = raster_3D(3,1,:);
raster_2D(10,:) = raster_3D(3,2,:);
raster_2D(11,:) = raster_3D(3,3,:);
raster_2D(12,:) = raster_3D(3,4,:);
raster_2D(13,:) = raster_3D(4,1,:);
raster_2D(14,:) = raster_3D(4,2,:);
raster_2D(15,:) = raster_3D(4,3,:);
raster_2D(16,:) = raster_3D(4,4,:);

fontsize=24;
linewidth = 10;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster_2D,1);
    for tt=1:30051; 
        if raster_2D(k,tt)~=0; 
            plot(tt,raster_2D(k,tt)-0.1*k,'ko','MarkerFaceColor','k','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Bin # (500 Hz)')
ylabel('Channel #')
ylim([-10 10])
xlim([0 30051])
axis off
set(gcf,'position',[366    20   619   898])


%D-AP5
load Rat_MEA2Pharma_DIV22_raster_3D_A7_trial6.mat %use for paper
raster_2D(1,:) = raster_3D(1,1,:);
raster_2D(2,:) = raster_3D(1,2,:);
raster_2D(3,:) = raster_3D(1,3,:);
raster_2D(4,:) = raster_3D(1,4,:);
raster_2D(5,:) = raster_3D(2,1,:);
raster_2D(6,:) = raster_3D(2,2,:);
raster_2D(7,:) = raster_3D(2,3,:);
raster_2D(8,:) = raster_3D(2,4,:);
raster_2D(9,:) = raster_3D(3,1,:);
raster_2D(10,:) = raster_3D(3,2,:);
raster_2D(11,:) = raster_3D(3,3,:);
raster_2D(12,:) = raster_3D(3,4,:);
raster_2D(13,:) = raster_3D(4,1,:);
raster_2D(14,:) = raster_3D(4,2,:);
raster_2D(15,:) = raster_3D(4,3,:);
raster_2D(16,:) = raster_3D(4,4,:);

fontsize=24;
linewidth = 10;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster_2D,1);
    for tt=1:30051; 
        if raster_2D(k,tt)~=0; 
            plot(tt,raster_2D(k,tt)-0.1*k,'ko','MarkerFaceColor','k','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Bin # (500 Hz)')
ylabel('Channel #')
ylim([-10 10])
xlim([0 30051])
axis off
set(gcf,'position',[366    20   619   898])

%GABAZINE
load Rat_MEA2Pharma_DIV22_raster_3D_C8_trial21.mat %use for paper
raster_2D(1,:) = raster_3D(1,1,:);
raster_2D(2,:) = raster_3D(1,2,:);
raster_2D(3,:) = raster_3D(1,3,:);
raster_2D(4,:) = raster_3D(1,4,:);
raster_2D(5,:) = raster_3D(2,1,:);
raster_2D(6,:) = raster_3D(2,2,:);
raster_2D(7,:) = raster_3D(2,3,:);
raster_2D(8,:) = raster_3D(2,4,:);
raster_2D(9,:) = raster_3D(3,1,:);
raster_2D(10,:) = raster_3D(3,2,:);
raster_2D(11,:) = raster_3D(3,3,:);
raster_2D(12,:) = raster_3D(3,4,:);
raster_2D(13,:) = raster_3D(4,1,:);
raster_2D(14,:) = raster_3D(4,2,:);
raster_2D(15,:) = raster_3D(4,3,:);
raster_2D(16,:) = raster_3D(4,4,:);

fontsize=24;
linewidth = 10;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster_2D,1);
    for tt=1:30051; 
        if raster_2D(k,tt)~=0; 
            plot(tt,raster_2D(k,tt)-0.1*k,'ko','MarkerFaceColor','k','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Bin # (500 Hz)')
ylabel('Channel #')
ylim([-10 10])
xlim([0 30051])
axis off
set(gcf,'position',[366    20   619   898])
