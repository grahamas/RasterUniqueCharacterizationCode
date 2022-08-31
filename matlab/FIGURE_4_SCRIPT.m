%% Triple Correlation Script - Figure 4

close all;
clearvars
clc;

set(0,'defaultAxesFontSize',30)

C=150;  % columns
R=150;  % rows

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

str = 'Local Dynamics';
for n=9:50:Rr;
for t=9:50:Cr;
     p(n,t)=1;
     p(n,t+2) = 1;
     p(n,t+4) = 1;
end;
end;


%create spike rate matched noise raster
raster = p;

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

for i = 1:n_iterations
    [pre_noise_raster{i}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i} = cat(2,pre_noise_raster{i},curr_noise_raster{i},post_noise_raster{i});
end

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    noisy_raster = noise_raster{i};
    [class_contribution_noise{i},class_count,contribution_noise] = triple_correlation_class_contributions_2D_periodic(noisy_raster, neuron_window, time_window);
    scale = numel(noisy_raster(:,max_time_lag+1:N_times-max_time_lag));
    actual_noise{i} = class_contribution_noise{i}./scale;
    [expectation_noise{i}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i},neuron_window, time_window);
    nt{i} = actual_noise{i}./expectation_noise{i};
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise');

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise)');
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 20;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Time')
ylabel('Neuron #')
title(str,'fontsize',60)
axis('off')

x = 0:13;
figure;
hold on; set(gca,'fontsize',30)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',60)
xlim([-1 14])
set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
%ylim([0 1])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
ax = gca;
ax.FontSize = 25;
xlabel('Motif-Class #','fontsize',30)
title(str,'fontsize',40)

%% Synchrony
clearvars
clc;

set(0,'defaultAxesFontSize',30)

C=150;  % columns
R=150;  % rows

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

str = 'Synchrony';
for n=9:50:Rr;
    for t=9:50:Cr;
         p(n,t)=1;
         p(n+2,t) = 1;
         p(n+4,t) = 1;
    end;
end;


%create spike rate matched noise raster
raster = p;

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

for i = 1:n_iterations
    [pre_noise_raster{i}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i} = cat(2,pre_noise_raster{i},curr_noise_raster{i},post_noise_raster{i});
end

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    noisy_raster = noise_raster{i};
    [class_contribution_noise{i},class_count,contribution_noise] = triple_correlation_class_contributions_2D_periodic(noisy_raster, neuron_window, time_window);
    scale = numel(noisy_raster(:,max_time_lag+1:N_times-max_time_lag));
    actual_noise{i} = class_contribution_noise{i}./scale;
    [expectation_noise{i}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i},neuron_window, time_window);
    nt{i} = actual_noise{i}./expectation_noise{i};
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise');

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise)');
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 20;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Time')
ylabel('Neuron #')
title(str,'fontsize',60)
axis('off')

x = 0:13;
figure;
hold on; set(gca,'fontsize',30)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',60)
xlim([-1 14])
set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
%ylim([0 1])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
ax = gca;
ax.FontSize = 25;
xlabel('Motif-Class #','fontsize',30)
title(str,'fontsize',40)

%% Feedforward
clearvars
clc;

set(0,'defaultAxesFontSize',30)

C=150;  % columns
R=150;  % rows

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

str = 'Feedforward';
for n=9:50:R;
    for t=9:50:C;
         p(n,t)=1;
         p(n+3,t+1) = 1;
         p(n+1,t+3) = 1;
    end;
end;


%create spike rate matched noise raster
raster = p;

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

for i = 1:n_iterations
    [pre_noise_raster{i}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i} = cat(2,pre_noise_raster{i},curr_noise_raster{i},post_noise_raster{i});
end

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    noisy_raster = noise_raster{i};
    [class_contribution_noise{i},class_count,contribution_noise] = triple_correlation_class_contributions_2D_periodic(noisy_raster, neuron_window, time_window);
    scale = numel(noisy_raster(:,max_time_lag+1:N_times-max_time_lag));
    actual_noise{i} = class_contribution_noise{i}./scale;
    [expectation_noise{i}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i},neuron_window, time_window);
    nt{i} = actual_noise{i}./expectation_noise{i};
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise');

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise)');
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 20;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Time')
ylabel('Neuron #')
title(str,'fontsize',60)
axis('off')

x = 0:13;
figure;
hold on; set(gca,'fontsize',30)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',60)
xlim([-1 14])
set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
%ylim([0 1])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
ax = gca;
ax.FontSize = 25;
xlabel('Motif-Class #','fontsize',30)
title(str,'fontsize',40)

%% Feedback
clearvars
clc;

set(0,'defaultAxesFontSize',30)

C=150;  % columns
R=150;  % rows

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

str = 'Feedback';
for n=9:50:Rr;
    for t=9:50:Cr;
         p(n,t)=1;
         p(n+3,t+1) = 1;
         p(n,t+2) = 1;
    end;
end;


%create spike rate matched noise raster
raster = p;

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

for i = 1:n_iterations
    [pre_noise_raster{i}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i} = cat(2,pre_noise_raster{i},curr_noise_raster{i},post_noise_raster{i});
end

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    noisy_raster = noise_raster{i};
    [class_contribution_noise{i},class_count,contribution_noise] = triple_correlation_class_contributions_2D_periodic(noisy_raster, neuron_window, time_window);
    scale = numel(noisy_raster(:,max_time_lag+1:N_times-max_time_lag));
    actual_noise{i} = class_contribution_noise{i}./scale;
    [expectation_noise{i}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i},neuron_window, time_window);
    nt{i} = actual_noise{i}./expectation_noise{i};
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise');

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise)');
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 20;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Time')
ylabel('Neuron #')
title(str,'fontsize',60)
axis('off')

x = 0:13;
figure;
hold on; set(gca,'fontsize',30)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',60)
xlim([-1 14])
set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
%ylim([0 1])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
ax = gca;
ax.FontSize = 25;
xlabel('Motif-Class #','fontsize',30)
title(str,'fontsize',40)

%% Divergence
clearvars
clc;

set(0,'defaultAxesFontSize',30)

C=150;  % columns
R=150;  % rows

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

str = 'Divergence';
for n=9:50:Rr;
    for t=9:50:Cr;
         p(n,t)=1;
         p(n+1,t+2) = 1;
         p(n+3,t+2) = 1;
    end;
end;


%create spike rate matched noise raster
raster = p;

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

for i = 1:n_iterations
    [pre_noise_raster{i}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i} = cat(2,pre_noise_raster{i},curr_noise_raster{i},post_noise_raster{i});
end

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    noisy_raster = noise_raster{i};
    [class_contribution_noise{i},class_count,contribution_noise] = triple_correlation_class_contributions_2D_periodic(noisy_raster, neuron_window, time_window);
    scale = numel(noisy_raster(:,max_time_lag+1:N_times-max_time_lag));
    actual_noise{i} = class_contribution_noise{i}./scale;
    [expectation_noise{i}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i},neuron_window, time_window);
    nt{i} = actual_noise{i}./expectation_noise{i};
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise');

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise)');
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 20;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Time')
ylabel('Neuron #')
title(str,'fontsize',60)
axis('off')

x = 0:13;
figure;
hold on; set(gca,'fontsize',30)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',60)
xlim([-1 14])
set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
%ylim([0 1])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
ax = gca;
ax.FontSize = 25;
xlabel('Motif-Class #','fontsize',30)
title(str,'fontsize',40)

%% Convergence
clearvars
clc;

set(0,'defaultAxesFontSize',30)

C=150;  % columns
R=150;  % rows

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

str = 'Convergence';
for n=9:50:Rr;
    for t=9:50:Cr;
         p(n,t+2)=1;
         p(n+1,t) = 1;
         p(n+3,t) = 1;
    end;
end;


%create spike rate matched noise raster
raster = p;

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

for i = 1:n_iterations
    [pre_noise_raster{i}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i} = cat(2,pre_noise_raster{i},curr_noise_raster{i},post_noise_raster{i});
end

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    noisy_raster = noise_raster{i};
    [class_contribution_noise{i},class_count,contribution_noise] = triple_correlation_class_contributions_2D_periodic(noisy_raster, neuron_window, time_window);
    scale = numel(noisy_raster(:,max_time_lag+1:N_times-max_time_lag));
    actual_noise{i} = class_contribution_noise{i}./scale;
    [expectation_noise{i}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i},neuron_window, time_window);
    nt{i} = actual_noise{i}./expectation_noise{i};
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise');

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise)');
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 20;

figure; 
hold on; ylabel('')
hold on;
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);
        end;
    end;
end;
xlabel('Time')
ylabel('Neuron #')
title(str,'fontsize',60)
axis('off')

x = 0:13;
figure;
hold on; set(gca,'fontsize',30)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',60)
xlim([-1 14])
set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
%ylim([0 1])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
ax = gca;
ax.FontSize = 25;
xlabel('Motif-Class #','fontsize',30)
title(str,'fontsize',40)
