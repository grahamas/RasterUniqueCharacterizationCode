%% Triple Correlation Script - Figure 3

clearvars
clc;
close all;

C=150;  % columns
R=150;  % rows

%define the frequency
f = 0.12;

% deterministic component
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R  
    ph=0;
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end

%addition of noise
y=rand(R,C);
x=1*y+1*S;                        
x = x/(1+1); %scale the noise

% create raster
for n=1:R;
    for t=1:C;
        if x(n,t)>.5;       
            p(n,t)=1;
        else;
            p(n,t)=0;
        end;
    end;
end;

%create spike rate matched noise raster
raster = p;
[n,t] = size(raster);

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

% triple correlation function

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    [pre_noise_raster{i,1}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i,1}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i,1}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i,1} = cat(2,pre_noise_raster{i,1},curr_noise_raster{i,1},post_noise_raster{i,1});
    noisy_raster = noise_raster{i,1};
    [class_contribution_noise{i,1},class_count,contribution_noise] = triple_correlation_class_contributions_2D_periodic(noisy_raster, neuron_window, time_window);
    scale = numel(curr_noise_raster{i,1});
    actual_noise{i,1} = class_contribution_noise{i,1}./scale;
    [expectation_noise{i,1}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i,1},neuron_window, time_window);
    nt{i,1} = (actual_noise{i,1})./(expectation_noise{i,1});
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise);

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise));
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);
a_minust = actual - conditioned_expectation;
n_minust = actual_noise - (cell2mat(expectation_noise));

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 10;
ylabelfontsize = 40;
xlabelfontsize = 40;
titlefontsize = 50;

x = 0:13;
figure; hold on; set(gca,'fontsize',fontsize)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',ylabelfontsize)
xlim([-1 14])
%ylim([-0.4 1.2])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
xlabel('Motif-Class #')
title('Figure 3B: SNR = 0 dB')

figure; 
ylabel('')
hold on;
for k=1:N_neurons;
    for tt=1:N_times; 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',10);
        end;
    end;
end;
title('Figure 3A: SNR = 0 dB','fontsize',50)
axis('off');


%% Figure 3C-D
clearvars
clc

C=150;  % columns
R=150;  % rows

%define the frequency
f = 0.12;

% deterministic component
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R  
    ph=0;
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end

%addition of noise
y=rand(R,C);
x=2.8*y+1*S;                        
x = x/(1+2.8); %scale the noise

% create raster
for n=1:R;
    for t=1:C;
        if x(n,t)>.5;       
            p(n,t)=1;
        else;
            p(n,t)=0;
        end;
    end;
end;

%create spike rate matched noise raster
raster = p;
[n,t] = size(raster);

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

% triple correlation function

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    [pre_noise_raster{i,1}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i,1}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i,1}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i,1} = cat(2,pre_noise_raster{i,1},curr_noise_raster{i,1},post_noise_raster{i,1});
    noisy_raster = noise_raster{i,1};
    [class_contribution_noise{i,1},class_count,contribution_noise] = triple_correlation_class_contributions_continuous_2D_spatial_wr(noisy_raster, neuron_window, time_window);
    scale = numel(curr_noise_raster{i,1});
    actual_noise{i,1} = class_contribution_noise{i,1}./scale;
    [expectation_noise{i,1}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i,1},neuron_window, time_window);
    nt{i,1} = (actual_noise{i,1})./(expectation_noise{i,1});
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise);

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise));
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);
a_minust = actual - conditioned_expectation;
n_minust = actual_noise - (cell2mat(expectation_noise));

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 10;
ylabelfontsize = 40;
xlabelfontsize = 40;
titlefontsize = 50;

x = 0:13;
figure; hold on; set(gca,'fontsize',fontsize)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',ylabelfontsize)
xlim([-1 14])
%ylim([-0.4 1.2])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
xlabel('Motif-Class #')
title('Figure 3D: SNR = -9 dB')

figure; 
ylabel('')
hold on;
for k=1:N_neurons;
    for tt=1:N_times; 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',10);
        end;
    end;
end;
title('Figure 3C: SNR = -9 dB','fontsize',50)
axis('off');

%% Figure 3E-F
clearvars
clc

C=150;  % columns
R=150;  % rows

%define the frequency
f = 0.12;

% deterministic component
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R  
    ph=0;
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end

%addition of noise
y=rand(R,C);
x=7*y+1*S;                        
x = x/(1+7); %scale the noise

% create raster
for n=1:R;
    for t=1:C;
        if x(n,t)>.5;       
            p(n,t)=1;
        else;
            p(n,t)=0;
        end;
    end;
end;

%create spike rate matched noise raster
raster = p;
[n,t] = size(raster);

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

% triple correlation function

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    [pre_noise_raster{i,1}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i,1}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i,1}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i,1} = cat(2,pre_noise_raster{i,1},curr_noise_raster{i,1},post_noise_raster{i,1});
    noisy_raster = noise_raster{i,1};
    [class_contribution_noise{i,1},class_count,contribution_noise] = triple_correlation_class_contributions_continuous_2D_spatial_wr(noisy_raster, neuron_window, time_window);
    scale = numel(curr_noise_raster{i,1});
    actual_noise{i,1} = class_contribution_noise{i,1}./scale;
    [expectation_noise{i,1}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i,1},neuron_window, time_window);
    nt{i,1} = (actual_noise{i,1})./(expectation_noise{i,1});
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise);

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise));
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);
a_minust = actual - conditioned_expectation;
n_minust = actual_noise - (cell2mat(expectation_noise));

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 10;
ylabelfontsize = 40;
xlabelfontsize = 40;
titlefontsize = 50;

x = 0:13;
figure; hold on; set(gca,'fontsize',fontsize)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',ylabelfontsize)
xlim([-1 14])
%ylim([-0.4 1.2])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
xlabel('Motif-Class #')
title('Figure 3F: SNR = -17 dB')

figure; 
ylabel('')
hold on;
for k=1:N_neurons;
    for tt=1:N_times; 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',10);
        end;
    end;
end;
title('Figure 3E: SNR = -17 dB','fontsize',50)
axis('off');

%% Figure G-H
clearvars
clc

C=150;  % columns
R=150;  % rows

%define the frequency
f = 0.12;

% deterministic component
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R  
    ph=0;
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end

%addition of noise
y=rand(R,C);
x=100*y+1*S;                        
x = x/(1+100); %scale the noise

% create raster
for n=1:R;
    for t=1:C;
        if x(n,t)>.5;       
            p(n,t)=1;
        else;
            p(n,t)=0;
        end;
    end;
end;

%create spike rate matched noise raster
raster = p;
[n,t] = size(raster);

neuron_window = 14;
time_window = 14;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 1;

% triple correlation function

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_2D_periodic(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    [pre_noise_raster{i,1}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i,1}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i,1}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    noise_raster{i,1} = cat(2,pre_noise_raster{i,1},curr_noise_raster{i,1},post_noise_raster{i,1});
    noisy_raster = noise_raster{i,1};
    [class_contribution_noise{i,1},class_count,contribution_noise] = triple_correlation_class_contributions_continuous_2D_spatial_wr(noisy_raster, neuron_window, time_window);
    scale = numel(curr_noise_raster{i,1});
    actual_noise{i,1} = class_contribution_noise{i,1}./scale;
    [expectation_noise{i,1}] = expectation_conditioned_on_constituent_parts_2D(actual_noise{i}, curr_noise_raster{i,1},neuron_window, time_window);
    nt{i,1} = (actual_noise{i,1})./(expectation_noise{i,1});
end

sum_raster_activity = sum(temp_raster);

[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_raster,neuron_window, time_window);
[unconditioned_expectation,unscaled_expected, p] = theoretical_expectation_2D(temp_raster, neuron_window, time_window);

actual_noise = cell2mat(actual_noise);

mean_actual_noise = mean(actual_noise,1);

anr_mean = actual./mean_actual_noise;

nt_noise = (actual_noise)./(cell2mat(expectation_noise));
mean_nt = mean(nt_noise,1);

at = actual./(conditioned_expectation);
a_minust = actual - conditioned_expectation;
n_minust = actual_noise - (cell2mat(expectation_noise));

for i = 1:14
    x = actual_noise(:,i);
    sem_noise(:,i) = std(x)/sqrt(length(x));
end

nt_noise1 = nt_noise-1;

fontsize=24;
linewidth = 10;
markersize = 10;
ylabelfontsize = 40;
xlabelfontsize = 40;
titlefontsize = 50;

x = 0:13;
figure; hold on; set(gca,'fontsize',fontsize)
stem(x,at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',ylabelfontsize)
xlim([-1 14])
%ylim([-0.4 1.2])
xticks(0:1:13)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
xlabel('Motif-Class #')
%title('Figure 3G: SNR = -40 dB')

figure; 
ylabel('')
hold on;
for k=1:N_neurons;
    for tt=1:N_times; 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',10);
        end;
    end;
end;
title('SNR = -40 dB','fontsize',50)
axis('off');