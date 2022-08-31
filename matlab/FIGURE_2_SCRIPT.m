%% Triple Correlation Script - Figure 2

clearvars
clc
close all
                    
C=150;  % columns
R=150;  % rows

%define the frequency
f = 0.08;

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
x=0*y+1*S;                        

str = 'Raster: f=0.08 AU';

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
neuron_window = 20;
time_window = 20;
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

n_iterations = 100;

% triple correlation function

[class_contribution,class_count,contribution] = triple_correlation_class_contributions_continuous_2D_spatial_wr(raster, neuron_window, time_window);
temp_raster = raster(:,max_time_lag+1:N_times-max_time_lag);
scale = numel(temp_raster);

actual=class_contribution./scale;

for i = 1:n_iterations
    disp(i)
    
    [pre_noise_raster{i,1}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i,1}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i,1}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    
    noise_raster{i,1} = cat(2,pre_noise_raster{i,1},curr_noise_raster{i,1},post_noise_raster{i,1});
    
    [class_contribution_noise{i,1},class_count,contribution_noise] = triple_correlation_class_contributions_continuous_2D_spatial_wr(noise_raster{i,1}, neuron_window, time_window);
    
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

x=1:1:14;
figure; hold on; set(gca,'fontsize',fontsize)
stem(at-1,'r','LineWidth',linewidth)
ylabel('A/T-1','fontsize',ylabelfontsize)
xlim([0 15])
ylim([-0.4 1.2])
xticks(1:1:14)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
xlabel('Motif-Class #')

x = 0:13;
figure; hold on; set(gca,'fontsize',fontsize)
yline(0,'k')
boxplot(nt_noise1)
ylabel('N/T-1','fontsize',ylabelfontsize)
xlim([0 15])
ylim([-0.06 0.06])
xticks(1:1:14)
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
xlabel('Motif-Class #','fontsize',40)

x = 1:14;
markersize = 120;
figure;hold on; set(gca,'fontsize',15)
scatter(x,actual,150,'r','filled','MarkerEdgeColor','r')
scatter(x,conditioned_expectation,120,'k','filled','MarkerEdgeColor','k','marker','d')
scatter(x,mean_actual_noise,'b','_')
boxplot(actual_noise,'colors','b')
legend('Actual Contributions (A)','Theoretical Expectation (T)','Noise Contributions (N, n=100)','Location','Northwest','Fontsize',20)
ylabel('Scaled Triplet Count','fontsize',30)
xlim([0 15])
xticks(1:1:14)
ylim([10e-2 10e4])
xticklabels({'0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII'})
xlabel('Motif-Class #','fontsize',30)
set(gca,'YScale','log')

figure; 
ylabel('')
hold on;
for k=1:N_neurons;
    for tt=1:N_times; 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',5);
        end;
    end;
end;
title('Raster: f=0.08 AU','fontsize',50)
axis('off');
hold on; 


noise_raster_figure = noise_raster{1,1};
figure; 
ylabel('')
hold on;
for k=1:N_neurons;
    for tt=1:N_times; 
        if noise_raster_figure(k,tt)~=0; 
            plot(tt,noise_raster_figure(k,tt)-2*k,'k.','MarkerSize',5);
        end;
    end;
end;
title('Noise','fontsize',titlefontsize)
axis('off');
hold on; 