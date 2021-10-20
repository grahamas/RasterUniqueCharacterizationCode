%% Triple Correlation Script

%% Figure 2
clear all;
clc;
close all;     

% Figures 2A-B
C=150;  % columns
R=150;  % rows

%define the frequency
f = 0.08;

% deterministic component
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R
    ph=0; %no phase shift
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end

%addition of noise
y=rand(R,C);
x=0*y+1*S;                        
 
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

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   
Wt = 7;
Wn = 7;

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

%define the lag window size
neuron_window = 15;
time_window = 15;

[class_contribution,class_count,contribution,nonzero_contributions] = triple_correlation_class_contributions(p, neuron_window, time_window);

Scale = ((C - 2*Wt)) * ((R - 2*Wn));
actual=nonzero_contributions./Scale;
potential=class_contribution./Scale;
actual_potential_ratio = actual./potential;

figure(2); hold on; subplot(1,2,1)
sgtitle('Figures 2A-B')
hold on;
for k=1:Rr;
    for tt=1:Cr; 
        if p(k,tt)~=0; 
            plot(tt,p(k,tt)-2*k,'k.','MarkerSize',12);
        end;
    end;
end;
axis([0 Cr -Rr 0]);
title({'Raster: f=0.08 AU'})
axis('off');
hold on; 

figure(2); hold on; subplot(1,2,2)
bar(actual)
set(gca, 'YScale', 'log')
xlim([0 15]);
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
ylim([10e-2 10e4])
set(gca,'YMinorTick','Off')

% Figures 2C-D
C=150;  % columns
R=150;  % rows

%define the frequency
f = 0.08;

% deterministic component
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R  
    ph=0; %no phase shift
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end

%addition of noise
y=rand(R,C);
x=1*y+0*S;                        
 
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

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   
Wt = 7;
Wn = 7;

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

%define the lag window size
neuron_window = 15;
time_window = 15;

[class_contribution,class_count,contribution,nonzero_contributions] = triple_correlation_class_contributions(p, neuron_window, time_window);

Scale = ((C - 2*Wt)) * ((R - 2*Wn));

%noise calculations
actual_noise=nonzero_contributions./Scale;
potential_noise=class_contribution./Scale;
actual_noise_ratio = actual./actual_noise;

figure(3); hold on; subplot(1,2,1)
sgtitle('Figures 2C-D')
hold on;
for k=1:Rr;
    for tt=1:Cr; 
        if p(k,tt)~=0; 
            plot(tt,p(k,tt)-2*k,'k.','MarkerSize',12);
        end;
    end;
end;
axis([0 Cr -Rr 0]);
title({'Noise'})
axis('off');
hold on; 

figure(3); hold on; subplot(1,2,2)
bar(actual_noise,'k')
set(gca, 'YScale', 'log')
xlim([0 15]);
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
ylim([10e-2 10e4])
set(gca,'YMinorTick','Off')
xlabel('Motif-Classes')

% Figures 2E-F
C=150;  % columns
R=150;  % rows

% deterministic component
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R  
    ph=0; %no phase shift
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end                    
 
% create all-firing raster
for n=1:R;
    for t=1:C;
        p(n,t)=1;
    end;
end;

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   
Wt = 7;
Wn = 7;

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

%define the lag window size
neuron_window = 15;
time_window = 15;

[class_contribution,class_count,contribution,nonzero_contributions] = triple_correlation_class_contributions(p, neuron_window, time_window);

Scale = ((C - 2*Wt)) * ((R - 2*Wn));
actual=nonzero_contributions./Scale;
potential=class_contribution./Scale;

figure(4); hold on; subplot(1,2,1)
sgtitle('Figures 2E-F')
hold on;
for k=1:Rr;
    for tt=1:Cr; 
        if p(k,tt)~=0; 
            plot(tt,p(k,tt)-2*k,'k.','MarkerSize',12);
        end;
    end;
end;
axis([0 Cr -Rr 0]);
title({'All-Firing Raster'})
axis('off');
hold on; 

figure(4); hold on; subplot(1,2,2)
bar(potential,'k')
set(gca, 'YScale', 'log')
xlim([0 15]);
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
ylim([10e-2 10e4])
set(gca,'YMinorTick','Off')
xlabel('Motif-Classes')
ylabel('Number of Triplets (Scaled Count)')

figure(5); hold on; subplot(1,2,1)
sgtitle('Figures 2G-H')
stem(actual_potential_ratio,'bo','LineWidth',7)
ylim([0 0.6])
yticks([0 0.2 0.4 0.6])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
set(gcf,'Position',[406   101   869   719])
ylabel('A/P')
xlabel('Motif-Classes')

figure(5); hold on; subplot(1,2,2)
stem(actual_noise_ratio,'ro','LineWidth',7)
ylim([0 5])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
set(gcf,'Position',[406   101   869   719])
ylabel('A/N')
xlabel('Motif-Classes')

%% Figure 3D
clear all;
clc;

C=150;  % columns
R=150;  % rows

% deterministic component
f = 0.12;
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R  
    ph=0;
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end

%addition of noise
y=rand(R,C);
x=1*y+0*S;                        
 
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

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   
Wt = 7;
Wn = 7;

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

sum_raster_activity = sum(p);

%define the lag window size
neuron_window = 15;
time_window = 15;

[class_contribution,class_count,contribution,nonzero_contributions] = triple_correlation_class_contributions(p, neuron_window, time_window);

Scale = ((C - 2*Wt)) * ((R - 2*Wn));
actual_noise=nonzero_contributions./Scale;
potential_noise=class_contribution./Scale;
actual_noise_ratio = actual_noise./actual_noise;
actual_potential_ratio = actual_noise./potential_noise;

figure(6); hold on; 
sgtitle('Figure 3D')
subplot(7,1,[1 2 3 4]); 
ylabel('')

hold on;
for k=1:Rr;
    for tt=1:Cr; 
        if p(k,tt)~=0; 
            plot(tt,p(k,tt)-2*k,'k.','MarkerSize',12);
        end;
    end;
end;
axis([0 Cr -Rr 0]);
title({'Noise'})
axis('off');
hold on; 
subplot(7,1,5);
plot(sum_raster_activity); hold on; area(sum_raster_activity,'FaceColor','k','FaceAlpha',1)
ylim([0 Rr])
xlim([0 Cr])
axis off
box off
subplot(7,1,6)
stem(actual_potential_ratio,'bo','LineWidth',7)
ylim([0 0.6])
yticks([0 0.2 0.4 0.6])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
set(gcf,'Position',[406   101   869   719])
ylabel('A/P')
xlabel('Motif-Classes')

subplot(7,1,7)
stem(actual_noise_ratio,'ro','LineWidth',7)
ylim([0 4])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
set(gcf,'Position',[406   101   869   719])
ylabel('A/N')
xlabel('Motif-Classes')

set(gcf,'Position',[406   101   869   719])

% Figure 3A
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
x=0*y+1*S;                        
 
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

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   
Wt = 7;
Wn = 7;

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

sum_raster_activity = sum(p);

%define the lag window size
neuron_window = 15;
time_window = 15;

[class_contribution,class_count,contribution,nonzero_contributions] = triple_correlation_class_contributions(p, neuron_window, time_window);

Scale = ((C - 2*Wt)) * ((R - 2*Wn));
actual=nonzero_contributions./Scale;
potential=class_contribution./Scale;
actual_potential_ratio = actual./potential;

figure(7); hold on; 
sgtitle('Figure 3A')
subplot(7,1,[1 2 3 4]); 
ylabel('')

hold on;
for k=1:Rr;
    for tt=1:Cr; 
        if p(k,tt)~=0; 
            plot(tt,p(k,tt)-2*k,'k.','MarkerSize',12);
        end;
    end;
end;
axis([0 Cr -Rr 0]);
title({'Raster: f=0.12 AU'})
axis('off');
hold on; 
subplot(7,1,5);
stem(sum_raster_activity,'k','LineWidth',4,'Marker','none')
ylim([0 Rr])
xlim([0 Cr])
axis off
box off
subplot(7,1,6)
stem(actual_potential_ratio,'bo','LineWidth',7)
ylim([0 0.6])
yticks([0 0.2 0.4 0.6])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
ylabel('A/P')
xlabel('Motif-Classes')

subplot(7,1,7)
stem(actual_noise_ratio,'ro','LineWidth',7)
ylim([0 4])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
ylabel('A/N')
xlabel('Motif-Classes')
set(gcf,'Position',[406   101   869   719])

% Figure 3B
C=150;  % columns
R=150;  % rows

%define the frequency
f = 0.12;

% deterministic component
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R  
    ph=6*pi*(ch/R); %add a phase shift
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end

%addition of noise
y=rand(R,C);
x=0*y+1*S;                        
 
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

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   
Wt = 7;
Wn = 7;

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

sum_raster_activity = sum(p);

%define the lag window size
neuron_window = 15;
time_window = 15;

[class_contribution,class_count,contribution,nonzero_contributions] = triple_correlation_class_contributions(p, neuron_window, time_window);

Scale = ((C - 2*Wt)) * ((R - 2*Wn));
actual=nonzero_contributions./Scale;
potential=class_contribution./Scale;
actual_potential_ratio = actual./potential;

figure(8); hold on; 
sgtitle('Figure 3B')
subplot(7,1,[1 2 3 4]); 
ylabel('')

hold on;
for k=1:Rr;
    for tt=1:Cr; 
        if p(k,tt)~=0; 
            plot(tt,p(k,tt)-2*k,'k.','MarkerSize',12);
        end;
    end;
end;
axis([0 Cr -Rr 0]);
title({'Phase Shift'})
axis('off');
hold on; 
subplot(7,1,5);
plot(sum_raster_activity); hold on; area(sum_raster_activity,'FaceColor','k','FaceAlpha',1)
ylim([0 Rr])
xlim([0 Cr])
axis off
box off
subplot(7,1,6)
stem(actual_potential_ratio,'bo','LineWidth',7)
ylim([0 0.6])
yticks([0 0.2 0.4 0.6])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
ylabel('A/P')
xlabel('Motif-Classes')

subplot(7,1,7)
stem(actual_noise_ratio,'ro','LineWidth',7)
ylim([0 4])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
ylabel('A/N')
xlabel('Motif-Classes')
set(gcf,'Position',[406   101   869   719])

% Figure 3C
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
x = x/(1+1); %scale signal to account for the addition of noise
 
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

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   
Wt = 7;
Wn = 7;

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

sum_raster_activity = sum(p);

%define the lag window size
neuron_window = 15;
time_window = 15;

[class_contribution,class_count,contribution,nonzero_contributions] = triple_correlation_class_contributions(p, neuron_window, time_window);

Scale = ((C - 2*Wt)) * ((R - 2*Wn));
actual=nonzero_contributions./Scale;
potential=class_contribution./Scale;
actual_potential_ratio = actual./potential;

figure(9); hold on; 
sgtitle('Figure 3C')
subplot(7,1,[1 2 3 4]); 
ylabel('')

hold on;
for k=1:Rr;
    for tt=1:Cr; 
        if p(k,tt)~=0; 
            plot(tt,p(k,tt)-2*k,'k.','MarkerSize',12);
        end;
    end;
end;
axis([0 Cr -Rr 0]);
title({'SNR = 0dB'})
axis('off');
hold on; 
subplot(7,1,5);
plot(sum_raster_activity); hold on; area(sum_raster_activity,'FaceColor','k','FaceAlpha',1)
ylim([0 Rr])
xlim([0 Cr])
axis off
box off
subplot(7,1,6)
stem(actual_potential_ratio,'bo','LineWidth',7)
ylim([0 0.6])
yticks([0 0.2 0.4 0.6])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
ylabel('A/P')
xlabel('Motif-Classes')

subplot(7,1,7)
stem(actual_noise_ratio,'ro','LineWidth',7)
ylim([0 4])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
ylabel('A/N')
xlabel('Motif-Classes')
set(gcf,'Position',[406   101   869   719])
