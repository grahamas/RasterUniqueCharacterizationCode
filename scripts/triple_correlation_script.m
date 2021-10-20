clear all;
clc;
close all;

load ihot.mat                       

C=150;  % columns
R=150;  % rows

%define the frequency
f = 0.08;

% deterministic component
t=0:C-1;
s=1*sin(2*pi*f*t);
S=zeros(R,C);

for ch=1:R
    %ph=6*pi*(ch/R);
    % OR ... pick one phase shift   
    ph=0;
    S(ch,:)=1*(((sin(2*pi*f*t+ph)+1)/2));
end

%addition of noise
y=rand(R,C);
x=0*y+1*S;                        
%x = x/(1+1); %scale the noise, if needed
 
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

% create all-firing raster (Fig 2E)
% for n=1:R;
%     for t=1:C;    
%         p(n,t)=1;
%     end;
% end;

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   

Wt = 7;
Wn = 7;
C = 150; 
R = 150;
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

%noise calculations
% actual_noise=nonzero_contributions./Scale;
% potential_noise=class_contribution./Scale;
%actual_noise_ratio = actual./actual_noise;

%% Generate figures
figure; hold on;
subplot(5,1,[1 2 3 4]); 
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
title({'Raster: f=0.08 AU'})
axis('off');
hold on; 
subplot(5,1,5);
stem(sum_raster_activity,'k','LineWidth',4,'Marker','none')
%plot(sum_raster_activity); hold on; area(sum_raster_activity,'FaceColor','k','FaceAlpha',1)
ylim([0 Rr])
xlim([0 Cr])
axis('off')
box off
set(gcf,'Position',[406   101   869   719])

figure; hold on;
hold on;
bar(potential,'k')
set(gca, 'YScale', 'log')
xlim([0 15]);
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
xlabel('Motif-Class #')
ylabel('Number of Triplets (Scaled Count)')
title('Potential Contributions')
ylim([10e-2 10e4])
set(gca,'YMinorTick','Off')

figure; hold on; 
stem(actual_potential_ratio,'bo','LineWidth',7)
ylim([0 0.6])
yticks([0 0.2 0.4 0.6])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
set(gcf,'Position',[406   101   869   719])
ylabel('A/P')
xlabel('Motif-Classes')
