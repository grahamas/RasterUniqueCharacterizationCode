%% Triple Correlation Script

clearvars
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

% Only analyze part of the image
% determined by Rr and Cr
T1=7; 
N1=7;   
T2=7;  
N2=7;   

Cr=C-max(T1,T2);
Rr=R-max(N1,N2);
 
neuron_window = 15;
time_window = 15;

[class_contribution,class_count,contribution,nonzero_contributions,actual_potential_ratio] = triple_correlation_class_contributions(p, neuron_window, time_window);

sum_raster_activity = sum(p);

figure(99); hold on;
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
figure(99); hold on
set(gcf,'Position',[406   101   869   719])


Wt = 7;
Wn = 7;
C = 150; 
R = 150;

Scale = ((C - Wt) -(1+Wt) + 1) * ((R - Wn) -(1+Wn) + 1);
actual=nonzero_contributions./Scale;
potential=class_contribution./Scale;
actual_potential_ratio = actual./potential;

% actual_noise=nonzero_contributions./Scale;
% potential_noise=class_contribution./Scale;

actual_noise_ratio = actual./actual_noise;

figure; hold on; subplot(2,1,1)
stem(actual_potential_ratio,'bo','LineWidth',7)
ylim([0 0.6])
yticks([0 0.2 0.4 0.6])
xlim([0 15])
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})

% title('Actual Contribution / Potential Contribution')

set(gcf,'Position',[406   101   869   719])

subplot(2,1,2);
stem(actual_noise_ratio,'ro','MarkerEdgeColor','r','LineWidth',7)
ylim([0 4.2])
yticks([0 2 4])
xlim([0 15]);
xlim([0 15]);
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
set(gcf,'Position',[406   101   869   719])

% title('Actual Contribution / Noise Contribution')
% ylabel('Ratio')
% xlabel('Motif-Classes')

%legend('Actual Contribution / Noise','Actual Contribution / Potential Contribution','Location','NorthEast')

% 
% figure(2); hold on;
% subplot(2,2,2);
% hold on;
% bar(potential,'k')
% grid on
% set(gca, 'YScale', 'log')
% xlim([0 15]);
% xticks([1 5 10 14])
% xticklabels({'I','V','X','XIV'})
% xlabel('Motif-Class #')
% ylabel('Number of Triplets (Scaled)')
% title('Potential Contributions')
% grid on
% ylim([10e-2 10e5])
% % 
% figure(2); hold on;
% subplot(2,2,3);
% hold on;
% bar(actual)
% set(gca, 'YScale', 'log')
% xlim([0 15]);
% xticks([1 5 10 14])
% xticklabels({'I','V','X','XIV'})
% xlabel('Motif-Class #')
% ylabel('Number of Triplets (Scaled)')
% title('Actual Contributions')
% grid on
% ylim([10e-2 10e5])
% 
% %%
% % 
figure; hold on;
ylabel('')

axis([0 Cr 0 Rr]);
hold on;
for k=1:Rr;
    for tt=1:Cr; 
        if p(k,tt)~=0; 
            plot(tt,p(k,tt)-2*k,'k.','MarkerSize',15);
        end;
    end;
end;
% xlabel('Time')
% ylabel('Neuron #')
title({'Noise'})
axis('off')
axis([0 Cr -Rr 0]);
% 
% % 
% figure(2); hold on;
% subplot(2,2,4);
% hold on;
% stem(actual_potential_ratio,'bo','LineWidth',5);
% ylim([0 1])
% xlim([0 15]);
% xticks([1 5 10 14])
% xticklabels({'I','V','X','XIV'})
% xlabel('Motif-Class #')
% ylabel('Ratio Actual/Potential Contributions')
% grid on
% 
% set(0,'defaultAxesFontSize',30)
% 
% figure; hold on;
% hold on;
% bar(potential,'k')
% set(gca, 'YScale', 'log')
% xlim([0 15]);
% xticks([1 5 10 14])
% xticklabels({'I','V','X','XIV'})
% xlabel('Motif-Class #')
% ylabel('Number of Triplets (Scaled Count)')
% title('Potential Contributions')
% ylim([10e-2 10e4])
% set(gca,'YMinorTick','Off')
% 
% % 
figure; hold on;
hold on;
bar(actual,'r')
set(gca, 'YScale', 'log')
xlim([0 15]);
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
xlabel('Motif-Class #')
ylabel('Number of Triplets (Scaled Count)')
title('Actual Contributions')
ylim([10e-2 10e4])
set(gca,'YMinorTick','Off')
yticks([10e0 10e1 10e2 10e3 10e4 10e5])
% 
% 
% % motif_class_num = 1:1:14;
% figure; hold on;
% hold on;
% stem(actual_potential_ratio,'bo','LineWidth',11);
% ylim([0 0.6])
% xlim([0 15]);
% xticks([1 5 10 14])
% xticklabels({'I','V','X','XIV'})
% xlabel('Motif-Class #')
% ylabel('A/P')

figure; hold on;
hold on;
stem(actual_noise_ratio,'ro','LineWidth',11);
ylim([0 5])
xlim([0 15]);
xticks([1 5 10 14])
xticklabels({'I','V','X','XIV'})
xlabel('Motif-Class #')
ylabel('A/N')


%% Table 2 Data
LevelAct = sum(class_contribution([1 2 4]));
Sync=sum(class_contribution([4 5 7 8 12 13]));
LocDyn=sum(class_contribution([2 3 7 8 9 10 11]));
FeedF = sum(class_contribution([6 7 8 9 10 11 12 13 14]));
Div=sum(class_contribution([8 11 12 14]));
Conv=sum(class_contribution([7 9 13 14]));
%FeedF=sum(class_contribution([14]));
FeedB=sum(class_contribution([10]));
Tot=sum(class_contribution);
% Report %%
LevelAct=LevelAct*100/Tot
Sync=Sync*100/Tot
LocDyn=LocDyn*100/Tot
%NetProp = NetProp*100/Tot
Div=Div*100/Tot
Conv=Conv*100/Tot
FeedF=FeedF*100/Tot
FeedB=FeedB*100/Tot