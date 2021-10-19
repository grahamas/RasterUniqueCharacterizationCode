clear;
close all;

 
load ihot.mat                       % MUST BE AVAILABLE IN path or dir

 
C=150;  %300;  % columns
R=150;  %100;  % rows

 
% deterministic component
t=0:C-1;
s=1*sin(2*pi*.5*t); % used .04 and .08
S=zeros(R,C);
step=10;  % in paper Fig. 3, step=1 
for ch=1:step:R;
    ph=2*pi*(ch/R);
    % OR ... pick one pase shift    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SELECT
    ph=0;
    S(ch,:)=1*(((sin(2*pi*.1*t+ph)+1)/2));  % used .04 and .08
end;

 
% addition of a random component
x=rand(R,C);
% OR --- pick one of the x's        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SELECT
%x=zeros(R,C);
x=.51*x+.49*S;                        %%%%%%%%% change factors for x and S

 
% create raster
for n=1:R;
    for t=1:C;
        if x(n,t)>.5;       % was .5 for the figs in the draft
            p(n,t)=1;
        else;
            p(n,t)=0;
        end;
    end;
end;

 
% Only analyze part of the image
% determined by Rr and Cr
T1=15;  % was 20 for the draft
N1=15;   % was 5 for the draft
T2=15;  % was 20 for the draft
N2=15;   % was 5 for the draft

 
Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

 
NN=zeros(max(N1+1,N2+1));
TT=zeros(max(T1+1,T2+1));

 
%figure;
%pcolor(p(1:Rr,1:Cr));
%colormap(ihot);

 
figure;
axis([0 Cr 0 Rr]);
hold;
for k=1:Rr;
    for tt=1:Cr; 
        if p(k,tt)~=0; 
            plot(tt,p(k,tt)-2*k,'k.');
        end;
    end;
end;
xlabel('Time')
ylabel('Unit#')
title('Raster Plot of Multi-Unit Activity')
axis([0 Cr -Rr 0]);
axis('off');
count=0;
contributions=zeros(1,25);
nonzero=zeros(1,25);

 
% list of the lags
for t1=0:T1;
for n1=0:N1;
for t2=0:T2;
for n2=0:N2;

 
% compute Eq (15)
cd(n1+1,t1+1,n2+1,t2+1)=0;
count=count+1;
for c=1:Cr;
    for r=1:Rr;
        flag=0;
        addition=p(r,c)*p(r+n1,c+t1)*p(r+n2,c+t2);
        if addition ~=0; 
            flag=1; 
        end;
        cd(n1+1,t1+1,n2+1,t2+1)=cd(n1+1,t1+1,n2+1,t2+1)+addition;
        % count the scenario contributions
        if n1==0 & n2==0 & t1==0 & t2==0
            contributions(1)=contributions(1)+1;
            if flag==1;
            nonzero(1)=nonzero(1)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1==0 & n2==0 & t1==0 & t2~=0
            contributions(2)=contributions(2)+1;
            if flag==1;
            nonzero(2)=nonzero(2)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1==0 & n2==0 & t1~=0 & t2==0
            contributions(3)=contributions(3)+1;
            if flag==1;
            nonzero(3)=nonzero(3)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1==0 & n2==0 & t1~=0 & t2~=0 & t1~=t2
            contributions(4)=contributions(4)+1;
            if flag==1;
            nonzero(4)=nonzero(4)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1==0 & n2==0 & t1~=0 & t2~=0 & t1==t2
            contributions(5)=contributions(5)+1;
            if flag==1;
            nonzero(5)=nonzero(5)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1==0 & n2~=0 & t1==0 & t2==0 
            contributions(6)=contributions(6)+1;
            if flag==1;
            nonzero(6)=nonzero(6)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1==0 & n2~=0 & t1==0 & t2~=0 
            contributions(7)=contributions(7)+1;
            if flag==1;
            nonzero(7)=nonzero(7)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1==0 & n2~=0 & t1~=0 & t2==0 
            contributions(8)=contributions(8)+1;
            if flag==1;
            nonzero(8)=nonzero(8)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1==0 & n2~=0 & t1~=0 & t2~=0 & t1~=t2
            contributions(9)=contributions(9)+1;
            if flag==1;
            nonzero(9)=nonzero(9)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1==0 & n2~=0 & t1~=0 & t2~=0 & t1==t2
            contributions(10)=contributions(10)+1;
            if flag==1;
            nonzero(10)=nonzero(10)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2==0 & t1==0 & t2==0 
            contributions(11)=contributions(11)+1;
            if flag==1;
            nonzero(11)=nonzero(11)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2==0 & t1==0 & t2~=0 
            contributions(12)=contributions(12)+1;
            if flag==1;
            nonzero(12)=nonzero(12)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2==0 & t1~=0 & t2==0 
            contributions(13)=contributions(13)+1;
            if flag==1;
            nonzero(13)=nonzero(13)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2==0 & t1~=0 & t2~=0 & t1~=t2
            contributions(14)=contributions(14)+1;
            if flag==1;
            nonzero(14)=nonzero(14)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2==0 & t1~=0 & t2~=0 & t1==t2
            contributions(15)=contributions(15)+1;
            if flag==1;
            nonzero(15)=nonzero(15)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1==0 & t2==0 & n1~=n2
            contributions(16)=contributions(16)+1;
            if flag==1;
            nonzero(16)=nonzero(16)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1==0 & t2==0 & n1==n2
            contributions(17)=contributions(17)+1;
            if flag==1;
            nonzero(17)=nonzero(17)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1==0 & t2~=0 & n1~=n2
            contributions(18)=contributions(18)+1;
            if flag==1;
            nonzero(18)=nonzero(18)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1==0 & t2~=0 & n1==n2
            contributions(19)=contributions(19)+1;
            if flag==1;
            nonzero(19)=nonzero(19)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1~=0 & t2==0 & n1~=n2
            contributions(20)=contributions(20)+1;
            if flag==1;
            nonzero(20)=nonzero(20)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1~=0 & t2==0 & n1==n2
            contributions(21)=contributions(21)+1;
            if flag==1;
            nonzero(21)=nonzero(21)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1~=0 & t2~=0 & n1~=n2 & t1~=t2
            contributions(22)=contributions(22)+1;
            if flag==1;
            nonzero(22)=nonzero(22)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1~=0 & t2~=0 & n1~=n2 & t1==t2
            contributions(23)=contributions(23)+1;
            if flag==1;
            nonzero(23)=nonzero(23)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1~=0 & t2~=0 & n1==n2 & t1~=t2
            contributions(24)=contributions(24)+1;
            if flag==1;
            nonzero(24)=nonzero(24)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;
        if n1~=0 & n2~=0 & t1~=0 & t2~=0 & n1==n2 & t1==t2
            contributions(25)=contributions(25)+1;
            if flag==1;
            nonzero(25)=nonzero(25)+1;
            TT(t1+1,t2+1)=TT(t1+1,t2+1)+1;
            NN(n1+1,n2+1)=NN(n1+1,n2+1)+1;
            end;
        end;

    
    end;
end;
test(count)=cd(n1+1,t1+1,n2+1,t2+1);
end;
end;
end;
end;
% scale cd
Scale=1/((Rr-max(n1,n2))*(Cr-max(t1,t2)));
cds=cd.*Scale;
nonzeros=nonzero.*Scale;
tests=test.*Scale;

 
NNN=hist(tests);
figure;bar((100*NNN)./sum(NNN));
title('Distribution of Contributions to Triple Correlation')
xlabel('Contribution Category')
ylabel('%')
figure;hist(tests);
title('Distribution of Contributions to Triple Correlation')
xlabel('Contribution')
ylabel('n')
figure;bar(contributions)
title('All Potential Scenario Contributions')
xlabel('Scenario#')
ylabel('All Potential Contributions')
figure;bar(nonzeros)
title('Actual Scenario Contributions')
xlabel('Scenario#')
ylabel('Actual Contributions')
figure;stem(nonzeros./contributions,'o')
axis([0 26 0 2.5*10^(-5)]);
title('Ratio Actual/Potential Scenario Contributions')
xlabel('Scenario#')
ylabel('Ratio Actual/Potential Contributions')
figure;hold;
for tmp=1:25;
    tmpt=num2str(tmp);
    plot(contributions(tmp),nonzeros(tmp),'.');
    text(contributions(tmp),nonzeros(tmp),tmpt);
end;
title('Actual vs Potential Scenario Contributions by Scenario#')
xlabel('Potential Scenario Contributions')
ylabel('Actual Scenario Contributions')
% Report Results
Scale
max(max(max(max(cds))))
mean(mean(mean(mean(cds))))

 
%TEST
%figure;histogram2(TT3(:,1),TT3(:,2),100)
%figure;histogram2(TT9(:,1),TT9(:,2),100)
%figure;histogram2(TT14(:,1),TT14(:,2),100)
%figure;histogram2(TT22(:,1),TT22(:,2),100)
%figure;histogram2(TT24(:,1),TT24(:,2),100)

 
% Table 2 Data
Activity=sum(contributions([1]));
Sync=sum(contributions([6 8 10 11 15 16 17 18 19 20 21 23]));
LocDyn=sum(contributions([2 3 4 5 8 9 10 12 14 15 19 21 24]));
Div=sum(contributions([22 23]));
Conv=sum(contributions([18 20 22]));
FeedF=sum(contributions([22 23]));
FeedB=sum(contributions([9 14 24]));
Tot=sum(contributions);
% Report %%
Activity=Activity*100/Tot
Sync=Sync*100/Tot
LocDyn=LocDyn*100/Tot
Div=Div*100/Tot
Conv=Conv*100/Tot
FeedF=FeedF*100/Tot
FeedB=FeedB*100/Tot

 
figure;surf(TT)
figure;surf(NN)
comment=['Note that the numbers for t and n are augmented by 1 in the surf plots']