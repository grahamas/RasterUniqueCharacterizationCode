#using CairoMakie

C=1000;  # columns
R=1000;  # rows

 
# deterministic component
t = 0:C-1;
s = @. 1*sin(2*pi*.04*t);
S = zeros(R,C);
for ch=1:R;
    ph=2*pi*(ch/R);
    # OR ... pick one pase shift    ############################## SELECT
    ph=0;
    @. S[ch,:] = 1*(((sin(2*pi*.04*t+ph)+1)/2));
end;

 
# addition of a random component
x=rand(R,C);
# OR --- pick one of the x"s        ############################## SELECT
#x=zeros(R,C);
x=0*x+1*S;                        ######### change factors for x and S

p = Matrix{eltype(x)}(undef, size(x)...) 

# create raster
for n=1:R;
    for t=1:C;
        if x[n,t]>.5;       # was .5 for the figs in the draft
            p[n,t]=1;
        else;
            p[n,t]=0;
        end;
    end;
end;

 
# Only analyze part of the image
# determined by Rr and Cr
T1=999;  # was 20 for the draft
N1=999;   # was 5 for the draft
T2=999;  # was 20 for the draft
N2=999;   # was 5 for the draft


# raster_fig = Figure()
# raster_ax = Axis(raster_fig, title = "Raster Plot of Multi-Unit Activity")
# for k=1:Rr;
#     for tt=1:Cr; 
#         if p[k,tt] != 0; 
#             lines!(ax, tt,p[k,tt]-2*k,color=:black);
#         end;
#     end;
# end;
# xlabel!(ax, "Time")
# ylabel!(ax, "Unit#")
# # axis([0 Cr -Rr 0]);
# # axis("off");
count=0;
nonzero=zeros(25);


function compute_potential_contributions(T1,N1,T2,N2)

contributions=zeros(25);
Cr=C-max(T1,T2);
Rr=R-max(N1,N2);

# list of the lags
for t1=0:T1, n1=0:N1, t2=0:T2, n2=0:N2
    N_possible_base_nodes = Cr * Rr#max(0, C - t1) * max(0, C - t2) * max(0, R - n1) * max(0, R - n2)
 
# compute Eq (15)
#cd[n1+1,t1+1,n2+1,t2+1]=0;
#count=count+1;
    # addition=p(r,c)*p(r+n1,c+t1)*p(r+n2,c+t2);
    # if addition  != 0; 
    #     flag=1; 
    # end;
    #cd[n1+1,t1+1,n2+1,t2+1]=cd[n1+1,t1+1,n2+1,t2+1]+addition;
    # count the scenario contributions
    if ((n1==0) && (n2==0) && t1==0 && t2==0)
        contributions[1] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[1]=nonzero[1] + N_possible_base_nodes;
        #end;
    elseif (n1==0 && n2==0 && t1==0 && t2 != 0)
        contributions[2] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[2]=nonzero[2] + N_possible_base_nodes;
        #end;
    elseif (n1==0 && n2==0 && t1 != 0 && t2==0)
        contributions[3] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[3]=nonzero[3] + N_possible_base_nodes;
        #end;
    elseif (n1==0 && n2==0 && t1 != 0 && t2 != 0 && t1 != t2)
        contributions[4] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[4]=nonzero[4] + N_possible_base_nodes;
        #end;
    elseif (n1==0 && n2==0 && t1 != 0 && t2 != 0 && t1==t2)
        contributions[5] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[5]=nonzero[5] + N_possible_base_nodes;
        #end;
    elseif (n1==0 && n2 != 0 && t1==0 && t2==0 )
        contributions[6] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[6]=nonzero[6] + N_possible_base_nodes;
        #end;
    elseif (n1==0 && n2 != 0 && t1==0 && t2 != 0 )
        contributions[7] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[7]=nonzero[7] + N_possible_base_nodes;
        #end;
    elseif (n1==0 && n2 != 0 && t1 != 0 && t2==0 )
        contributions[8] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[8]=nonzero[8] + N_possible_base_nodes;
        #end;
    elseif (n1==0 && n2 != 0 && t1 != 0 && t2 != 0 && t1 != t2)
        contributions[9] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[9]=nonzero[9] + N_possible_base_nodes;
        ##TT9(nonzero[9],:)=[t1 t2];
        #end;
    elseif (n1==0 && n2 != 0 && t1 != 0 && t2 != 0 && t1==t2)
        contributions[10] += N_possible_base_nodes;
        #if flag==1;
        #nonzero(10)=nonzero(10)+1;
        #end;
    elseif (n1 != 0 && n2==0 && t1==0 && t2==0 )
        contributions[11] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[11]=nonzero[11] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2==0 && t1==0 && t2 != 0 )
        contributions[12] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[12]=nonzero[12] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2==0 && t1 != 0 && t2==0 )
        contributions[13] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[13]=nonzero[13] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2==0 && t1 != 0 && t2 != 0 && t1 != t2)
        contributions[14] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[14]=nonzero[14] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2==0 && t1 != 0 && t2 != 0 && t1==t2)
        contributions[15] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[15]=nonzero[15] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1==0 && t2==0 && n1 != n2)
        contributions[16] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[16]=nonzero[16] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1==0 && t2==0 && n1==n2)
        contributions[17] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[17]=nonzero[17] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1==0 && t2 != 0 && n1 != n2)
        contributions[18] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[18]=nonzero[18] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1==0 && t2 != 0 && n1==n2)
        contributions[19] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[19]=nonzero[19] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1 != 0 && t2==0 && n1 != n2)
        contributions[20] += N_possible_base_nodes;
        #if flag==1;
        #nonzero(20)=nonzero(20)+1;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1 != 0 && t2==0 && n1==n2)
        contributions[21] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[21]=nonzero[21] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1 != 0 && t2 != 0 && n1 != n2 && t1 != t2)
        contributions[22] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[22]=nonzero[22] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1 != 0 && t2 != 0 && n1 != n2 && t1==t2)
        contributions[23] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[23]=nonzero[23] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1 != 0 && t2 != 0 && n1==n2 && t1 != t2)
        contributions[24] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[24]=nonzero[24] + N_possible_base_nodes;
        #end;
    elseif (n1 != 0 && n2 != 0 && t1 != 0 && t2 != 0 && n1==n2 && t1==t2)
        contributions[25] += N_possible_base_nodes;
        #if flag==1;
        #nonzero[25]=nonzero[25] + N_possible_base_nodes;
        #end;
    end;    
end;
return contributions
end;

potential_contributions = compute_potential_contributions(T1, N1, T2, N2)
# scale cd
Scale=1/((Rr-max(N1,N2))*(Cr-max(T1,T2)));
#cds=cd.*Scale;
nonzeros=nonzero.*Scale;
#tests=test.*Scale;

potential_contributions
 
# NN=hist(tests);
# figure;bar((100*NN)./sum(NN));
# title("Distribution of Contributions to Triple Correlation")
# xlabel("Contribution Category")
# ylabel("#")
# figure;hist(tests);
# title("Distribution of Contributions to Triple Correlation")
# xlabel("Contribution")
# ylabel("n")
# figure;bar(contributions)
# title("All Potential Scenario Contributions")
# xlabel("Scenario#")
# ylabel("All Potential Contributions")
# figure;bar(nonzeros)
# title("Actual Scenario Contributions")
# xlabel("Scenario#")
# ylabel("Actual Contributions")
# figure;stem(nonzeros./contributions,"o")
# axis([0 26 0 2.5*10^(-5)]);
# title("Ratio Actual/Potential Scenario Contributions")
# xlabel("Scenario#")
# ylabel("Ratio Actual/Potential Contributions")
# figure;hold;
# for tmp=1:25;
#     tmpt=num2str(tmp);
#     plot(contributions(tmp),nonzeros(tmp),".");
#     text(contributions(tmp),nonzeros(tmp),tmpt);
# end;
# title("Actual vs Potential Scenario Contributions by Scenario#")
# xlabel("Potential Scenario Contributions")
# ylabel("Actual Scenario Contributions")
# # Report Results
# Scale
# max(max(max(max(cds))))
# mean(mean(mean(mean(cds))))

 
#TEST
#figure;histogram2(TT3(:,1),TT3(:,2),100)
#figure;histogram2(TT9(:,1),TT9(:,2),100)
#figure;histogram2(TT14(:,1),TT14(:,2),100)
#figure;histogram2(TT22(:,1),TT22(:,2),100)
#figure;histogram2(TT24(:,1),TT24(:,2),100)

 
# Table 2 Data
Sync=sum(potential_contributions[[6 8 10 11 15 16 17 18 19 20 21 23]]);
LocDyn=sum(potential_contributions[[2 3 4 5 8 9 10 12 14 15 19 21 24]]);
Div=sum(potential_contributions[[22 23]]);
Conv=sum(potential_contributions[[18 20 22]]);
FeedF=sum(potential_contributions[[22 23]]);
FeedB=sum(potential_contributions[[9 14 24]]);
Tot=sum(potential_contributions);
# Report ##
Sync=Sync*100/Tot; @show Sync
LocDyn=LocDyn*100/Tot; @show LocDyn
Div=Div*100/Tot; @show Div
Conv=Conv*100/Tot; @show Conv
FeedF=FeedF*100/Tot; @show FeedF
FeedB=FeedB*100/Tot; @show FeedB

for (i, contribution) in enumerate(potential_contributions)
    print("Scenario $i:  $(contribution / Tot)\n")
end