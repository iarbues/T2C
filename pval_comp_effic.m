
% With this program, we do not nest the loops used to estimate the true
% distribution of the statistic and the Bootstrap distribution.

% parameters of the simulation

P=floor(T*lambda);

M=400; % number of simulations of the test and the bootstrap
MB=400;

cont=0;

MSEmM_pval=zeros(M,1);  % minimax
MSEMm_pval=zeros(M,1);  % maximin
MSEnaive_0_pval=zeros(M,1); % best OOS benchmark (not reported in article)
MSEnaive_1a_pval=zeros(M,1); % max_j MSE-t (i,j) with i <- AIC
MSEnaive_1b_pval=zeros(M,1); % max_j MSE-t (i,j) with i <- BIC
MSEnaive_2_pval=zeros(M,1); % largest benchmark

error_cont=0;

all_reg_uni=[];             % largest benchmark
for k=1:length(sA);
    all_reg_uni=union(all_reg_uni,sA{k});
end

% determine the index of the largest benchmark
for k=1:length(sA);
    if isempty(setxor(sA{k},all_reg_uni))
       largest_index=k; 
    end
end

nA=length(sA);
nB=length(sB);

im=true(nA,nB);

for i=1:nA              % nested/nonnested matrix
    for j=1:nB
        im(i,j)=isequal(intersect(sA{i},sB{j}),sA{i}); % subset
    end
end

scheme='rolling';

sample_test_mM=zeros(M,1);
sample_test_Mm=zeros(M,1);
sample_test_naive_0=zeros(M,1);
sample_test_naive_1a=zeros(M,1);
sample_test_naive_1b=zeros(M,1);
sample_test_naive_2=zeros(M,1);


cont_selec_naive_0=zeros(length(sA),1);
cont_selec_naive_1a=zeros(length(sA),1);
cont_selec_naive_1b=zeros(length(sA),1);
cont_selec_naive_2=zeros(length(sA),1);

for cont=1:M
    
    cont
    
    [y, X]=feval(simDGP,T,b,tau,1); % simulate data
    
    ylags=zeros(T-maxlag,maxlag+1);
    
    % lags of y
    
    for j=0:maxlag                  
        
        ylags(:,j+1)=y((1+maxlag-j):(T-j));
        
    end
    
    Xr=[ylags X(1+maxlag:end,:)]; % append the prediction variable and its lags to the regressors
    
    if intercept,
        Xr=[ones(size(Xr,1),1) Xr];     % if there is intercept, append it
    end
    
    yr=y(1+maxlag:end);    % truncation of first values
    
    % we do the test
    
    [MSEmM, MSEMm, MSEnaive, iopt, MSEi, meanA, meanB]=...
        grc_comp(yr,Xr,tau,sA,sB,lambda,im,scheme);
    
    sample_test_mM(cont)=MSEmM;
    sample_test_Mm(cont)=MSEMm;
    
    % Method 0
    
    modB=[];
    for j=1:length(sB)
        if isequal(intersect(sB{j},sA{iopt}),sA{iopt}),
            modB=[modB;j];
        end
    end
    
    
    [MSEnaive_0, ~, ~, ~, ~, meanA_0, meanB_0]=...
        grc_comp(y(1+maxlag:end),Xr,tau,sA(iopt),sB(modB),lambda,im(iopt,modB),scheme);
    

    sample_test_naive_0(cont)=MSEnaive_0;
    cont_selec_naive_0(iopt)=cont_selec_naive_0(iopt)+1;

    
    % Method 1a
    
    modA=modselec(yr(1:R),Xr(1:R,:),tau,sA,'AIC');
    cont_selec_naive_1a(modA)=cont_selec_naive_1a(modA)+1;

    
    modB=[];
    
    for j=1:nB
        if isequal(intersect(sB{j},sA{modA}),sA{modA}),
            modB=[modB;j];
        end
    end
    
    [MSEnaive_1a, ~, ~, ~, ~, meanA_0, meanB_0]=...
        grc_comp(yr,Xr,tau,sA(modA),sB(modB),lambda,im(modA,modB),scheme);
    
    sample_test_naive_1a(cont)=MSEnaive_1a;
    
    %  Method 1b
    
    modA=modselec(yr(1:R),Xr(1:R,:),tau,sA,'BIC');
    cont_selec_naive_1b(modA)=cont_selec_naive_1b(modA)+1;

    modB=[];
    
    for j=1:nB
        if isequal(intersect(sB{j},sA{modA}),sA{modA}),
            modB=[modB;j];
        end
    end
    
    [MSEnaive_1b, ~, ~, ~, ~, meanA_1b, meanB_1b]=...
        grc_comp(y(1+maxlag:end),Xr,tau,sA(modA),sB(modB),lambda,im(modA,modB),scheme);

    sample_test_naive_1b(cont)=MSEnaive_1b;
    
    % Method 2
    
    modB=[];
    for j=1:length(sB)
        if all(ismember(all_reg_uni,sB{j}))
            modB=[modB;j];
        end
    end
    
    [MSEnaive_2, ~, ~, ~, wer, meanA_2, meanB_2]=...
        grc_comp(y(1+maxlag:end),Xr,tau,{all_reg_uni},sB(modB),lambda,ones(1,length(modB)),scheme);

    sample_test_naive_2(cont)=MSEnaive_2;
    
end
    
% Bootstrap Mm and mM

[cMSEmM, cMSEMm, cMSEi, sample_MSEmM_boot, sample_MSEMm_boot, sample_MSEi_boot]=...
        critval_bootstrap_comp(yr,Xr,sA,sB,lambda,tau,0.9,MB,dist,im,scheme);
    
   
% Bootstrap single benchmark

samples_single_boot=zeros(MB,length(sA));

for i=1:length(sA)
    
    % calculate bootstrap samples of single-benchmark for each i
    
    modA=i;
    
    modB=[];
    
    for j=1:nB
        if isequal(intersect(sB{j},sA{modA}),sA{modA}),
            modB=[modB;j];
        end
    end
    
    [~, ~, ~, sample_single_boot_i]=...
        critval_bootstrap_comp(y(1+maxlag:end),Xr,sA(modA),sB(modB),lambda,tau,0.9,MB,dist,im(i,modB),scheme);
    
    samples_single_boot(:,i)=sample_single_boot_i;
    
end

% selection probabilities

prob_selec_naive_0=cont_selec_naive_0/M;
prob_selec_naive_1a=cont_selec_naive_1a/M;
prob_selec_naive_1b=cont_selec_naive_1b/M;
prob_selec_naive_2=zeros(length(sA),1);
prob_selec_naive_2(largest_index)=1;

% CURVES

sorted_mM=sort(sample_test_mM);
sorted_Mm=sort(sample_test_Mm);
sorted_naive_0=sort(sample_test_naive_0);
sorted_naive_1a=sort(sample_test_naive_1a);
sorted_naive_1b=sort(sample_test_naive_1b);
sorted_naive_2=sort(sample_test_naive_2);
    
p_gen=(0:M+1)'/(M+1);
q_test_mM=zeros(M+2,1);
q_test_Mm=zeros(M+2,1);
q_naive_0=zeros(M+2,1);
q_naive_1a=zeros(M+2,1);
q_naive_1b=zeros(M+2,1);
q_naive_2=zeros(M+2,1);


for k=1:M
    
    % curves of minimax and maximin
    
    q_test_mM(k+1)=mean(sample_MSEmM_boot<sorted_mM(k));
    
    q_test_Mm(k+1)=mean(sample_MSEMm_boot<sorted_Mm(k));
        
    % curves of naïve methods
    
    q_naive_0(k+1)=mean(samples_single_boot<sorted_naive_0(k))*prob_selec_naive_0;
    
    q_naive_1a(k+1)=mean(samples_single_boot<sorted_naive_1a(k))*prob_selec_naive_0;
    
    q_naive_1b(k+1)=mean(samples_single_boot<sorted_naive_1b(k))*prob_selec_naive_1b;
    
    q_naive_2(k+1)=mean(samples_single_boot<sorted_naive_2(k))*prob_selec_naive_2;
           
end

q_test_mM(M+2)=1;
q_test_Mm(M+2)=1;
q_naive_0(M+2)=1;
q_naive_1a(M+2)=1;
q_naive_1b(M+2)=1;
q_naive_2(M+2)=1;

% simulate pvalues to fit with the rest of the code


pd_mM = makedist('PiecewiseLinear','x',p_gen,'Fx',q_test_mM);
pd_Mm = makedist('PiecewiseLinear','x',p_gen,'Fx',q_test_Mm);
pd_naive_0 = makedist('PiecewiseLinear','x',p_gen,'Fx',q_naive_0);
pd_naive_a1 = makedist('PiecewiseLinear','x',p_gen,'Fx',q_naive_1a);
pd_naive_1b = makedist('PiecewiseLinear','x',p_gen,'Fx',q_naive_1b);
pd_naive_2 = makedist('PiecewiseLinear','x',p_gen,'Fx',q_naive_2);

MSEmM_pval = random(pd_mM,M,1);
MSEMm_pval = random(pd_Mm,M,1);
MSEnaive_0_pval = random(pd_naive_0,M,1);
MSEnaive_1a_pval = random(pd_naive_a1,M,1);
MSEnaive_1b_pval = random(pd_naive_1b,M,1);
MSEnaive_2_pval = random(pd_naive_2,M,1);

       