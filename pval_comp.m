
% parameters of the simulation

P=floor(T*lambda);

M=400; % number of simulations of the test
MB=5000; % number of bootstrap simulations for each test

cont=0;

MSEmM_pval=zeros(M,1);  % minimax
MSEMm_pval=zeros(M,1);  % maximin
MSEnaive_0_pval=zeros(M,1); % best OOS benchmark (not reported in article)
MSEnaive_1a_pval=zeros(M,1); % max_j MSE-t (i,j) with i <- AIC
MSEnaive_1b_pval=zeros(M,1); % max_j MSE-t (i,j) with i <- BIC
MSEnaive_2_pval=zeros(M,1); % largest benchmark

error_cont=0;

all_reg_uni=[];             % greatest benchmark
for k=1:length(sA);
    all_reg_uni=union(all_reg_uni,sA{k});
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
    
    % get p-values or critical values
    
    [cMSEmM, cMSEMm, cMSEi, MSEmM_boot_k, MSEMm_boot_k, MSEi_boot_k]=...
        critval_bootstrap_comp(yr,Xr,sA,sB,lambda,tau,0.9,MB,dist,im,scheme);
    
    MSEmM_pval(cont)=mean(MSEmM_boot_k>MSEmM);
    MSEMm_pval(cont)=mean(MSEMm_boot_k>MSEMm);
    
    % Method 0
    
    modB=[];
    for j=1:length(sB)
        if isequal(intersect(sB{j},sA{iopt}),sA{iopt}),
            modB=[modB;j];
        end
    end
    
    [MSEnaive_0, ~, ~, ~, ~, meanA_0, meanB_0]=...
        grc_comp(y(1+maxlag:end),Xr,tau,sA(iopt),sB(modB),lambda,im(iopt,modB),scheme);
    [~, ~, ~, MSE_naive0_boot_k]=...
        critval_bootstrap_comp(y(1+maxlag:end),Xr,sA(iopt),sB(modB),lambda,tau,0.9,MB,dist,im(iopt,modB),scheme);
    
    MSEnaive_0_pval(cont)=mean(MSE_naive0_boot_k>MSEnaive_0);
    
    % Method 1a
    
    modA=modselec(yr(1:R),Xr(1:R,:),tau,sA,'AIC');
    
    modB=[];
    
    for j=1:nB
        if isequal(intersect(sB{j},sA{modA}),sA{modA}),
            modB=[modB;j];
        end
    end
    
    [MSEnaive_1a, ~, ~, ~, ~, meanA_0, meanB_0]=...
        grc_comp(yr,Xr,tau,sA(modA),sB(modB),lambda,im(modA,modB),scheme);
    [~, ~, ~, MSEnaive_boot_1a]=...
        critval_bootstrap_comp(yr,Xr,sA(modA),sB(modB),lambda,tau,0.9,MB,dist,im(modA,modB),scheme);
    
    MSEnaive_1a_pval(cont)=mean(MSEnaive_boot_1a>MSEnaive_1a);
    
    %  Method 1b
    
    modA=modselec(yr(1:R),Xr(1:R,:),tau,sA,'BIC');
    
    modB=[];
    
    for j=1:nB
        if isequal(intersect(sB{j},sA{modA}),sA{modA}),
            modB=[modB;j];
        end
    end
    
    [MSEnaive_1b, ~, ~, ~, ~, meanA_1b, meanB_1b]=...
        grc_comp(y(1+maxlag:end),Xr,tau,sA(modA),sB(modB),lambda,im(modA,modB),scheme);
    [~, ~, ~, MSEnaive_boot_1b]=...
        critval_bootstrap_comp(y(1+maxlag:end),Xr,sA(modA),sB(modB),lambda,tau,0.9,MB,dist,im(modA,modB),scheme);
    
    MSEnaive_1b_pval(cont)=mean(MSEnaive_boot_1b>MSEnaive_1b);
    
    % Method 2
    
    modB=[];
    for j=1:length(sB)
        if all(ismember(all_reg_uni,sB{j}))
            modB=[modB;j];
        end
    end
    
    [MSEnaive_2, ~, ~, ~, wer, meanA_2, meanB_2]=...
        grc_comp(y(1+maxlag:end),Xr,tau,{all_reg_uni},sB(modB),lambda,ones(1,length(modB)),scheme);
    [qwe, asd, zxc, MSEnaive_boot_2]=...
        critval_bootstrap_comp(y(1+maxlag:end),Xr,{all_reg_uni},sB(modB),lambda,tau,0.9,MB,dist,ones(1,length(modB)),scheme);
    
    MSEnaive_2_pval(cont)=mean(MSEnaive_boot_2>MSEnaive_2);
    
    
    
end
