
function [cMSEmM, cMSEMm, cMSEi, varargout]=critval_bootstrap_comp(y,X,sA,sB,lambda,tau,alpha,M,dist,im,scheme)

if dist==0,
    p=(sqrt(5)+1)/(2*sqrt(5));
    x1=-(sqrt(5)-1)/2;
    x2=(sqrt(5)+1)/2;
else
    p=0.5;
    x1=-1;
    x2=1;
end

dif=x2-x1;

T=size(y,1);
R=floor(lambda*T);
P=T-R;

k=400;
m=floor(M/k);

nA=length(sA);
nB=length(sB);

a=[];
for i=1:nA
    a=union(a,sA{i});
end

b=[];
for i=1:nB
    b=union(b,sB{i});
end

model.ar={};
model.ma={1:(tau-1)};
model.dif=[];
model.nreg=0;

% estimate with all regressors

y_tau=y((1+tau):T);
X_tau=X(1:(T-tau),b);

beta=(X_tau(1:R,:)'*X_tau(1:R,:))\X_tau(1:R,:)'*y_tau(1:R);

v_tau=y_tau-X_tau*beta;

MSEmM_=zeros(m*k,1);
MSEMm_=zeros(m*k,1);
MSEi_=zeros(m*k,nA);

Xa0=X(1:(R-tau),a);
y0=y((1+tau):R);
Xa=X(1:(T-tau),a);

beta0=(Xa0'*Xa0)\Xa0'*y0;

for j=1:m
    
    j;
    
    if tau>1,

        [theta, res]=clsarmaest(v_tau,model);

        u=res;

        eta=x1+(rand(T-tau,k)>p)*dif;

        v_=eta.*repmat(v_tau,1,k);

        u_=v_;

        for h=1:tau-1

            u_=u_-theta(h)*[zeros(h,k); v_(1:(end-h),:)];

        end

    else

        eta=x1+(rand(T-tau,k)>p)*dif; % Two-point distribution
        %eta=randn(T-tau,k);     % Normal distribution
                
        u_=eta.*repmat(abs(v_tau),1,k);

    end

    X_=X((1+tau):end,:);
    
    y_=repmat(Xa*beta0,1,k)+u_; % Wild Bootstrap CMC
    %y_=repmat(Xa*beta0,1,k)+randn(size(Xa,1),k)*mean(mean(u_));     
 
    [MSEmM, MSEMm, qwe, asd, MSEi , zxc, meanB]=grc_comp(y_,X_,tau,sA,sB,lambda,im,scheme);

    MSEmM_(((j-1)*k+1):(j*k))=MSEmM;
    MSEMm_(((j-1)*k+1):(j*k))=MSEMm;
    MSEi_(((j-1)*k+1):(j*k),:)=MSEi';

end

cMSEmM=quantile(MSEmM_,alpha);
cMSEMm=quantile(MSEMm_,alpha);
cMSEi=quantile(MSEi_,alpha);

if nargout>3,
    varargout{1}=MSEmM_;
    varargout{2}=MSEMm_;
    varargout{3}=MSEi_;
end
