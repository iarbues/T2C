function [MSEtmM, MSEtMm, MSEnaive, iopt, MSEi, meanA, meanB]=grc(y,X,tau,sA,sB,lambda,scheme)

[T, k]=size(y);

R=floor(lambda*T);
P=T-R-tau+1;

mA=length(sA);
mB=length(sB);

errorA=zeros(P,mA,k);
errorB=zeros(P,mB,k);

fact=sqrt(P-tau+1);

for i=1:mA
    
    iset=sA{i};
    
    y_tauR=y((1+tau):(R-1),:);
    Xv=X(1:(R-1-tau),iset);
    M=(Xv'*Xv);
    b=Xv'*y_tauR;
    
    for t=R:(T-tau)
        
        M=M+X(t-tau,iset)'*X(t-tau,iset);
        b=b+X(t-tau,iset)'*y(t,:);
        
        beta=M\b;
        errorA(t-R+1,i,:)=y(t+tau,:)-X(t,iset)*beta;
        
        if ~isequal(scheme,'recursive') && t<T-tau
            
            M=M-X(t-R+1,iset)'*X(t-R+1,iset);
            b=b-X(t-R+1,iset)'*y(t-R+tau+1,:);
            
        end
                
    end
    
end


for i=1:mB
    
    iset=sB{i};
    
    y_tauR=y((1+tau):(R-1),:);
    Xv=X(1:(R-1-tau),iset);
    M=(Xv'*Xv);
    b=Xv'*y_tauR;

    for t=R:(T-tau)
        
        M=M+X(t-tau,iset)'*X(t-tau,iset);
        b=b+X(t-tau,iset)'*y(t,:);
        
        beta=M\b;
        errorB(t-R+1,i,:)=y(t+tau,:)-X(t,iset)*beta;
        
        if ~isequal(scheme,'recursive') && t<T-tau
            
            M=M-X(t-R+1,iset)'*X(t-R+1,iset);
            b=b-X(t-R+1,iset)'*y(t-R+tau+1,:);
            
        end

    end

end

LA=errorA.^2;
LB=errorB.^2;

% TEST 1

tij=zeros(mA,mB,k);

for i=1:mA
    for j=1:mB
        
        dif=squeeze(LA(:,i,:)-LB(:,j,:));
        
        S=vest(dif,tau-1);
        
        tij(i,j,:)=fact*mean(dif)./S;
        
    end
end

MSEi=max(tij,[],2);

MSEtmM=squeeze(min(MSEi,[],1));
MSEtMm=squeeze(max(min(tij,[],1),[],2));

MSEi=squeeze(MSEi);

% TEST 2

[~, iopt]=min(mean(LA,1),[],2);

for l=1:k
    maxtij(l)=max(tij(iopt(l),:,l));
end

MSEnaive=maxtij';


meanA=mean(LA);
meanB=mean(LB);


