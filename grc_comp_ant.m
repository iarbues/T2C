function [MSEtmM, MSEtMm, MSEnaive, iopt, MSEi, meanA, meanB]=grc_comp(y,X,tau,sA,sB,lambda,im)

[T, k]=size(y);

R=floor(lambda*T);
P=T-R;

mA=length(sA);
mB=length(sB);

errorA=zeros(P,mA,k);
errorB=zeros(P,mB,k);

fact=sqrt(P-tau+1);
fact2=P-tau+1;

for i=1:mA
    
    iset=sA{i};
    
    y_tauR=y((1+tau):(R+tau),:);
    Xv=X(1:R,iset);
    M=(Xv'*Xv);
    b=Xv'*y_tauR;

    for t=(R+1):(T-tau)
        
        M=M+X(t,iset)'*X(t,iset);
        b=b+X(t,iset)'*y(t+tau,:);
        
        beta=M\b;
        
        errorA(t-R,i,:)=y(t+tau,:)-X(t,iset)*beta;

    end

end


for i=1:mB
    
    iset=sB{i};
    
    y_tauR=y((1+tau):(R+tau),:);
    Xv=X(1:R,iset);
    M=(Xv'*Xv);
    b=Xv'*y_tauR;

    for t=(R+1):(T-tau)
        
        M=M+X(t,iset)'*X(t,iset);
        b=b+X(t,iset)'*y(t+tau,:);
        
        beta=M\b;
        
        errorB(t-R,i,:)=y(t+tau,:)-X(t,iset)*beta;

    end

end

LA=errorA.^2;
LB=errorB.^2;

% TEST 1

tij=zeros(mA,mB,k);

for i=1:mA
    for j=1:mB
                       
        if im(i,j)
            
            u0=squeeze(errorA(:,i,:));
            uj=squeeze(errorB(:,j,:));
            c=u0.*(u0-uj);
                
            tij(i,j,:)=fact2*mean(c)./mean(uj.^2);
            
        else
            
            dif=squeeze(LA(:,i,:)-LB(:,j,:));
            S=vest(dif,tau-1);            
            tij(i,j,:)=fact*mean(dif)./S;
            
        end
        
    end
end

MSEi=max(tij,[],2);

MSEtmM=squeeze(min(MSEi,[],1));
MSEtMm=squeeze(max(min(tij,[],1),[],2));

MSEi=squeeze(MSEi);

% TEST 2

[qwe, iopt]=min(mean(LA,1),[],2);

for l=1:k
    maxtij(l)=max(tij(iopt(l),:,l));
end

MSEnaive=maxtij';


meanA=mean(LA);
meanB=mean(LB);


