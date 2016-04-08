
% this function performs the minimax and maximin tests. INPUT:
% y - predictand, X - all the regressors, tau - forecasting horizon, sA -
% cell array with all the models of first class, sB - cell array with all
% the models of the second class, lambda - parameter for the partition, im
% - matrix (i,j)=true when model i is nested in j, false otherwise, scheme
% - string ="recursive" for the recursive method, rolling otherwise.
%
% OUTPUT: MSEtmM, MSEtmM - test statistics, MSEnaive - obseolete, iopt -
% best OOS benchmark, meanA, meanB - vectors with MSFE of models. 

function [MSEtmM, MSEtMm, MSEnaive, iopt, MSEi, meanA, meanB]=grc_comp(y,X,tau,sA,sB,lambda,im,scheme)

[T, k]=size(y);

R=floor(lambda*T);
P=T-R-tau+1;

mA=length(sA);
mB=length(sB);

errorA=zeros(P,mA,k);
errorB=zeros(P,mB,k);

% factors for the statistics

fact2=P-tau+1;  
fact=sqrt(fact2);

for i=1:mA
    
    iset=sA{i};
    
    y_tauR=y((1+tau):(R-1),:); % first estimate
    Xv=X(1:(R-1-tau),iset);
    M=(Xv'*Xv);
    b=Xv'*y_tauR;
    
    for t=R:(T-tau)
        
        M=M+X(t-tau,iset)'*X(t-tau,iset);   % udpate design matrix
        b=b+X(t-tau,iset)'*y(t,:);          % update vector RHS
        
        beta=M\b;                                       % get new parameters
        errorA(t-R+1,i,:)=y(t+tau,:)-X(t,iset)*beta;    % get errors
        
        % if rolling, we remove values at the beginning of the window
        
        if ~isequal(scheme,'recursive') && t<T-tau      
                                                        
            M=M-X(t-R+1,iset)'*X(t-R+1,iset);
            b=b-X(t-R+1,iset)'*y(t-R+tau+1,:);
            
        end
                
    end
    
end

% same as above, for second class of models

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

LA=errorA.^2;       % MSFE
LB=errorB.^2;

tij=zeros(mA,mB,k);

for i=1:mA
    for j=1:mB
        
        % depending on whether nested or not
        
        if im(i,j)
            
            u0=squeeze(errorA(:,i,:));
            uj=squeeze(errorB(:,j,:));
            c=u0.*(u0-uj);
                
            tij(i,j,:)=fact2*mean(c)./mean(uj.^2); % ENC-F
                       
        else
            
            dif=squeeze(LA(:,i,:)-LB(:,j,:));
            S=vest(dif,tau-1);            
            tij(i,j,:)=fact*mean(dif)./S;           % MSE-T
            
        end
        
    end
end

[MSEi, modjM]=max(tij,[],2);

[aux, modi]=min(MSEi,[],1);

MSEtmM =squeeze(aux);

MSEtMm=squeeze(max(min(tij,[],1),[],2));

MSEi=squeeze(MSEi);

% TEST 2 (obsolete)

[~, iopt]=min(mean(LA,1),[],2);

%for l=1:k
%    maxtij(l)=max(tij(iopt(l),:,l));
%end

%MSEnaive=maxtij';

MSEnaive=nan;


meanA=mean(LA);
meanB=mean(LB);


