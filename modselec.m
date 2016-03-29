function modind=modselec(y,X,tau,sA,method)

[T, k]=size(y);

mA=length(sA);
nA=zeros(mA,1);
eA=zeros(mA,1);

for i=1:mA
    
    iset=sA{i};
    nA(i)=length(iset);
    
    y_tau=y((1+tau):T,:);
    Xv=X(1:(T-tau),iset);
    M=(Xv'*Xv);
    b=Xv'*y_tau;
    
    beta=M\b;
    
    res=y_tau-Xv*beta;
    
    eA(i)=mean(res.^2);    

end

if isequal(method,'AIC')

    [qwe, modind]=min(log(eA)+2*nA/T);
    
else
    
    [asd, modind]=min(log(eA)+nA*log(T)/T);
    
end
    

