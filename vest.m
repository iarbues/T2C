function S=vest(X,h);

[T,n]=size(X);

if h==0,
    S=std(X);
else
    
    X=X-repmat(mean(X),size(X,1),1);
    
    G=zeros(h+1,n);
    
    for k=0:h
        
        G(k+1,:)=sum(X((1+k):T,:).*X(1:(T-k),:))/(T-k);
                
    end
    
    % Uniform Kernel
    
    S=sum(G,1).^.5;
    
end

        
        