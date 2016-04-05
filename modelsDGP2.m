
% parameters of the simulation

P=floor(T*lambda);
maxlag=0;
intercept=true(1);

nx=8;

sA={};
sB={};

for j=0:(2^nx-1)
    
    s=[];
    
    for k=1:nx
        if bitget(uint8(j),k),
            s=[s k];
        end
    end
    
    if any(ismember([1 6 7 8],s))
        sB{end+1}=[1 2 2+s];       % intercept, lag of y, x1,x2, ...
    else        
        sA{end+1}=[1 2 2+s];       
    end
        
    
end
    
simDGP='simulaDGP2';
    