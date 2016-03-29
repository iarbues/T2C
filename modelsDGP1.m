
% parameters of the simulation

P=floor(T*lambda);
maxlag=4;
intercept=false(1);

sA={};

for j=1:maxlag+1
    sA{j}=1:j;
end

sB={};

for j=1:maxlag+1
    for i=1:maxlag+1
        
        sB{end+1}=[1:i maxlag+1+(1:j)];       
        
    end
end

simDGP='simulaDGP1';
    