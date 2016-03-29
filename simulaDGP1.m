function [y, X]=simulaDGP1(T,b,tau,M)

phi1=[0.6 -0.3 0.2];
phi2=[0.4 0.2 0.1];

r=6;
p=3;

y=randn(T+r,1);
x=randn(T+r,1);

for t=p+1:T+r
    
    y(t)=phi1(end:-1:1)*y((t-p):(t-1))+b*x(t-1)+randn;
    x(t)=phi2(end:-1:1)*x((t-p):(t-1))+randn;
    
end

X=[];

for i=0:(r-1)
    
    X=[X x((r+1-i):(end-i))];
    
end

y=y(r+1:end);
