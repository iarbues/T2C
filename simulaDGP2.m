function [y, X]=simulaDGP2(T,b,tau,M)

n=8;

y=[randn(1,M); zeros(T-1,M)];

C=sqrtm(toeplitz(1:-0.05:0.65));

X=randn(T,n)*C;

u=randn(T,1);

g=zeros(1,7);

for t=1:T-tau

    y(t+tau,:)=0.3*y(t+tau-1,:)+b*X(t,1)+g*X(t,2:end)'+u(t+tau);

end
