function [y, X]=simulaDGP2(T,b,tau,M)

n=8;

y=[randn(1,M); zeros(T-1,M)];

%C1=sqrtm(toeplitz(1:-0.01:0.97));
C1=sqrtm(toeplitz(1:-0.0025:0.9925));

X=nan(T,8);

X(:,2:5)=randn(T,4)*C1;

C2=sqrtm(toeplitz(1:-0.1:0.7));

X(:,[1 6:8])=randn(T,4)*C2;

u=randn(T,1);

g=zeros(1,7);

for t=1:T-tau

    y(t+tau,:)=0.3*y(t+tau-1,:)+b*X(t,1)+g*X(t,2:end)'+u(t+tau);

end
