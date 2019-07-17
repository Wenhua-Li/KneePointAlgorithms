function z=DEB2DK(x)

n=numel(x);
K=4;
% s=1;
% g=1+(9/(n-1)) * sum(x(2:end));
% r=5+10*(x(1)-0.5)^2+(1/K) * cos(2*K*pi*x(1));
% f1=g*r*(sin(pi*x(1)/2));
% f2=g*r*(cos(pi*x(1)/2));
% z=[f1 f2]';

g=1+(9/(n-1)) * sum(x(2:end));
r=5+10*(x(1)-0.5)^2+(1/K) * cos(2*K*pi*x(1));
f1=g*r*(sin(pi*x(1)/2));
f2=g*r*(cos(pi*x(1)/2));
z=[f1 f2]';