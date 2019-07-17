function z=DO2DK(x)
n=numel(x);
K=3;
s=1;
g=1+9*sum(x(2:end))/(n-1);
r=5+10*(x(1)-0.5)^2+cos(2*K*pi*x(1))*2^(s/2)/K;
f1=g*r*(sin(pi*x(1)/2^(s+1)+(1+(2^s-1)/2^(s+2))*pi)+1);
f2=g*r*(cos(pi*x(1)/2+pi)+1);

z=[f1;f2];