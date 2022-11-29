x=[0:0.0001:1];
A=1;
B=1;
S=2;
l =12;
u=x;
m1 = 2 + min(sin(2.*A.*pi.*power(u,B)), cos(2.*A.*pi.*power(u,B)-pi./l))./(power(2,S)); %% A.*power(2,S)
plot(x,m1);
