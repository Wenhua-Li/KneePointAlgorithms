x=[0:0.0001:1];
A=6;
B=1;
S=-1;
u=x;
m1 = 2 + abs(sin(A.*power(u, B))-cos(A.*power(u, B)-pi./4))./(A.*power(2,S)); %% A.*power(2,S)
plot(x,m1)