x=[0:0.0001:1];
A=4;
B=1;
S=2;
u=x;
m1 = 1 + exp(cos(A.*power(u, B).*pi+pi./2))./(power(2, S).* A);
plot(x,m1);