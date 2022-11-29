x=[0:0.0001:1];
A=4;
B=1;
S=2;
u=x;
m1 = 1 + exp(sin(A.*power(u, B).*pi+pi./2))./(power(2, S).* A);


f1=m1.*(1-cos(u.*pi./2));
f2=m1.*(1-sin(u.*pi./2));

% plot(f1,f2)
hold on;
plot(x,m1);
hold on
