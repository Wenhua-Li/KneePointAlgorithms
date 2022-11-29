x=[0:0.0001:1];
A=4;
B=1;
S=2;
u=x;
m1 = 2 - exp(cos(A.*power(u, B).*pi)+0.5.*(cos(A.*power(u, B).*pi)-0.5).^4)./(power(2, S).* A);
plot(x,m1);
hold on
