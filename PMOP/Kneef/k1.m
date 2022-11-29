x=[0:0.0001:1];
A=2;
B=1;
S=-2;
u=x;
m1 = 5+10.*(u-0.5).^2+ cos(A.*u.^B.*pi)./(A.*2.^S);
plot(x,m1)
[row,col] = size(x);
%% the minima
index = round(col/2)
[mm,mxid] = min(m1(1,1:index));
[mm2,mxid2] = min(m1(1,index:col));
%% corresponding x 
xknee = [];
minx = x(mxid);
minx2 = x(mxid2+index-1);
xknee =[xknee;minx;minx2];