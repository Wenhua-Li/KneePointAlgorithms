

Costs=[PF.Cost];
n=50;
[cj, so]=sort(Costs(1,:));
%% 确定窗口大小
WinSize=(Costs(so(end))-Costs(so(1)));
WinSize=WinSize/10;

p1=Costs(:,so(1));p2=Costs(:,so(end));
A=p2(2)-p1(2);B=p2(1)-p1(1);
C=p2(1)*p1(2)-p1(1)*p2(2);

d=zeros(1,n);
for i=1:1:n
    d(i)=A*Costs(1,i)-B*Costs(2,i)+C;
end
% d=-d;
dis=d(so);

%% 寻找局部最大值
index=[];
for i=1:1:n
%     currpoint=Costs(:,i);
tmpresult=[];
    for j=1:1:n

%         if i~=j
           newdis =  (Costs(1,i)-Costs(1,j))^2+(Costs(2,i)-Costs(2,j))^2;
           newdis = newdis^0.5;
           if newdis<2
               tmpresult=[tmpresult j];
           end
%         end
    end
    if max(dis(tmpresult))==dis(i)
        index=[index i];
    end
end


% [value index]=findpeaks(dis,'MinPeakDistance',5);
index=so(index);

% index=[];
% for i=1:1:numel(value)
%     index(i) = find(d==value(i));
% end

knee=pop(index);
num = numel(index);