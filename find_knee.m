function [knee,num]=find_knee(pop,F,it)


Costs=[pop(F{1}).Cost];
n=numel(F{1});
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
% dis=d(so);
dis=d;
figure(1)
plot(Costs(1,so),dis(so)/5);
hold on

%% 寻找局部最大值
index=[];
for i=1:1:n
    tmpresult=[];
    tmpdis=[];
    for j=1:1:n
        newdis =  sum((Costs(:,i)-Costs(:,j)).^2);
        newdis = newdis^0.5;
        if newdis<1
            tmpresult=[tmpresult j];
            tmpdis=[tmpdis newdis];
        end
    end
    
    if max(dis(tmpresult))<=dis(i)
        if numel(index)==0
            index=i;
        else
            flag=0;
            for j=1:1:numel(index)
                knee_dis=sum((Costs(:,i)-Costs(:,index(j))).^2);
                knee_dis=knee_dis^0.5;
                if knee_dis<1
                    flag=1;
                    if dis(index(j))<dis(i)
                        index(j)=i;
                    end
                end
            end
            index=unique(index);
            if ~flag
                index=[index i];
            end
        end
    end
end

% [value index]=findpeaks(dis,'MinPeakDistance',2);
% index=find(x==index);
%     index=so(index);

knee=pop(index);
num = numel(index);