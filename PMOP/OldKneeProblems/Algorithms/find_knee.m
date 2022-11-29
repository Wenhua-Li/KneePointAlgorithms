function [knee,num]=find_knee(pop,F,rate,kind)

% index=find(F==1);
% Costs=pop.objs';
n=numel(F);
Current=1:n;
PopObj=pop.objs;
M=size(PopObj,2);

%% 确定窗口大小

% Find the extreme points
[~,Rank]   = sort(PopObj(Current,:),'descend');
Extreme    = zeros(1,M);
Extreme(1) = Rank(1,1);
for j = 2 : length(Extreme)
    k = 1;
    Extreme(j) = Rank(k,j);
    while ismember(Extreme(j),Extreme(1:j-1))
        k = k+1;
        Extreme(j) = Rank(k,j);
    end
end
EdgePoints = sort(PopObj(Current(Extreme),:));
WinSize = sum((EdgePoints(end,:)-EdgePoints(1,:)).^2);
WinSize = 0.04*WinSize.^0.5*M;
% WinSize = 1.1*(1-exp(-M*rate));
% WinSize =1;

switch kind
    case 1
        % 计算hyperplane
        Hyperplane = PopObj(Current(Extreme),:)\ones(length(Extreme),1);
        % 计算距离
        Distance(Current) = -(PopObj(Current,:)*Hyperplane-1)./sqrt(sum(Hyperplane.^2));
        dis=Distance;
    case 2
        T=zeros(n);
        minobj=min(PopObj);
        maxobj=max(PopObj);
        for i=1:n
            for j=1:n
                tmp1=(PopObj(j,:)-PopObj(i,:))./(maxobj-minobj);
                tmp2=(PopObj(i,:)-PopObj(j,:))./(maxobj-minobj);
                T(i,j)=sum(max(0,tmp1))/sum(max(0,tmp2));
            end
        end
        dis=min(T');
    case 3
        [~,so]=sort(PopObj(:,1));
        f=PopObj(so,:);
        lambda=zeros(n);
        for i=1:n
            for j=1:n
                lambda(i,j)=(f(j,2)-f(i,2))/(f(i,1)-f(j,1)+ f(j,2)-f(i,2));
            end
        end
        
        U=zeros(1,n);
        for i=2:n-1
            tmp=(f(i,:)-f(i-1,:)).*[lambda(i-1,i) 1-lambda(i-1,i); lambda(i-1,i+1) 1-lambda(i-1,i+1)];
            tmp=sum(sum(tmp));
            tmp2=(f(i,:)-f(i-1,:)).*[lambda(i-1,i+1) 1-lambda(i-1,i+1);lambda(i,i+1) 1-lambda(i,i+1)];
            tmp2=sum(sum(tmp2));
            U(i)=tmp+tmp2;
        end
        dis=-U;
end

%% 寻找局部最大值
index=[];
for i=1:1:n
    tmpresult=[];
    tmpdis=[];
    for j=1:1:n
        newdis =  sum((PopObj(i,:)-PopObj(j,:)).^2);
        newdis = newdis^0.5;
        if newdis<WinSize
            tmpresult=[tmpresult j];
            tmpdis=[tmpdis newdis];
        end
    end
    
    if max(dis(tmpresult))<=dis(i)
        if numel(index)==0
            index=i;
        else
            index=[index i];
        end
    end
end

knee=pop(index);
num = numel(index);