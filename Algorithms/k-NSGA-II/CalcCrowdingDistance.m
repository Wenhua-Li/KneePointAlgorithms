%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA120
% Project Title: Non-dominated Sorting Genetic Algorithm II (NSGA-II)
% Publisher: Yarpiz (www.yarpiz.com)
%
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
%
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function CrowdDis=CalcCrowdingDistance(pop,F,kc,n_kc,rate)

nF=max(F);
cedu=0;

objs=pop.objs';
CrowdDis=zeros(size(F));
for k=1:nF
    index=find(F==k);
    
    Costs=objs(:,index);
    
    nObj=size(Costs,1);
    n=numel(index);
    d=zeros(n,nObj);
    
    for j=1:nObj
        [cj, so]=sort(Costs(j,:));
        d(so(1),j)=inf;
        cedu=cedu+(Costs(j,so(end))-Costs(j,so(1))).^2;
        for i=2:n-1
            d(so(i),j)=abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));
        end
        d(so(end),j)=inf;
    end
    for i=1:n  
        CrowdDis(index(i))=sum(d(i,:));
    end
    cedu=cedu^0.5;
    cedu=cedu*0.008*nObj;

    if n_kc~=0
        for i=1:n
            knee_center=kc.objs';
            dis = zeros(n_kc,1);
            for j=1:n_kc %给knee中心的赋值为最大，确保能保留knee
                dis(j)=sum((Costs(:,i)-knee_center(:,j)).^2);%+(Costs(2,i)-knee_center(2,j))^2;
                dis(j) = dis(j)^0.5;
            end
            d = min(dis);
            
            if d<cedu
                CrowdDis(index(i))=CrowdDis(index(i))+3;
            else
                CrowdDis(index(i))=CrowdDis(index(i))+1/(d+2);  
            end
        end
    end
    
end


end

