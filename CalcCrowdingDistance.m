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

function pop=CalcCrowdingDistance(pop,F,kc,n_kc,it)

nF=numel(F);

for k=1:nF
    
    Costs=[pop(F{k}).Cost];
    
    nObj=size(Costs,1);
    
    n=numel(F{k});
    
    d=zeros(n,nObj);
    
    for j=1:nObj
        
        [cj, so]=sort(Costs(j,:));
        
        d(so(1),j)=inf;
        
        for i=2:n-1
            
            d(so(i),j)=abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));
            
        end
        
        d(so(end),j)=inf;
        
    end
    
    caldis=0;
    for i=1:n
        
        pop(F{k}(i)).CrowdingDistance=sum(d(i,:));
        caldis=caldis+sum(d(i,:));
        
    end
    %caldis=n/caldis;
    
    cedu = sum((Costs(1,so(end))-Costs(1,so(1))).^2);%+(Costs(2,so(end))-Costs(2,so(1)))^2;
    cedu=cedu^0.5;
    cedu=cedu*0.02;
%     if it<
%     cedu=cedu*0.1;
%     elseif it<50
%         cedu=cedu*0.08;
%     elseif it<60
%             cedu=cedu*0.04;    
%     end
    if n_kc~=0
        for i=1:n
            knee_center=[kc.Cost];
            dis = zeros(n_kc,1);
            for j=1:n_kc %给knee中心的赋值为最大，确保能保留knee
                kc(j).CrowdingDistance = inf;
                dis(j)=sum((Costs(1,i)-knee_center(1,j)).^2);%+(Costs(2,i)-knee_center(2,j))^2;
                dis(j) = dis(j)^0.5;
                
            end
            d = min(dis);
            
            if d<cedu
                pop(F{k}(i)).CrowdingDistance=pop(F{k}(i)).CrowdingDistance+5;
            else
                pop(F{k}(i)).CrowdingDistance=pop(F{k}(i)).CrowdingDistance+1/(d+2);  
            end
%             pop(F{k}(i)).CrowdingDistance=0;
        end
    end
    
end


end

