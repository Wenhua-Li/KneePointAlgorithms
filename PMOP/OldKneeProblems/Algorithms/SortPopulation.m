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

function [newpop, F,CrowdDis]=SortPopulation(pop,kc,n_kc,FrontNo,CrowdDis,N)

% Sort Based on Crowding Distance
[~, CDSO]=sort(CrowdDis,'descend');
pop=pop(CDSO);
FrontNo=FrontNo(CDSO);
CrowdDis=CrowdDis(CDSO);

% Sort Based on Rank
[~, RSO]=sort(FrontNo);
pop=pop(RSO);
FrontNo=FrontNo(RSO);
CrowdDis=CrowdDis(RSO);

%%  加入knee region的保留策略

if 0%n_kc~=0
    pop=pop(F{1});
    Costs=[pop.Cost];
    [cj, so]=sort(Costs(1,:));
    length = Costs(1,so(end))-Costs(1,so(1)); %目标跨度
    length = length*8/100;
    f1=Costs(1,so);
    n=size(Costs,1);
    n=round(0.12*n);
    save_point = [];
    for i=1:1:n_kc
        kneePoint = kc(i).Cost;
        [~,index]=find(abs(f1-kneePoint(1))<length);
        while numel(index)<n
            index=[index index];
        end
        save_point=[save_point index(randperm(numel(index),n))];
    end
    pop=[pop(so(save_point));oldpop];
end

newpop=pop(1:N);
F=FrontNo(1:N);
CrowdDis=CrowdDis(1:N);

end