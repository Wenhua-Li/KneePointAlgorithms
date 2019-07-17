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

function [pop, F]=SortPopulation(pop,kc,n_kc)

% Sort Based on Crowding Distance
[~, CDSO]=sort([pop.CrowdingDistance],'descend');
pop=pop(CDSO);

% Sort Based on Rank
[~, RSO]=sort([pop.Rank]);
pop=pop(RSO);

% Update Fronts
Ranks=[pop.Rank];
MaxRank=max(Ranks);
F=cell(MaxRank,1);
for r=1:MaxRank
    F{r}=find(Ranks==r);
end

oldpop = pop;

%%  加入knee region的保留策略

if n_kc~=0
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
        %         if numel(index)>n
        %             save_point=[save_point index(randperm(numel(index),n))];
        %         else
        %             save_point=[save_point index];
        %         end
        while numel(index)<n
            index=[index index];
        end
        save_point=[save_point index(randperm(numel(index),n))];
    end
    
    %     pop = [pop(so(save_point));pop];
    pop=[pop(so(save_point));oldpop];
    for i=1:1:numel(save_point)
        pop(i).Reserved=1;
    end
end


end