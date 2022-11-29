function Score = KGD(PopObj,PF)
% <metric> <min>
% Knee-based Generational distance
    PF = unique(PF,'rows');
    Distance = min(pdist2(PopObj,PF),[],2);
    Score    = norm(Distance) / length(Distance);
end