function Score = IGD(PopObj,PF)
% <metric> <min>
% Knee-based Inverted generational distance
    PF = unique(PF,'rows');
    Distance = min(pdist2(PF,PopObj),[],2);
    Score    = mean(Distance);
end