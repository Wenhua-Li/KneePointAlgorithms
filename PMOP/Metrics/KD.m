function Score = KDissimilarity_Final(PopObj,KneeSet)
% <metric> <min>
% The Knee-based dissimilarity between solutions and knee points.
    PopObj  = unique(PopObj,'rows');
    KneeSet = unique(KneeSet,'rows');
    [nr1,~] = size(KneeSet);
    KD = 0;
    for i=1:nr1
       dist = 0;
       dist = min(pdist2(KneeSet(i,:),PopObj),[],2);
       KD = KD+dist;
    end
    Score    = KD/nr1;
end