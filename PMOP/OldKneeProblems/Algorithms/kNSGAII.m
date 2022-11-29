function kNSGAII(Global)
% <algorithm> <H-N>
% kind --- 1 --- The type of knee detect function
% yita --- 0.4 --- knee detect function stage
% Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
%% Parameter setting
[kind,yita] = Global.ParameterSet(1,0.4);
% kind=2;

%% init
kc=[];n_kc=0;

%% Generate random population
Population = Global.Initialization();
[FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,Global.N);
CrowdDis = CalcCrowdingDistance(Population,FrontNo,kc,n_kc,0);

%% Optimization
saveCount = 0;
while Global.NotTermination(Population)
    
    MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
    Offspring  = Global.Variation(Population(MatingPool));
    newPopulation = [Population,Offspring];
    [FrontNo,MaxFNo] = NDSort(newPopulation.objs,newPopulation.cons,numel(newPopulation));
    CrowdDis = CalcCrowdingDistance(newPopulation,FrontNo,kc,n_kc,Global.evaluated/Global.evaluation);
    [Population,FrontNo,CrowdDis] = SortPopulation(newPopulation,kc,n_kc,FrontNo,CrowdDis,Global.N);
    
    if Global.evaluated/Global.evaluation > yita
        [kc,n_kc]=find_knee(Population,FrontNo,Global.evaluated/Global.evaluation,kind);
    end
%     saveDataM(Global,saveCount); %% save data
%     saveCount = saveCount+1;
end
end