%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA120A
% Project Title: Non-dominated Sorting Genetic Algorithm II (NSGA-II)
% Publisher: Yarpiz (www.yarpiz.com)
%
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
%
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;
load DEB2DK4_10

%% Problem Definition

CostFunction=@(x) DEB2DK(x);      % Cost Function

nVar=5;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables

% Number of Objective Functions
nObj=numel(CostFunction(unifrnd(VarMin,VarMax,VarSize)));

n_kc = 0;%初始化knee个数为0
kc=[];

%% NSGA-II Parameters

MaxIt=80;      % Maximum Number of Iterations

nPop=50;        % Population Size

pCrossover=0.7;                         % Crossover Percentage
nCrossover=2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

pMutation=0.4;                          % Mutation Percentage
nMutation=round(pMutation*nPop);        % Number of Mutants

mu=0.02;                    % Mutation Rate

sigma=0.1*(VarMax-VarMin);  % Mutation Step Size


%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];
empty_individual.Reserved=[];

pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
end
% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F,kc,n_kc,0);

% Sort Population
[pop, F]=SortPopulation(pop,kc,n_kc);


%% NSGA-II Main Loop

for it=1:MaxIt
    
    
    %% Crossover 与离得最近的knee交叉
    popc=repmat(empty_individual,nCrossover/2,2);
    for k=1:nCrossover/2
        
        i1=randi([1 nPop]);
        p1=pop(i1);
        
        if 0%n_kc~=0
            Costs=[p1.Cost];
            knee_center=[kc.Cost];
            dis = zeros(n_kc,1);
            for j=1:n_kc
                dis(j)=(Costs(1)-knee_center(1,j))^2+(Costs(2)-knee_center(2,j))^2;
                dis(j) = dis(j)^0.5;
            end
            [d,index] = min(dis);
            p2=kc(index);
        else
            i2=randi([1 nPop]);
            p2=pop(i2);
        end
        [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position);
        
        popc(k,1).Cost=CostFunction(popc(k,1).Position);
        popc(k,2).Cost=CostFunction(popc(k,2).Position);
        
    end
    popc=popc(:);
    
    %% Mutation
    popm=repmat(empty_individual,nMutation,1);
    for k=1:nMutation
        
        i=randi([1 nPop]);
        p=pop(i);
        
        popm(k).Position=Mutate(p.Position,mu,sigma);
        
        %% 越界处理
        chrom=[popm(k).Position];
        chrom(find(chrom>1))=1;
        chrom(find(chrom<0))=0;
        popm(k).Position=chrom;
        
        popm(k).Cost=CostFunction(popm(k).Position);
        
    end
    
    % Merge
    pop=[pop
        popc
        popm]; %#ok
    
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);
    

    
    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F,kc,n_kc,it);
    
    % Sort Population
    pop=SortPopulation(pop,kc,n_kc);
    
    % Truncate
    pop=pop(1:nPop);
    
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);
    
        %% find knee center
    if it/MaxIt>0.3
        [kc,n_kc] = find_knee(pop,F,it);
    end
    
    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F,kc,n_kc,it);
    
    % Sort Population
    [pop, F]=SortPopulation(pop,kc,n_kc);
    
    % Store F1
    F1=pop(F{1});
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    
    % Plot F1 Costs
    figure(1);
    PlotCosts(PF,F1,kc);
    
    pause(0.01);
    
end

%% Results

