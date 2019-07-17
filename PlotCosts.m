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

function PlotCosts(PF,pop,kc)
Costs=[pop.Cost];
pf=[PF.Cost];
plot(pf(1,:),pf(2,:),'go','MarkerSize',2);


plot(Costs(1,:),Costs(2,:),'ro','MarkerSize',2);


if ~isempty(kc)
    kc_cost=[kc.Cost];
    plot(kc_cost(1,:),kc_cost(2,:),'b+','MarkerSize',8);
end

xlabel('1^{st} Objective');
ylabel('2^{nd} Objective');
title('Non-dominated Solutions (F_{1})');
grid on;
hold off
end