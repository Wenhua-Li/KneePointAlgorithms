function varargout = DEB2DK(Operation,Global,input)
% <problem> <KneeProblem>
% Comparison of Multiobjective Evolutionary Algorithms for knee
% identification
% operator --- EAreal
    K = 4;
    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 7;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            g = 1+9.*sum(PopDec(:,2:end),2)./(Global.D -1);
            r = 5+10.*(PopDec(:,1)-0.5).*(PopDec(:,1)-0.5)+cos(2.*K.*pi.*PopDec(:,1))./K;
            PopObj(:,1) = g.*r.*sin(PopDec(:,1).*pi./2);
            PopObj(:,2) = g.*r.*cos(PopDec(:,1).*pi./2);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            M = Global.M;
            h = UniformPoint(input,M);            
            X1 = M-1;
            lower    = zeros(1,X1);
            upper    = ones(1,X1);  
            PopDec   = rand(input,X1).*repmat(upper-lower,input,1) + repmat(lower,input,1);
            N = input;
            r = zeros(N,M-1);
            k = zeros(N,1);
            for i=1:M-1     
                r(:,i) = 5+10.*(PopDec(:,i)-0.5).*(PopDec(:,i)-0.5)+cos(2.*K.*pi.*PopDec(:,i))./K;
            end
            k = sum(r,2)./(M-1); %
            
            PopObj(:,1) = k.*sin(PopDec(:,1).*pi./2);
            PopObj(:,2) = k.*cos(PopDec(:,1).*pi./2);
  
            [FrontNo,MaxFNo] = NDSort(PopObj,N);
            f= [];
            t=1;
            for i=1:N   
                if FrontNo(i) ==1
                  f(t,:) = PopObj(i,:);
                  t = t+1;
                end
            end
%             saveData('DEB2DK-',f,K);
            varargout = {f}; 
    end
end