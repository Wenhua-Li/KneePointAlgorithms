function varargout = DEB4DK(Operation,Global,input)
% <problem> <KneeProblem>
% Comparison of Multiobjective Evolutionary Algorithms for knee
% identification
% operator --- EAreal
    K = 1; %% control the number of knees.
    switch Operation
        case 'init'
            Global.M        = 4;
            Global.D        = 12;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            g = 1+9.*sum(PopDec(:,2:end),2)./(Global.D -1);
            g=1;
            
            r1 = 5+10.*(PopDec(:,1)-0.5).*(PopDec(:,1)-0.5)+3.*cos(2.*K.*pi.*PopDec(:,1))./K;
            r2 = 5+10.*(PopDec(:,2)-0.5).*(PopDec(:,2)-0.5)+3.*cos(2.*K.*pi.*PopDec(:,2))./K;
            r3 = 5+10.*(PopDec(:,3)-0.5).*(PopDec(:,3)-0.5)+3.*cos(2.*K.*pi.*PopDec(:,3))./K;
            r = (r1+r2+r3)./3;
       
            PopObj(:,1) = g.*r.*sin(PopDec(:,1).*pi./2).*sin(PopDec(:,2).*pi./2).*sin(PopDec(:,3).*pi./2);
            PopObj(:,2) = g.*r.*sin(PopDec(:,1).*pi./2).*sin(PopDec(:,2).*pi./2).*cos(PopDec(:,3).*pi./2);
            PopObj(:,3) = g.*r.*sin(PopDec(:,1).*pi./2).*cos(PopDec(:,2).*pi./2);
            PopObj(:,4) = g.*r.*cos(PopDec(:,1).*pi./2);
            
            PopCon = [];
            varargout = {input,PopObj,PopCon};
        case 'PF'
             M = Global.M;
            h = UniformPoint(input,M);            
            X1 = M-1;

            lower    = zeros(1,X1);
            upper    = ones(1,X1);  
            PopDec   = rand(input,X1).*repmat(upper-lower,input,1) + repmat(lower,input,1);
            
%             g = 1+9.*sum(PopDec(:,2:end),2)./(Global.D -1);
%             
%             r1 = 5+10.*(PopDec(:,1)-0.5).*(PopDec(:,1)-0.5)+2.*cos(2.*K.*pi.*PopDec(:,1))./K;
%             r2 = 5+10.*(PopDec(:,2)-0.5).*(PopDec(:,2)-0.5)+2.*cos(2.*K.*pi.*PopDec(:,2))./K;
%             r = (r1+r2)./2;
%        
%             PopObj(:,1) = g.*r.*sin(PopDec(:,1).*pi./2).*sin(PopDec(:,2).*pi./2);
%             PopObj(:,2) = g.*r.*sin(PopDec(:,1).*pi./2).*cos(PopDec(:,2).*pi./2);
%             PopObj(:,3) = g.*r.*cos(PopDec(:,1).*pi./2);
   
            N = input;
            r = zeros(N,M-1);
            k = zeros(N,1);
            for i=1:M-1
             r(:,i) = 5+10.*(PopDec(:,i)-0.5).*(PopDec(:,i)-0.5)+3.*cos(2.*K.*pi.*PopDec(:,i))./K;
            end
            k = sum(r,2)./(M-1); %% k1...
            
%             k=k.*(1+9.*sum(PopDec(:,2:end),2)./(Global.D -1));

            PopObj(:,1) = k.*sin(PopDec(:,1).*pi./2).*sin(PopDec(:,2).*pi./2).*sin(PopDec(:,3).*pi./2);
            PopObj(:,2) = k.*sin(PopDec(:,1).*pi./2).*sin(PopDec(:,2).*pi./2).*cos(PopDec(:,3).*pi./2);
            PopObj(:,3) = k.*sin(PopDec(:,1).*pi./2).*cos(PopDec(:,2).*pi./2);
            PopObj(:,4) = k.*cos(PopDec(:,1).*pi./2);

            [FrontNo,MaxFNo] = NDSort(PopObj,N);
            f= [];
            t=1;
            for i=1:N   
                if FrontNo(i) ==1
                  f(t,:) = PopObj(i,:);
                  t = t+1;
                end
            end
%             saveData('DEB3DK-',f,K);
            varargout = {f};
    end
end