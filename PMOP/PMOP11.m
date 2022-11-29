function varargout = PMOP11(Operation,Global,input)
% <problem> <PMOP>
% Benchmarks for Knee identification and optimization
% operator --- EAreal
% g1 & g2 l2 h2 k2

    A = 4; %% control the number of knees, and width of knee region.
    B = 1; %% control the bias
    S = 2; %% control the depth of knee
    P = 1;  %% bias in shape function
    Linkage = 0; %% control the linkage function. 

switch Operation
        case 'init'
            Global.M        = 3;
            Global.D        = 12;
            Global.lower    = zeros(1,Global.D);
            X1 = Global.M-1;
            X2 = Global.D - Global.M + 1;
            Global.upper    = [ones(1,X1),10.*ones(1,X2)];  
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
            Global.operator = @EAreal;
            
        case 'value'
            PopDec = input;
            [N,D]  = size(PopDec);
            M      = Global.M;
                       
           %% construct g function:
            g1 = zeros(N,1);
            g2 = zeros(N,1);
            if Linkage == 0  
                %% calc g1
                [g1,~] = max(abs(PopDec(:,M:end)),[],2); %% g1...
                %% calc g2
                temp1 = zeros(N,D-M+1);
                temp1 = PopDec(:,M:end).^2;
                g2 = sum(temp1,2); %% g2...    
            end
                 
            if  Linkage ~= 0   %% l2 linkage function in g1
                temp = zeros(N,D-M+1);
                for i=1:D-M+1
                   temp(:,i) =  (1+ cos(0.5.*pi.*i./(D-M+1))).*(PopDec(:,M-1+i) - Global.lower(M-1+i)) - PopDec(:,1).*(Global.upper(M-1+i) - Global.lower(M-1+i));   
                end 
                [g1,~] = max(abs(temp(:,1:end)),[],2); %% g1...
                temp = temp.^2;
                g2 = sum(temp,2); %% g2...       
            end
            
           %% construct the knee functions:
            r = zeros(N,M-1);
            k = zeros(N,1);
            for i=1:M-1
                r(:,i) = 1 + exp(cos(A.*power(PopDec(:,i), B).*pi+pi./2))./(power(2, S).* A);
            end
            Tr = r'; %% transposition
             if M>2
              k = prod(Tr(:,1:N))./(M-1); %% multiplicative
             end
             if M==2
              k = Tr(:,1:N)./(M-1); %% multiplicative
             end
             k =log(1./(k')+1); %% transposition
            
           %% construct the shape functions:           
            PopObj = ones(N,M);    
            for i=1:M
                if mod(i,2) == 0
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g1).* k;
                end
                if mod(i,2) == 1 
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g2).* k;
                end
            end
            for i=1:M
               for j=1:M-i
                  PopObj(:,i) =  PopObj(:,i).*cos(power(PopDec(:,j), P).*pi./2);
               end
               if i~=1
                   aux = M - i +1;
                   PopObj(:,i)= PopObj(:,i).*sin(power(PopDec(:,aux), P).*pi./2);
               end
            end

            %% construct the knee functions: f = (1+g)*k*h:             
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
                r(:,i) = 1 + exp(cos(A.*power(PopDec(:,i), B).*pi+pi./2))./(power(2, S).* A);
            end
            Tr = r'; %% transposition
             if M>2
              k = prod(Tr(:,1:N))./(M-1); %% multiplicative
             end
             if M==2
              k = Tr(:,1:N)./(M-1); %% multiplicative
             end
             k =log(1./(k')+1); %% transposition
            
            PopObj = ones(N,M);    
            for i=1:M
                  PopObj(:,i) = PopObj(:,i).* k;
            end
            for i=1:M
               for j=1:M-i
                  PopObj(:,i) =  PopObj(:,i).*cos(power(PopDec(:,j), P).*pi./2);
               end
               if i~=1
                   aux = M - i +1;
                   PopObj(:,i)= PopObj(:,i).*sin(power(PopDec(:,aux), P).*pi./2);
               end
            end

            [FrontNo,MaxFNo] = NDSort(PopObj,N);
            f= [];
            t=1;
            for i=1:N   
                if FrontNo(i) ==1
                  f(t,:) = PopObj(i,:);
                  t = t+1;
                end
            end
            
            varargout = {f};  
    end
end

