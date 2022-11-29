function varargout = PMOP3(Operation,Global,input)
% <problem> <PMOP>
% Benchmarks for Knee identification and optimization
% operator --- EAreal
% g3 l1 h3 k3

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
            g = zeros(N,1);
            if Linkage == 0   
                g = PopDec(:,M:end).^2 - 10.*cos(4.*pi.*PopDec(:,M:end)); %% g3...
                g = 1+10.*(D-M+1)+sum(g,2);
            end
           if  Linkage ~= 0   %% l1 linkage function in g3
                temp = zeros(N,D-M+1);
                for i=1:D-M+1
                   temp(:,i) =  (1+ i./(D-M+1)).*(PopDec(:,M-1+i) - Global.lower(M-1+i)) - PopDec(:,1).*(Global.upper(M-1+i) - Global.lower(M-1+i));   
                end 
                g = temp(:,1:end).^2 - 10.*cos(4.*pi.*temp(:,1:end)); %% g3...
                g = 1+10.*(D-M+1)+sum(g,2);
            end
            
           %% construct the knee functions:
            r = zeros(N,M-1);
            k = zeros(N,1);
            for i=1:M-1
                r(:,i) = 1 + exp(sin(A.*power(PopDec(:,i), B).*pi+pi./2))./(power(2, S).* A);
            end
             Tr = r'; %% transposition
             if M>2
              k = prod(Tr(:,1:N))./(M-1); %% multiplicative
             end
             if M==2
              k = Tr(:,1:N)./(M-1); %% multiplicative
             end
             k = 2.^(k'); %% transposition
            
           %% construct the shape functions:           
            PopObj = ones(N,M);    
            for i=1:M
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g).* k;
            end
            for i=1:M
               for j=1:M-i
                  PopObj(:,i) =  PopObj(:,i).*(1-cos(power(PopDec(:,j), P).*pi./2));
               end
               if i~=1
                   aux = M - i +1;
                   PopObj(:,i)= PopObj(:,i).*(1-sin(power(PopDec(:,aux), P).*pi./2));
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
                r(:,i) = 1 + exp(sin(A.*power(PopDec(:,i), B).*pi+pi./2))./(power(2, S).* A);
            end
           Tr = r'; %% transposition
             if M>2
              k = prod(Tr(:,1:N))./(M-1); %% multiplicative
             end
             if M==2
              k = Tr(:,1:N)./(M-1); %% multiplicative
             end
             k = 2.^(k'); %% transposition
            
           %% construct the shape functions:           
            PopObj = ones(N,M);    
            for i=1:M
                 PopObj(:,i) = PopObj(:,i).*(1.0+1.0).* k; %% min g3 = 1
            end
            for i=1:M
               for j=1:M-i
                  PopObj(:,i) =  PopObj(:,i).*(1-cos(power(PopDec(:,j), P).*pi./2));
               end
               if i~=1
                   aux = M - i +1;
                   PopObj(:,i)= PopObj(:,i).*(1-sin(power(PopDec(:,aux), P).*pi./2));
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

