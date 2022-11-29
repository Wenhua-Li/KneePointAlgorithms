function varargout = PMOP14(Operation,Global,input)
% <problem> <PMOP>
% Benchmarks for Knee identification and optimization
% operator --- EAreal
% g6 & g8 l1 h3 k3

    A = 2; %% control the number of knees, and width of knee region.
    B = 1; %% control the bias
    S = -1; %% control the depth of knee
    P = 1;  %% bias in shape function
    Linkage = 0; %% control the linkage function. 

switch Operation
        case 'init'
            Global.M        = 3;
            Global.D        = 12;
            Global.lower    = zeros(1,Global.D);
            if Global.M==2
                disp('it is degenerated knee problem, and M must be larger than 2.');
            end
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
            g6 = zeros(N,1);
            g8 = zeros(N,1);
            if Linkage == 0  
                %% calc g6
                temp1 = zeros(N,D-M+1);
                temp1 = PopDec(:,M:end).^2 - 10.*cos(2.*pi.*PopDec(:,M:end)) + 10;
                g6 = sum(temp1,2); % g6
                %% calc g8
                temp2 = zeros(N,D-M+1);
                temp3 = zeros(N,D-M+1); 
                temp2 = sum(PopDec(:,M:end).^2,2)./(D-M+1);
                temp3 = sum(cos(2.*pi.*PopDec(:,M:end)),2)./(D-M+1);
                g8 = -20.*exp(-0.2.*sqrt(temp2)) - exp(temp3) + 20 + exp(1); %% g8... 
            end
           if  Linkage ~= 0   %% l1 linkage function in g6 g8
                temp = zeros(N,D-M+1);
                for i=1:D-M+1
                   temp(:,i) =  (1+ i./(D-M+1)).*(PopDec(:,M-1+i) - Global.lower(M-1+i)) - PopDec(:,1).*(Global.upper(M-1+i) - Global.lower(M-1+i));   
                end 
                %% calc g6
                temp1 = zeros(N,D-M+1);
                temp1 = temp(:,1:end).^2 - 10.*cos(2.*pi.*temp(:,1:end)) + 10;
                g6 = sum(temp1,2); % g6
                %% calc g8
                temp2 = zeros(N,D-M+1);
                temp3 = zeros(N,D-M+1); 
                temp2 = sum(temp(:,1:end).^2,2)./(D-M+1);
                temp3 = sum(cos(2.*pi.*temp(:,1:end)),2)./(D-M+1);
                g8 = -20.*exp(-0.2.*sqrt(temp2)) - exp(temp3) + 20 + exp(1); %% g8...
            end
            
           %% construct the -- degenerated -- knee functions: 
            r = zeros(N,M-2);
            k = zeros(N,1);
            for i=1:M-2
                r(:,i) = 1 + exp(sin(A.*power(PopDec(:,i), B).*pi+pi./2))./(power(2, S).* A);
            end
            Tr = r'; %% transposition
             if M>3
              k = prod(Tr(:,1:N))./(M-2); %% multiplicative
             end
             if M==3
              k = Tr(:,1:N)./(M-2); %% multiplicative
             end
             k = sqrt(k'); %% transposition
            
           %% construct the shape functions:           
            PopObj = ones(N,M);    
            for i=1:M
                if mod(i,2) == 1
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g6).* k;
                end
                if mod(i,2) == 0 
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g8).* k;
                end
            end
             for i=1:M
               for j=1:M-i
                  PopObj(:,i) =  PopObj(:,i).*power(PopDec(:,j), P);
               end
               if i~=1
                   aux = M - i +1;
                   PopObj(:,i)= PopObj(:,i).*(1 - power(PopDec(:,aux), P));
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
                     
            %% construct the -- degenerated -- knee functions: 
            r = zeros(N,M-2);
            k = zeros(N,1);
            for i=1:M-2
                r(:,i) = 1 + exp(sin(A.*power(PopDec(:,i), B).*pi+pi./2))./(power(2, S).* A);
            end
            Tr = r'; %% transposition
             if M>3
              k = prod(Tr(:,1:N))./(M-2); %% multiplicative
             end
             if M==3
              k = Tr(:,1:N)./(M-2); %% multiplicative
             end
             k = sqrt(k'); %% transposition
            
           %% construct the shape functions:           
            PopObj = ones(N,M);    
            for i=1:M
                 PopObj(:,i) = PopObj(:,i).*(1.0).* k;
            end
             for i=1:M
               for j=1:M-i
                  PopObj(:,i) =  PopObj(:,i).*power(PopDec(:,j), P);
               end
               if i~=1
                   aux = M - i +1;
                   PopObj(:,i)= PopObj(:,i).*(1 - power(PopDec(:,aux), P));
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
