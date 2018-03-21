%%  Initialization
% Claudia Raducanu and Luka Van de Sype
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

% clearvars
% clear all

%%  Determine input
%   Select input file and sheet

input      =  'Input_Example.xlsx';

%% Inputs
[P, L, fare, demand , capacity, pathflights, flightnrs] ... 
            = matrixsetup1P2_EXAMPLE(input) ;

Bpr = xlsread('Bpr_EXAMPLE',1,'A1:I8') ;

col = 1;

% The binary value: Delta
delta = cell(P,1); % For each Path all Flights are checked: delta{p,1}(i)
for p = 1:P
    delta{p,1} = zeros(L,1);
    for i = 1:L
        if pathflights(p,1) == flightnrs(i) || pathflights(p,2) == flightnrs(i)
            delta{p,1}(i) = 1;            
        end
    end
end

% The daily unconstrained demand on flight(i): Q
Q = zeros(L,1);
for i = 1:L
    a = zeros(P,1);
    for p = 1:P
        a(p) = demand(p)*delta{p,1}(i);
    end
    Q(i) = sum(a);
end

% adding the ficticious fare and element
fare_r = [0;fare] ; % adding the ficticious recapture fare
R = P +1;


costfull = zeros(P,R);      % full cost matrix
for p = 1:P
    for r = 1:R
        costfull(p,r) = fare(p) - Bpr(p,r)*fare_r(r);
    end
end

iter=0;
%% B: Solve RMP
 while iter < 2
    %% Parameters
    
    iter = iter+1;
    disp('-------------------------------------------------');
    disp(['Iteration: ',num2str(iter)]);   
    disp('-------------------------------------------------');

    
cost = costfull(:,col);                % select only costs you need
cost = reshape(cost,numel(col)*P,1);   % reshape for cplex


  %%  Initiate CPLEX model
        %   Create model 
        model                 =   'Initial';  % name of model
        RMP                   =    Cplex(model); % define the new model
        RMP.Model.sense       =   'minimize';

        %   Decision variables
        DV                      =  numel(col)*P;

   %% Objective function

        obj                     =   cost ;
        lb                      =   zeros(DV, 1);                                 %Lower bounds
        ub                      =   inf(DV, 1);                                   %Upper bounds              

        RMP.addCols(obj, [], lb, ub);
        
   %%  Constraints
    % 1. Capacity constraint
        for i = 1:L
            C11 = zeros(1,DV);
            C12 = zeros(1,DV);
            for p = 1:P
                for r = 1:col
                    if delta{p,1}(i) ~= 0 && p ~= r-1
                        C11(Tindex(p,r)) = 1;
                        if r == 1
                            C12 = 0;
                        else
                            C12(Tindex(r-1,p+1)) = Bpr(r-1,p+1); % as the first column is 'ficticious'
                        end
                    end
                end
            end
            C1 = C11 - C12;
            RMP.addRows(Q(i)-capacity(i), C1, inf, sprintf('Capacity_%d',i));
        end
        
%     % 2. Demand constraint
%         for p = 1:P
%             C2 = zeros(1,DV);
%             for r = 1:R
%                 C2(Tindex(r,p)) = 1;
%             end
%             RMP.addRows(0, C2, demand(p), sprintf('Demand_%d',p));
%         end            
%         
       %%  Execute model
% %     RMP.Param.mip.limits.nodes.Cur    = 1e+8;        %max number of nodes to be visited (kind of max iterations)
% %     RMP.Param.timelimit.Cur           = 120;         %max time in seconds
 
    %   Run CPLEX
    RMP.solve();
    RMP.writeModel([model '.lp']);
    
    if RMP.Solution.status == 1 
        primal_feasibility = 1; 
    end
    
    disp('-------------------------------------------------');
    disp(['Primal Feasibility: ',num2str(primal_feasibility)]);  
    
    %   Get dual variables
    primal  = RMP.Solution.x;   
    dual    = RMP.Solution.dual;
    pi   = dual; % dual variables of capacity constraints (slack)
%     sigma_k = dual(1:L,1);     % dual variables of demand constraint
    sigma   = zeros(P,1);
        %% Pricing Problem
    x = zeros(P,P);
    R = P;
    for p = 1:P
        for r = 1:R
            pi_j = 0;
            pi_i = 0;
            
            for i = 1:L
            pi_i = pi_i + pi(i)*delta{p,1}(i);
            pi_j = pi_j + pi(i)*delta{r,1}(i);
                
            end
            if (fare(p) - pi_i) - (Bpr(p,r+1)*(fare(r) - pi_j )) - sigma(p) < 0
                col = [col, r];
            end
        end
    end
    
    
    
    
 end
    
        %%  Function to return index of decision variables
function out = Tindex(p,r)
        P=8;
        out = P*(r-1) + p;
end
    
    
    
    
    
