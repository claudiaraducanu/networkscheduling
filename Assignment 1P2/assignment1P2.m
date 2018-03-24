%%  Initialization
% Claudia Raducanu and Luka Van de Sype
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
input      =  'Input_AE4424_Ass1P2.xlsx';

% Inputs
[P, R, L, fare, fare_r, demand, capacity, col, delta, Q, costfull, Bpr] ... 
            = matrixsetup1P2(input) ;


iter=0;
%% SUPER BIG LOOP
not_opt_col = 1;
not_opt_row = 1;

titer =0;

while not_opt_col == 1 || not_opt_row == 1

    titer = titer +1;
    disp('-------------------------------------------------');
    disp(['Iteration: ',num2str(titer)]);   
    disp('-------------------------------------------------');
    
    
%% COLUMN GENERATION

   while not_opt_col == 1
   %% Objective function
   iter = iter +1;
    disp('-------------------------------------------------');
    disp(['Column iteration: ',num2str(iter)]);   
    disp('-------------------------------------------------');

  %%  Initiate CPLEX model
        %   Create model 
        model                 =   'Initiate';  % name of model
        RMP                   =    Cplex(model); % define the new model
        RMP.Model.sense       =   'minimize';


        %   Decision variables
        DV                      =  numel(col(:,1));
   
        cost = zeros(numel(col(:,1)),1);
        for p = 1:numel(col(:,1))    
            cost(p) = costfull(col(p,1),col(p,2));                % select only costs you need
        end

        obj                     =   cost ;
        lb                      =   zeros(DV, 1);                                 %Lower bounds
        ub                      =   inf(DV, 1);                                   %Upper bounds             
        
        RMP.addCols(obj, [], lb, ub);
        
   %%  Constraints
    % 1. Capacity constraint
        for i = 1:L
            C11 = zeros(1,DV);
            C12 = zeros(1,DV);
            for pr = 1:DV
                if delta{col(pr,1),1}(i) ~= 0 
                    C11(Tindex(pr) )= 1;
                end
                if delta{col(pr,2),1}(i) ~= 0 && ismember(col(pr,2), col(:,1))
                    C12(Tindex(pr)) = Bpr(col(pr,1),col(pr,2)); 
                end
            end
            C1 = C11 - C12;
            RMP.addRows(Q(i)-capacity(i), C1, inf, sprintf('Capacity_%d',i));
        end
    
   
        
       %%  Execute model 
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
    pi   = dual(1:L);                 % dual variables of capacity constraints (slack)
    if size(dual) <= L
        sigma = zeros(P,1);
    else
        sigma   = dual(L+1:end);
    end
    
        %% Pricing Problem
%    x = zeros(P,P);
    checkcol = 0; 
    R = P;
    %PR = [0 0];
    for p = 1:P
        for r = 1:R
            pi_j = 0;
            pi_i = 0;
            for i = 1:L
                pi_i = pi_i + pi(i)*delta{p,1}(i);
                pi_j = pi_j + pi(i)*delta{r,1}(i);    
            end
            if (fare(p) - pi_i) - (Bpr(p,r)*(fare(r) - pi_j )) - sigma(p) < 0 
                D = find(r==col(:,2));
                if isempty(D)  
                    checkcol = checkcol + 1;
                    col = [col; [p r]]; 
                elseif p ~= col(D,1) 
                    checkcol = checkcol + 1;
                    col = [col; [p r]]; 
                end
            end
        end
    end
    
    if checkcol == 0
        not_opt_col = 0;
    else 
        not_opt_col = 1;
    end    
    checkcol
    
   end
   
   iter
 
%% ROW GENERATION
     iteration = 0;
   while not_opt_row == 1
  %%  Initiate CPLEX model
  iteration = iteration +1;
        %   Create model 
        model                 =   'Initial';  % name of model
        RMP                   =    Cplex(model); % define the new model
        RMP.Model.sense       =   'minimize';
   %% Objective function
        DV                      =  numel(col(:,1));
   
        cost = zeros(numel(col(:,1)),1);
        for p = 1:numel(col(:,1))    
            cost(p) = costfull(col(p,1),col(p,2));                % select only costs you need
        end

        obj                     =   cost ;
        lb                      =   zeros(DV, 1);                                 %Lower bounds
        ub                      =   inf(DV, 1);                                   %Upper bounds             

        RMP.addCols(obj, [], lb, ub);
        
   %%  Constraints
    % 1. Capacity constraint
        for i = 1:L
            C11 = zeros(1,DV);
            C12 = zeros(1,DV);
            for pr = 1:DV
                if delta{col(pr,1),1}(i) ~= 0 
                    C11(Tindex(pr) )= 1;
                end
                if delta{col(pr,2),1}(i) ~= 0 && ismember(col(pr,2), col(:,1))
                    C12(Tindex(pr)) = Bpr(col(pr,1),col(pr,2)); 
                end
            end
            C1 = C11 - C12;
            RMP.addRows(Q(i)-capacity(i), C1, inf, sprintf('Capacity_%d',i));
        end
           
        
       %%  Execute model 
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
%    pi   = dual;                 % dual variables of capacity constraints (slack)
%     sigma_k = dual(1:L,1);     % dual variables of demand constraint
%    sigma   = zeros(P,1);
    
        %% Separation Problem
    % 2. Demand constraint
    check = 0;
    for p = 1:P
        B = find(col(:,1)==p);
        if sum(primal(B)) > demand(p)
            check = check + 1;
            C2 = zeros(1,DV);
            for b = B
                C2(Tindex(b)) = 1;
            end
            RMP.addRows(0, C2, demand(p), sprintf('Demand_%d',p));
        end
    end
    if check == 0
        not_opt_row = 0;
    else
        not_opt_row = 1;
        not_opt_col = 1;
    end
    check
        
        
    
   end    
   
    
 end
    
        %%  Function to return index of decision variables
function out = Tindex(p)
        out = p;
end
    
    
    
    
    
