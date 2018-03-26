%%  Initialization
% Claudia Raducanu and Luka Van de Sype
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka
clc
clearvars
clear all

%%  Determine input
input      =  'Input_AE4424_Ass1P2.xlsx';

% Inputs
[P, R, L, fare, fare_r, demand, capacity, col, delta, Q, costfull, Bpr, ... 
     recap_p, recap_r, recaprate] = matrixsetup1P2(input) ;

not_opt_col = 1;
not_opt_row = 1;
colz = [];
rowz = [];

titer =0;

 %%  Initiate CPLEX model
 %%
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
   %%
    % 1. Capacity constraint
        for i = 1:L
            C11 = zeros(1,DV);
            C12 = zeros(1,DV);
            for pr = 1:DV
                if delta{col(pr,1),1}(i) ~= 0 
                    C11(pr)= 1;
                end
                if delta{col(pr,2),1}(i) ~= 0 && ismember(col(pr,2), col(:,1))
                    C12(pr) = Bpr(col(pr,1),col(pr,2)); 
                end
            end
            C1 = C11 - C12;
            RMP.addRows(Q(i)-capacity(i), C1, inf, sprintf('Capacity_%d',i));
        end

    %%  Execute model 
    %%
    sol = RMP.solve();
    ObjVals = [];
    ObjVals = [ObjVals;RMP.Solution.objval];
    primal  = RMP.Solution.x;   
    dual    = RMP.Solution.dual;
    %RMP.writeModel([model '.lp']);
    

    
    
    
    
    
    

%% SUPER BIG LOOP
%%
while not_opt_col == 1 || not_opt_row == 1

    titer = titer +1;
    disp('-------------------------------------------------');
    disp(['Iteration: ',num2str(titer)]);   
    disp('-------------------------------------------------');
    
    dual    = RMP.Solution.dual;
    
    
%% COLUMN GENERATION
%%
iter=0;
   while not_opt_col == 1
   %% Iterations
   %%
   iter = iter +1;
    disp('-------------------------------------------------');
    disp(['Column iteration: ',num2str(iter)]);   
    disp('-------------------------------------------------');
    size(col)
    
   %% Pricing Problem
   %%
   % Get dual variables 
    pi   = dual(1:L);                 % dual variables of capacity constraints (slack)
    sigma = zeros(P,1);
    if size(dual) > L
        for ro = 1:numel(rowz)
            sigma(rowz(ro))  = dual(L+ro);
        end
    end
    
    checkcol = 0; 
    R = P;
    for re = 1:numel(recap_p)
            pi_j = 0;
            pi_i = 0;
            
            for i = 1:L
                if delta{recap_p(re),1}(i) ~= 0
                    pi_i = pi_i + pi(i);
                end
                if delta{recap_r(re),1}(i) ~= 0
                    pi_j = pi_j + pi(i);   
                end
            end
            
            
            if (fare(recap_p(re)) - pi_i) - (Bpr(recap_p(re),recap_r(re))*(fare(recap_r(re)) - pi_j )) - sigma(recap_p(re)) < 0 
                D = find(recap_r(re)==col(:,2));
                Dp = find(col(D,1)==recap_p(re));
                
                if isempty(Dp) == 1 || isempty(D) == 1
                    checkcol = checkcol + 1;
                    
                    colz = [colz; [recap_p(re) recap_r(re)]]; % only the added columns
                    col = [col; [recap_p(re) recap_r(re)]];   % all columns
                    
                    cost = [cost;costfull(col(recap_p(re),1),col(recap_p(re),2))];
                    
                    % New objective function
                    A       = RMP.Model.A(:,recap_p(re));
                    
                    % add a constraint to every L if the leg is used in p.
                    deltasp = find(delta{recap_p(re),1}(:)~=0);
                    if isempty(deltasp) ~= 0 
                         A(deltasp)= 1;
                    end
                    
                    % substract rr to every L if the leg is used in r.
                    deltasr = find(delta{recap_r(re),1}(:)~=0);
                    if isempty(deltasr) ~= 0
                        A(deltasr)= A(deltasr) - Bpr(recap_p(re),recap_r(re)); 
                    end
                    
                    obj     = costfull(col(recap_p(re),1),col(recap_p(re),2));
                    lb      = [0];
                    ub      = [inf];
                    ctype   = []; 

                    RMP.addCols(obj,A,lb,ub);
                end
            end
    end
   

    
    RMP.solve();
    ObjVals = [ObjVals;RMP.Solution.objval];
    dual = RMP.Solution.dual;
    pi   = dual(1:L);                 % dual variables of capacity constraints (slack)
    sigma = zeros(P,1);
    if size(dual) > L
        for ro = 1:numel(rowz)
            sigma(rowz(ro))  = dual(L+ro);
        end
    end

    
    
    if checkcol == 0
        not_opt_col = 0;
    else 
        not_opt_col = 1;
        not_opt_row = 1;       
    end    
    
   end
     
 
%% ROW GENERATION
%%
     iteration = 0;

   while not_opt_row == 1
  %%  Initiate CPLEX model
  iteration = iteration +1;

    disp('-------------------------------------------------');
    disp(['Row iteration: ',num2str(iteration)]);   
    disp('-------------------------------------------------');
    size(rowz)
  
  DV     =  numel(col(:,1));
  primal = RMP.Solution.x;
     
        %% Separation Problem
    % 2. Demand constraint
    check = 0;
    for p = 1:P
        B = find(col(:,1)==p);
        
        if sum(primal(B)) > demand(p)
            D = find(p==rowz);
            
            if isempty(D) == 1
                check = check + 1;
                C = zeros(1,DV);
                C(B) = 1;
                RMP.addRows(-inf, C, demand(p), sprintf('Demand_%03d',p));
                rowz = [rowz;p];
            end
        end
    end 

    
           %%  Execute model 
    %   Run CPLEX
    RMP.solve();
    ObjVals = [ObjVals;RMP.Solution.objval];

    
%     %   Get dual variables
%     primal  = RMP.Solution.x;   
%     dual    = RMP.Solution.dual;
    
    
    if check == 0
        not_opt_row = 0;
    else
        not_opt_row = 1;
        not_opt_col = 1;
    end       
        
   end
    
end
