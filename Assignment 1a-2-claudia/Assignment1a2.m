%%  Initialization
% Claudia Raducanu and Luka Van de Sype
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
%   Select input file and sheet

input      =  'Input_AE4424_Ass1P1.xlsx';

%% Inputs

    [Nodes, K, cost, capacity, origin, destination, demand,s,t] = ...
    read_arc(input);

    A = size(s,1);  % number of arcs
    
    % Set: P^k set of all paths for commodity k (for all in K)
    
    P   = K;         % All paths available at beginning for all k       
    P_k = ones(1,K); % Paths for commodity k at beginning  
   
    
    capacity    = nonzeros(capacity); % RHS
    cost_arc    = nonzeros(cost);     % Cost for each arc
    
%% A: Define initial set of columns

    %Set of shortest path for commodity k
    sp.dist = cell(P,1);
    sp.path = cell(K,1);
    sp.pred = cell(P,1);
    sp.arcs   = mat2cell([s t],ones(size(s,1),1));
    
    
    
    
    for k = 1:K
        [sp.dist{k,1}{P_k(k),1}, sp.path{k,1}{P_k(k),1}] = graphshortestpath(cost,origin(k),destination(k),...
                        'Directed', true);
    end
    
    GG = 
    
    kkk     = shortestpath(GG,origin,destination)
    
    
    f_p  = ones(K,1);
    
    dual_feasibility   = 0; 
    primal_feasibility = 0; 
    
    stop               = 0; 
    
i = 0;


%% B: Solve RMP
while i < 2
    %% Parameters
    
    
    i= i+1;
    disp('-------------------------------------------------');
    disp(['Iteration: ',num2str(i)]);   
    disp('-------------------------------------------------');
    
    
    % delta 
    
    param.delta = cell(K,1);
    
    for k=1:K
        param.delta{k,1} = zeros(P_k(k),A);
        for p = 1:P_k(k)
            for a = 1:A
                pathlength = size(sp.path{k,1}{p,1},2);
                if pathlength > 2 
                    ainp = repelem(sp.path{k,1}{p,1},2);
                    ainp = ainp(2:end-1);
                    ainp = transpose(reshape(ainp,2,pathlength-1));
                    for ac = 1:(pathlength-1) 
                        if sum(ainp(ac,:) == sp.arcs{a,1}) == 2
                            param.delta{k,1}(p,a)  = 1;
                        end
                    end
                else
                    if sum(sp.path{k,1}{p,1} == sp.arcs{a,1}) == 2
                     param.delta{k,1}(p,a)  = 1;
                    end
                end
            end
         end
    end 
    
    if i == 1
        transport = zeros(1,a); 
        
        for a = 1:A
            for k = 1:K
                for p = 1:P_k(k)
                    transport(a) = transport(a) + demand(k)*param.delta{k,1}(p,a);
                end
            end
        end
        
        cap_idx = find(transport > capacity');
        S       = size(cap_idx,2); % Number of slack variables in problem. 
    end
    
    %%  Initiate CPLEX model
        %   Create model 
        model                 =   'Initial';  % name of model
        RMP                   =    Cplex(model); % define the new model
        RMP.Model.sense       =   'minimize';

        %   Decision variables
        DV                      =  P+S;  % Number of Decision Var (xijk)
   %% Objective function
   
        obj.c_p = [];
        for k = 1:K
            for p = 1:P_k(k)
                obj.c_p       = [obj.c_p ; sp.dist{k,1}{p,1}*demand(k)];
            end
        end
        obj.M                      =   ones(S,1)*1000;
        obj.o                      =   [obj.c_p; obj.M] ;
        obj.lb                     =   zeros(DV, 1);                                 %Lower bounds
        obj.ub                     =   inf(DV, 1);                                   %Upper bounds         
        

        RMP.addCols(obj.o, [], obj.lb, obj.ub);

    %%  Constraints
    % 1. Commodity constraint
        for k = 1:K
            C1 = zeros(1,DV);
            for p = 1:P_k(k)
                C1(Findex(P_k,k,p)) = 1;   
            end
            RMP.addRows(1, C1, 1);
        end
    % 
    % 2. Bundle constraint
        for a = 1:A
            C2 = zeros(1,DV);
            for k = 1:K
                for p = 1:P_k(k)
                    C2(Findex(P_k,k,p)) = demand(k)*param.delta{k,1}(p,a);
                end
            end 
            if cap_idx == a 
                C2(end) = -1; 
            end
            RMP.addRows(0, C2, capacity(a));
        end

    %%  Execute model
    %RMP.Param.timelimit.Cur           = 120;         %max time in seconds
 
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
    %f_p     = primal(1:P,1);
    
    
    dual    = RMP.Solution.dual;
    pi_ij   = dual(K+1:end,1); % dual variables of bundle constraints (slack)
    sigma_k = dual(1:K,1);     % dual variables of 
       
    
%% C: Pricing Problem
  
    cost_arc     = cost_arc - pi_ij;
    
    cost       = sparse(s,t,cost_arc,Nodes,Nodes);
    
    for k = 1:K
        [sp.dist{k,i+1}, sp.path{k,i+1}] = ...
            graphshortestpath(sparse(cost),origin(k),destination(k),...
                        'Directed', true);
    end
    
%% D: Determine column(s) to add. 
    
    % Determine for which commodity to add path. 
    
    col_idx   = find(cell2mat(sp.dist(:,i+1)) < sigma_k./demand'); 
    
    
    if isempty(col_idx) 
        dual_feasibility = 1; 
    end
    
    disp(['Dual Feasibility: ',num2str(dual_feasibility)]);   
    
    for o = 1:size(col_idx,1) 
        sp.path{col_idx(o),1}{end+1,1} = sp.path{col_idx(o),i+1};
        sp.dist{col_idx(o),1}{end+1,1} = sp.dist{col_idx(o),i+1};
    end
    
    % Update set
    P             = P + size(col_idx,1); % total number of paths
    P_k(col_idx') = P_k(col_idx')+1;     % total number of paths for commodity k    
    
    stop =  dual_feasibility + primal_feasibility;
    
    disp(['Stop: ',num2str(stop)]);   
    disp('-------------------------------------------------');
    
    
end
   
    %%  Function to return index of decision variables
function out = Findex(P_k,k,p) 
    if k == 1 
        out = p;
    else
        out = sum(P_k(1:(k-1))) + p;  
    end
end



