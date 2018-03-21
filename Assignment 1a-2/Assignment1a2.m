%%  Initialization
% Claudia Raducanu and Luka Van de Sype
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
%   Select input file and sheet

input      =  'Input_AE4424_Ass1Verification.xlsx';

%% Inputs

    [nodes, K, arcs, origin, destination, demand,slack(:,1),slack(:,2),capacity] = ...
    read_arc_v(input);

    A = size(slack,1);  % number of arcs
    
    % Set: P^k set of all paths for commodity k (for all in K)
    
    P   = K;         % All paths available at beginning for all k       
    P_k = ones(1,K); % Paths for commodity k at beginning  
    
%% A: Define initial set of columns

    %Set of shortest path for commodity k
    sp.path = cell(K,1);
    sp.dist = cell(K,1);
    sp.delta = cell(K,1);
    sp.slack = zeros(A,P);
    
    for k = 1:K
        sp.path{k,1}{P_k(k),1} = shortestpath(arcs,origin(k),destination(k));
        eidx = findedge(arcs, sp.path{k,1}{P_k(k),1}(1:end-1), sp.path{k,1}{P_k(k),1}(2:end));
        sp.delta{k,1}                       = zeros(P_k(k),A); 
        sp.delta{k,1}(P_k(k),eidx)          = demand(k);
        sp.slack(eidx,k)                    = demand(k);
    end
    
    slack_idx  = find(sum(sp.slack,2) > capacity);
    
    S          = size(slack_idx,1);     
    
    %% Iteration stop conditions
    
    dual_feasibility   = 0; 
    primal_feasibility = 0; 
    
    stop               = 0; 
    
    i = 0;


%% B: Solve RMP
%while i < 4
    %% Parameters
    
    i= i+1;
    disp('-------------------------------------------------');
    disp(['Iteration: ',num2str(i)]);   
    disp('-------------------------------------------------');

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
                eidx = findedge(arcs, sp.path{k,1}{p,1}(1:end-1), sp.path{k,1}{p,1}(2:end));
                sp.dist{k,1}(p,1) = sum(arcs.Edges.Weight(eidx));
            end
            obj.c_p = [ obj.c_p ; sp.dist{k,1}(1:end,1) ];
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
 
%     % 2. Bundle constraint
        for a = 1:A
            C2 = zeros(1,DV);
            for k = 1:K
                for p = 1:P_k(k)
                    C2(Findex(P_k,k,p)) = sp.delta{k,1}(p,a);
                end
            end 
            if a == 6 
                C2(end) = -1; 
            end
            RMP.addRows(0, C2, capacity(a));
        end

%     %%  Execute model
%     %RMP.Param.timelimit.Cur           = 120;         %max time in seconds
%  
%     %   Run CPLEX
%     RMP.solve();
%     RMP.writeModel([model '.lp']);
%     
%     if RMP.Solution.status == 1 
%         primal_feasibility = 1; 
%     end
%     
%     disp('-------------------------------------------------');
%     disp(['Primal Feasibility: ',num2str(primal_feasibility)]);   
%     
%     %   Get dual variables
%     primal  = RMP.Solution.x;
%     %f_p     = primal(1:P,1);
%     
%     
%     dual    = RMP.Solution.dual;
%     pi_ij   = dual(K+1:end,1); % dual variables of bundle constraints (slack)
%     sigma_k = dual(1:K,1);     % dual variables of 
%        
%     
% %% C: Pricing Problem
%   
%     cost_arc     = cost_arc - pi_ij;
%     
%     cost       = sparse(s,t,cost_arc,Nodes,Nodes);
%     
%     for k = 1:K
%         [sp.dist{k,i+1}, sp.path{k,i+1}] = ...
%             graphshortestpath(sparse(cost),origin(k),destination(k),...
%                         'Directed', true);
%     end
%     
% %% D: Determine column(s) to add. 
%     
%     % Determine for which commodity to add path. 
%     
%     col_idx   = find(cell2mat(sp.dist(:,i+1)) < sigma_k./demand'); 
%     
%     
%     if isempty(col_idx) 
%         dual_feasibility = 1; 
%     end
%     
%     disp(['Dual Feasibility: ',num2str(dual_feasibility)]);   
%     
%     for o = 1:size(col_idx,1) 
%         sp.path{col_idx(o),1}{end+1,1} = sp.path{col_idx(o),i+1};
%         sp.dist{col_idx(o),1}{end+1,1} = sp.dist{col_idx(o),i+1};
%     end
%     
%     % Update set
%     P             = P + size(col_idx,1); % total number of paths
%     P_k(col_idx') = P_k(col_idx')+1;     % total number of paths for commodity k    
%     
%     stop =  dual_feasibility + primal_feasibility;
%     
%     disp(['Stop: ',num2str(stop)]);   
%     disp('-------------------------------------------------');
%     
%     
% end
%    
%%  Function to return index of decision variables

function out = Findex(P_k,k,p) 
    if k == 1 
        out = p;
    else
        out = sum(P_k(1:(k-1))) + p;  
    end
end



