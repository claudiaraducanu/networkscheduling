%%  Initialization
% Claudia Raducanu and Luka Van de Sype
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
%   Select input file and sheet

input      =  'Input_AE4424_Ass1Verification.xlsx';

%% Inputs

    [Nodes, K, cost, capacity, origin, destination, demand,s,t] = ...
    read_arc_v(input);

    A = size(s,1);  % number of arcs
    
    % Set: P^k set of all paths for commodity k (for all in K)
    
    P   = K;         % All paths available at beginning for all k       
    P_k = ones(1,K); % Paths for commodity k at beginning  
    
    capacity    = nonzeros(capacity); % RHS
    cost_arc    = nonzeros(cost);  % Cost for each arc
    
%% A: Define initial set of columns

    %Set of shortest path for commodity k
    sp.dist = cell(P,1);
    sp.path = cell(P,1);
    sp.pred = cell(P,1);
    sp.arcs   = mat2cell([s t],ones(size(s,1),1));
    
    for k = 1:K
        [sp.dist{k,1}, sp.path{k,1}, sp.pred{k,1}] = graphshortestpath(sparse(cost),origin(k),destination(k),...
                        'Directed', true);
    end

i = 0;
%% B: Solve RMP
while i <= 1   
    %% Parameters
    i= i+1;
    disp('-------------------------------------------------');
    disp(['Iteration: ',num2str(i)]);   
    disp('-------------------------------------------------');
    % delta 
    
    d_k_p_ij = cell(K,1);
    
    for k=1:K
        d_k_p_ij{k,1} = zeros(P_k(k),A);
        for p = 1:P_k(k)
            for a = 1:A
                pathlength = size(sp.path{k,p},2);
                if pathlength > 2 
                    ainp = repelem(sp.path{k,p},2);
                    ainp = ainp(2:end-1);
                    ainp = transpose(reshape(ainp,2,pathlength-1));
                    for ac = 1:(pathlength-1) 
                        if sum(ainp(ac,:) == sp.arcs{a,1}) == 2
                            d_k_p_ij{k,1}(p,a)  = 1;
                        end
                    end
                else
                    if sum(sp.path{k,1} == sp.arcs{a,1}) == 2
                     d_k_p_ij{k,1}(p,a)  = 1;
                    end
                end
            end
         end
    end 


    %%  Initiate CPLEX model
        %   Create model 
        model                 =   'Initial';  % name of model
        RMP                   =    Cplex(model); % define the new model
        RMP.Model.sense       =   'minimize';

        %   Decision variables
        DV                      =  P+A;  % Number of Decision Var (xijk)
   %% Objective function
        c_p = [];
        for k = 1:K
            for p = 1:P_k(k)
                c_p       = [c_p ; sp.dist{k,p}*demand(k)];
            end
        end
        M                       =   ones(A,1)*1000;
        obj                     =   [c_p; M] ;
        lb                      =   zeros(DV, 1);                                 %Lower bounds
        ub                      =   inf(DV, 1);                                   %Upper bounds              

        RMP.addCols(obj, [], lb, ub);

    %%  Constraints
    % 1. Commodity constraint
        for k = 1:K
            C1 = zeros(1,DV);
            for p = 1:P_k(k)
                C1(Findex(P_k,k,p)) = 1;   
            end
            RMP.addRows(1, C1, 1,sprintf('Commodity_%d',k));
        end
    % 
    % 2. Bundle constraint
        for a = 1:A
            C2 = zeros(1,DV);
            for k = 1:K
                for p = 1:P_k(k)
                    C2(Findex(P_k,k,p)) = demand(k)*d_k_p_ij{k,1}(p,a);
                end
            end
            C2(Sindex(a,P)) = -1; 
            RMP.addRows(0, C2, capacity(a),sprintf('Bundle_%d',a));
        end

    %%  Execute model
    RMP.Param.mip.limits.nodes.Cur    = 1e+8;        %max number of nodes to be visited (kind of max iterations)
    RMP.Param.timelimit.Cur           = 120;         %max time in seconds
 
    %   Run CPLEX
    RMP.solve();
    RMP.writeModel([model '.lp']);
    
    %   Get dual variables
    dual    = RMP.Solution.dual;
    pi_ij   = dual(K+1:end,1); % dual variables of bundle constraints (slack)
    sigma_k = dual(1:K,1);     % dual variables of 
    
    %  Modify costs at the arcs

    for k=1:K 
        for a = 1:A
            for p=1:P_k(k)
                sumarc(p,a,k) = demand(k)*(cost_arc(a)-pi_ij(a))*d_k_p_ij{k,1}(p,a);
            end
        end
    end
    
    arc_cost = find(sum(sumarc,2) <= sigma_k);
    
    
%% C: Pricing Problem
  
    cost_arc     = cost_arc - pi_ij;
    
    cost       = sparse(s,t,cost_arc,Nodes,Nodes);
    
    for k = 1:size(arc_cost)
        [sp.dist{arc_cost(k),i+1}, sp.path{arc_cost(k),i+1}, sp.pred{arc_cost(k),i+1}] = ...
            graphshortestpath(sparse(cost),origin(k),destination(k),...
                        'Directed', true);
    end
    
%% D: Determine column(s) to add. 
    
    % Determine for which commodity to add path. 
    col_idx  = find(cell2mat(sp.dist(:,i)) < sigma_k./demand'); 
    
    
    
%     for m = 1:size(col_idx,1)
%         for n = 1:size(arc_cost,1)
%             if col_idx(m,1) == arc_cost(n,1)
%                 col_idx(m,1) = 0;
%             end
%         end
%     end
    
    
    % Update set
    P             = P + size(col_idx,1); % total number of paths
    P_k(col_idx') = P_k(col_idx')+1;     % total number of paths for commodity k
  
    
end
    
    %%  Function to return index of decision variables
function out = Findex(P_k,k,p) 
    if k == 1 
        out = p;
    else
        out = sum(P_k(1:(k-1))) + p;  
    end
end

function out = Sindex(a,P)
    out = P + a;  
end



