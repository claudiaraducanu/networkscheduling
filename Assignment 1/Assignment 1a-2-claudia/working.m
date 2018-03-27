%%  Initialization
% Claudia Raducanu and Luka Van de Sype
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
%   Select input file and sheet     
     
%input      =  'Input_AE4424_Ass1Verification.xlsx';
input      =  'Input_AE4424_Ass1P1.xlsx';

%% Inputs

    [nodes, K, arcs, origin, destination, demand,od(:,1),od(:,2),capacity] = ...
    read_arc(input);

    A = size(od,1);  % number of arcs
    
    % Set: P^k set of all paths for commodity k (for all in K)
    
    P   = K;         % All paths available at beginning for all k       
    P_k = ones(1,K); % Paths for commodity k at beginning  
    
%% A: Define initial set of columns

    %Set of shortest path for commodity k
    sp.path = cell(K,1);
    sp.dist = cell(K,1);
    sp.delta = zeros(A,P);
    
    for k = 1:K
        sp.path{k,1}  = shortestpath(arcs,origin(k),destination(k));
        eidx = findedge(arcs, sp.path{k,1}(1:end-1), sp.path{k,1}(2:end));
        sp.delta(eidx,k)                    = demand(k);
    end
    
    
    slack_idx                          = find(sum(sp.delta,2) > capacity);
    S                                  = size(slack_idx,1);   
    sp.colslack                        = zeros(A,S);
    
    for s = 1:S
        sp.colslack(slack_idx(s),s) = -1;
    end
    
    c.ineq = [ sp.delta sp.colslack];
    c.eq   = [ eye(P,P) zeros(P,S)];
    
    
    obj.c_p = [];

    for k = 1:K
        for p = 1:P_k(k)
            eidx = findedge(arcs, sp.path{k,p}(1:end-1), sp.path{k,p}(2:end));
            sp.dist{k,1}(p,1) = sum(arcs.Edges.Weight(eidx));
        end
        obj.c_p = [ obj.c_p ; demand(k).*sp.dist{k,1}(1:end,1) ];
    end 
    
    obj.M                      =   ones(S,1)*1000;
    obj.o                      =   [obj.c_p; obj.M];
    
    %% Iteration stop conditions
    
    dual_feasibility   = 0; 
    primal_feasibility = 0; 
    
    stop               = 0; 
    
    i = 0;
   
%% B: Solve RMP

    %% Parameters
   
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
 
        RMP.addCols(obj.o);

    %%  Constraints
    
     %   1. Commodity constraint
    for k= 1:K
        C1 = c.eq(k,:);
        RMP.addRows(1, C1, 1);
    end

    % 2. Bundle constraint
    for a = 1:A
        C2 = c.ineq(a,:);
        RMP.addRows(0, C2, capacity(a,1));
    end
     %%  Execute model
    %RMP.Param.timelimit.Cur           = 120;         %max time in seconds  
    
while i < 5  
    i=i+1;

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
    pi_ij   = dual(K+1:end,1); % dual variables of bundle constraints (slack)
    sigma_k = dual(1:K,1);     % dual variables of 
    
 %% C: Pricing Problem
  
    pricing_arc  = arcs; 
    pricing_arc.Edges.Weight = pricing_arc.Edges.Weight - pi_ij;
    
    sp.coliq       = zeros(A,K);
    
    for k = 1:K
        sp.path{k,i+1}      = shortestpath(pricing_arc,origin(k),destination(k));
        eidx                = findedge(pricing_arc, sp.path{k,i+1}(1:end-1), sp.path{k,i+1}(2:end));             
        sp.dist{k,i+1}      = sum(arcs.Edges.Weight(eidx));
        sp.coliq(eidx,k)      = demand(k);
    end
    
%% D: Determine column(s) to add. 
    
    % Determine for which commodity to add path. 
    
    col_idx   = find(cell2mat(sp.dist(:,i+1)) - sigma_k./demand' < 0)
    
    P             = P + size(col_idx,1); % total number of paths
    P_k(col_idx') = P_k(col_idx')+1;     % total number of paths for commodity k  
    
    
    c.ineq    = sp.coliq(:,col_idx); 
    
    
    sp.coleq  = c.eq(:,col_idx);
    
    c.eq2   =  sp.coleq;
    
    A1       = [ c.eq2 ; c.ineq];
    
    if isempty(col_idx) 
        dual_feasibility = 1; 
    end
    
    RMP.addCols(demand(col_idx)'.*cell2mat(sp.dist(col_idx,i+1)),A1);
    stop =  dual_feasibility + primal_feasibility;
    
    disp(['Stop: ',num2str(stop)]);   
    disp('-------------------------------------------------');

    
end
%    
