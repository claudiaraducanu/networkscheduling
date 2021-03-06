3%%  Initialization
% Claudia Raducanu and Luka Van de Sype
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
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
    solution = []; 
    dual_sol = [];
    no_col   = [];
    no_col2  = [];
%% B: Solve RMP
while i < 5
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
    RMP.Param.timelimit.Cur           = 120;         %max time in seconds
 
    %   Run CPLEX
    RMP.solve();
    RMP.writeModel([model '.lp']);
    
    
    solution(i) = RMP.Solution.objval;
    
    
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
    col_idx   = find(cell2mat(sp.dist(:,i+1)) < sigma_k./demand');
    
    columns{i,1}   = col_idx;
    no_col(i)        = size(col_idx,1);
    dual_sol(i)      = dual'*[ ones(K,1) ; capacity ];

%     for k = 1:K
%         for n = 1:(size(sp.path,2)-1)
%             col_exists     = ismember(sp.path{k,i+1},sp.path{k,n},'legacy');
%             sp.repeat(k,n) = sum(col_exists) == size(sp.path{k,i+1},2);
%         end     
%     end
%     
%     idx_old = (sum(sp.repeat(:,1:i),2) > 0) == col_idx;
%     col_idx = find(col_idx - idx_old)
    
    

    P             = P + size(col_idx,1); % total number of paths
    P_k(col_idx') = P_k(col_idx')+1;     % total number of paths for commodity k  
    
    
    c.ineq = [c.ineq sp.coliq(:,col_idx')]; 
    
    sp.coleq      = c.eq(:,col_idx');
    
    c.eq   = [ c.eq sp.coleq];
    
    if isempty(col_idx) 
        dual_feasibility = 1; 
    end
  
    obj.o = [ obj.o ; demand(col_idx)'.*cell2mat(sp.dist(col_idx,i+1))];
    
    
    stop =  dual_feasibility + primal_feasibility;
    
    
    disp(['Stop: ',num2str(stop)]);   
    disp('-------------------------------------------------'); 
    
end
 data    = cell(K,1); 
 
 table2 = arcs.Edges(find(pi_ij),1);
 table2.marginal = pi_ij(find(pi_ij));
 
 
 it1 = primal(1:K,1);
 com = find(it1);
 it1 = it1(com);
 
 for i = 1:size(com,1)
     data{com(i),1} = sp.path{com(i),1};
     frac{com(i),1} = it1(i);
 end
 
  it2 = primal(54:79); 
  com = columns{1,1}(find(it2),1);
  it2 = it2(it2>0);
 for i = 1:size(com,1)
     data{com(i),2} = sp.path{com(i),2};
     frac{com(i),2} = it2(i);
 end
 
  it3 = primal(80:103); 
  com = columns{2,1}(find(it3),1);
  it3 = it3(it3>0);
 for i = 1:size(com,1)
     data{com(i),3} = sp.path{com(i),3};
     frac{com(i),3} = it3(i);
 end
 
   it4 = primal(104:130);
   com = columns{3,1}(find(it4),1);
   it4 = it4(it4>0);
 for i = 1:size(com,1)
     data{com(i),4} = sp.path{com(i),4};
     frac{com(i),4} = it4(i);
 end
 
   it5 = primal(131:156);
   com = columns{4,1}(find(it5),1);
   it5 = it5(it5>0);
 for i = 1:size(com,1)
     data{com(i),5} = sp.path{com(i),5};
     frac{com(i),5} = it5(i);
 end

    paths = []; 
    for k= 1:K
        dim = size(find(double(~cellfun('isempty',data(k,:)))),2);
        paths = [ paths; k*ones(dim,1)];
    end  
    


  dd(:,1) = num2cell(paths);
  data    = data';
  dd(:,2) = data(~cellfun('isempty',data));  
 
  
  
  for j = 1:size(dd,1)
      pathlength = size(dd{j,2},2);
      pathname   = arcs.Nodes.Name(dd{j,2}(1,1));
      n = 1; 
      while n < pathlength
        pathname = strcat( pathname, '-', arcs.Nodes.Name(dd{j,2}(1,n+1)) ); 
        n = n+1;
      end
      dd{j,2} = pathname;
  end

  
  
  frac    = frac';
  
  dd(:,3) = frac(~cellfun('isempty',data));  
  for j = 1:size(dd,1)
      dd{j,4} = demand(dd{j,1})*dd{j,3};
  end
  
   T = cell2table(dd,...
    'VariableNames',{'Commodity' 'Path' 'Fraction' 'Quantity'});
    writetable(T,'columngen.txt');
  