%%  Initialization
% Claudia Raducanu and Luka Van de Sype
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
%   Select input file and sheet

input      =   'Input_AE4424_Ass1P1.xlsx';

%% Inputs

[nw.N, nw.K, nw.A , ...
 nw.origin , nw.destination, nw.demand, network] ...
    = matrixsetuparc(input,1);

nodes    = network.Nodes;
%%  Initiate CPLEX model
%   Create model 

model                   =   'ARC_Model';    % name of model
cplex                   =    Cplex(model);  % define the new model
cplex.Model.sense       =   'minimize';

%   Decision variables
DV                      =  nw.A*nw.K;  % Number of Decision Variables (x_i_j_k)

%%  Objective Function


obj                     =   repelem(network.Edges.Weight,nw.K);
lb                      =   zeros(DV, 1);                                 %Lower bounds
ub                      =   inf(DV, 1);                                   %Upper bounds

l = 1;                                      % Array with DV names
for i = 1:nw.N
    for j = 1:nw.N                     % of the x_{ij}^k variables
        for k = 1:nw.K
            NameDV (l,:)  = ['X_' num2str(i,'%02d') ',' num2str(j,'%02d') '_' num2str(k,'%02d')];
            l = l + 1;
        end
    end
end

% Set up objective function
cplex.addCols(obj, [], lb, ub); %ctype, NameDV);

%%  Constraints
% 1. Demand Verification (#pax <= demand from i to j)
for i = 1:nw.N
    for k = 1:nw.K
        C1 = zeros(1,DV);
        for j = 1:nw.N
              a1 = findedge(network,i,j);
              a2 = findedge(network,j,i);
              
              if a1 ~= 0
                  C1(Aindex(nw.K,k,a1)) =  1;
              end
              if a2 ~= 0
                  C1(Aindex(nw.K,k,a2)) =  -1;
              end

        end
        if i == nw.origin(k)       % i is element of origin of k
            cplex.addRows(nw.demand(k), C1, nw.demand(k), ...
                sprintf('Quantity_Constraint_%d_%d',i,k));
        elseif i == nw.destination(k)  % i is element of destination of k
            cplex.addRows(-nw.demand(k), C1, -nw.demand(k), ...
                sprintf('Quantity_Constraint_%d_%d',i,k));
        else
            cplex.addRows(0, C1, 0, ...
                sprintf('Quantity_Constraint_%d_%d',i,k));
        end    
    end
end

% 2. Capacity constraint
for a = 1:nw.A
    C2 = zeros(1,DV);
    for k = 1:nw.K
        C2(Aindex(nw.K,k,a)) = 1;
    end
    cplex.addRows(0, C2, network.Edges.capacity(a,1),sprintf('Capacity_Constraint_%d_%d',i,j));
end

%%  Execute model

%   Run CPLEX
    cplex.solve();
    cplex.writeModel([model '.lp']);

%   Get solution
    sol.cost         = cplex.Solution.objval;
      
%%  Postprocessing
    
    sol.DV   = reshape(cplex.Solution.x,nw.K,nw.A);
    sol.path = cell(nw.K,1);
    nw.origin = network.Nodes.Name(nw.origin);
    
    for k = 1:nw.K
        sol.path{k,1} = find(sol.DV(k,:));
        sol.path{k,2} = table2cell(network.Edges(sol.path{k,1},1));
        
        for i = 1:size(sol.path{k,1},2)
            if strcmp(sol.path{k,2}{i,1}{1,1},nw.origin(k,1)) == 1
                sol.path{k,3} = sol.path{k,2}{i,1};
            end
        end
        
        pathlength = size(sol.path{k,1},2); % Path length in terms of arcs
        
        
        if pathlength > 1
            i = 1; 
            while i < (pathlength+1)
                for j = 1:pathlength
                    if strcmp(sol.path{k,2}{j,1}{1,1},sol.path{k,3}{1,end}) == 1
                        sol.path{k,3}{1,end+1} = sol.path{k,2}{j,1}{1,2};
                    end
                end
                i = i+1;
            end
        end
%         
     end
    
    
    
%% Functions 
% To return index of decision variables

function out = Aindex(K,k,a)
    out =  (a-1)*K + k;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
end

