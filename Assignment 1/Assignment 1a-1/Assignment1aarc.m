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
      
    multiplepaths    = [];
%%  Postprocessing
    
    sol.DV   = reshape(cplex.Solution.x,nw.K,nw.A);
    sol.path = cell(nw.K,1);
    nw.name  = network.Nodes.Name';
    %nw.destination   = network.Nodes.Name(nw.destination);
    path1 = 1;
    
     for k = 1:nw.K
        sol.path{k,1}          = find(sol.DV(k,:));
        [sOut,tOut]            = findedge(network,sol.path{k,1});
        sol.path{k,2}          = [ sOut tOut];
        
        
        idx_origin    = find(sol.path{k,2}(:,1) == nw.origin(k));
        idx_dest      = find(sol.path{k,2}(:,2) == nw.destination(k));
        
        nopaths       = size(idx_origin,1);
        pathlength    = size(sol.path{k,2},1);
        
        
        nodes       = unique(sol.path{k,2});
        node_counts = histc(sol.path{k,2}(:), nodes);
        
        nopaths     = nopaths + size(find(node_counts >= 3),1);
        
        if nopaths > 1
            multiplepaths = [ multiplepaths k];
        end
        
        
        if nopaths == 1  
            
            for i = 1:(pathlength)
                if isequal(sol.path{k,2}(i,1),nw.origin(k)) == 1
                    sol.path{k,3}      = sol.path{k,2}(i,:);
                end
            end
       
           
            while sol.path{k,3}(end) ~= nw.destination(k)
                for p = 1:pathlength
                    if isequal(sol.path{k,3}(end),sol.path{k,2}(p,1)) == 1
                        sol.path{k,3}(end+1) = sol.path{k,2}(p,2);
                    end
                end
            end
            
            
            sol.path{k,4} = nw.name(sol.path{k,3});
            pathname = [sol.path{k,4}{1,1}];
            j =1;
            while j < pathlength+1
                pathname = [ pathname, '-', sol.path{k,4}{1,j+1} ]; 
                j = j+1;
            end
          
            
            sol.p{path1,1} = k;
            sol.p{path1,2} = network.Nodes.Name(nw.origin(k));
            sol.p{path1,3} = network.Nodes.Name(nw.destination(k));
            sol.p{path1,4} = pathname;
            sol.p{path1,5} = sol.DV(k,sol.path{k,1}(1,1)); 
            sol.p{path1,6} = nw.demand(k); 
            sol.p{path1,7} = sum(network.Edges.Weight(sol.path{k,1}),1);
            path1 = path1+1;
        end
     end
     
     T = cell2table(sol.p,...
    'VariableNames',{'Commodity' 'Origin' 'Destination' 'Path' 'Quantity' 'Demand' 'Cost'});
    writetable(T,'onepath.txt')
    
    
    for j = 1:size(multiplepaths,2)
        sol.mpath{j,1} = network.Edges(sol.path{multiplepaths(j),1},1:2);
        sol.mpath{j,1}.Quantity = sol.DV(multiplepaths(j),sol.path{multiplepaths(j),1})';
        sol.mpath{j,2} = nw.name(nw.origin(multiplepaths(j)));
        sol.mpath{j,3} = nw.name(nw.destination(multiplepaths(j)));
        sol.mpath{j,4} = nw.demand(multiplepaths(j));
        
    end
    %nw.name(nw.destination(multiplepaths))'
    
    
%% Functions 
% To return index of decision variables

function out = Aindex(K,k,a)
    out =  (a-1)*K + k;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
end

