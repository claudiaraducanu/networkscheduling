%%  Initialization
% Claudia Raducanu and Luka Van de Sype
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
%   Select input file and sheet

%input     =   'Input_AE4424_Ass1Verification.xlsx';
input      =   'Input_AE4424_Ass1P1.xlsx';

%% Inputs

[nw.N, nw.K, nw.cost , ...
 nw.capacity , nw.origin , nw.destination, nw.demand, nodes] ...
    = matrixsetup(input,1);

%%  Initiate CPLEX model
%   Create model 

model                   =   'ARC_Model';    % name of model
cplex                   =    Cplex(model);  % define the new model
cplex.Model.sense       =   'minimize';

%   Decision variables
DV                      =  nw.N*nw.N*nw.K;  % Number of Decision Variables (x_i_j_k)

%%  Objective Function
cost_OF       =   reshape(nw.cost', nw.N*nw.N,1);

% As the cost is not dependent on k, but it should be Nodes*Nodes*K long,
% each element is repeated K times.

cost_OF       =   repelem(cost_OF,nw.K);        

% Initialize the objective function column
obj                     =   cost_OF ;
lb                      =   zeros(DV, 1);                                 %Lower bounds
ub                      =   inf(DV, 1);                                   %Upper bounds
ctype                   =   char(ones(1, (DV)) * ('I'));                  %Variable types 'C'=continuous; 'I'=integer; 'B'=binary

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
            C1(Xindex(i,j,k,nw.N,nw.K)) =  1;
            C1(Xindex(j,i,k,nw.N,nw.K)) = -1;
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
for i = 1:nw.N
    for j = 1:nw.N
        C2 = zeros(1,DV);
        for k = 1:nw.K
        C2(Xindex(i,j,k,nw.N,nw.K)) = 1;
        end
        cplex.addRows(0, C2, nw.capacity(i,j),sprintf('Capacity_Constraint_%d_%d',i,j));
    end
end

%%  Execute model

    cplex.Param.mip.limits.nodes.Cur    = 1e+8;        %max number of nodes to be visited (kind of max iterations)
    cplex.Param.timelimit.Cur           = 120;         %max time in seconds
    
%   Run CPLEX

    cplex.solve();
    cplex.writeModel([model '.lp']);

%   Get solution

    sol.cost         = cplex.Solution.objval;
    sol.DV           = cplex.Solution.x;
    
%%  Postprocessing

    [sol.x(:,1),sol.x(:,2),sol.x(:,3)] =  ...
        DVindex(find(sol.DV),nw.N,nw.K);
    sol.x(:,4) = sol.DV(find(sol.DV)); %#ok<FNDSB>
    sol.x      = sortrows(sol.x,3);
  
    sol.xname(:,1)  = nodes.Name(sol.x(:,1));
    sol.xname(:,2)  = nodes.Name(sol.x(:,2));
    sol.xname(:,[3 4])  = num2cell(sol.x(:,[3,4]));
    
    
     for k = 1:nw.K
        
        arcs_k     = sol.x(sol.x(:,3) == k,1:2);     % store arcs used by commodity
        pathlength = size(arcs_k,1);
        path       = reshape(transpose(arcs_k),1,pathlength*2);
        
        idx_o     = find(path == nw.origin(k));
        idx_d     = find(path == nw.destination(k));
        
        path_ord   = zeros(size(path));
        
        path_ord(1:2) = path(idx_o:idx_o+1);
        path_ord(end-1:end) = path(idx_d-1:idx_d);
        
       
        

%         if size(arcs_k,1) == 1
%             sol.paths{k,1} = nodes.Name(arcs_k(1:2));
%         else
%         % Initiate the path for commodity with the origin and next
%         % airport
%             path_k = arcs_k( arcs_k(:,1) == nw.origin(k),[1 2]);
%         %Determine the number of paths possible based on origin
%             num_path_k = size(path_k,1);
        
     end
    
%% Functions 
% To return index of decision variables

function out = Xindex(m,n,p,nd,k)
    out =  (m - 1)*nd*k + (n-1)*k + p;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
end

% To return indexing of DV in math model from CPLEX index

function [i,j,k] = DVindex(indx,nd,K)
    i  = floor((indx-1)/(nd*K))+1;
    j  = floor(((indx-(i-1)*nd*K)-1)/K)+1;
    k  = mod(indx,K);
    mm =  k == 0;
    k(mm) = K;  % Function given the DV index gives X(i,j,k) 
end
