%%  Initialization
% Claudia Raducanu and Luka Van de Sype

addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka


%%  Determine input
%   Select input file and sheet

inputA        =   [pwd '/Input_AE4424_Ass1P1.xlsx'];
filsol      =   'Solutions.xlsx';

%% Inputs
% These are the ouputs of the function matrixsetup;
%     Nodes = set of airports
%     K     = set of commodities
%     cost  = Nodes*Nodes matrix of arc costs
%     capa  = Nodes*Nodes matrix of arc capacities
%     quant = 1*K row vector of demand
%     origin = 1*K row vector of each commodity origin
%     destination = 1*K row vector of each commodity destination

[Nodes, K, cost , capa , origin , destination, demand] = matrixsetup(inputA);


%%  Initiate CPLEX model
%   Create model 
model                   =   'MCF_Model';  % name of model
cplex                   =   Cplex(model); % define the new model
cplex.Model.sense       =   'minimize';
%   Decision variables
DV                      =  Nodes*Nodes*K;  % Number of Decision Var (xijk)


%%  Objective Function
cost_OF       =   reshape(cost, Nodes*Nodes, 1);

% As the cost is not dependent on k, but it should be Nodes*Nodes*K long,
% each element is repeated K times.
cost_OF       =   repelem(cost_OF,K);        

obj                     =   cost_OF ;
lb                      =   zeros(DV, 1);                                 %Lower bounds
ub                      =   inf(DV, 1);                                   %Upper bounds
ctype                   =   char(ones(1, (DV)) * ('I'));                  %Variable types 'C'=continuous; 'I'=integer; 'B'=binary


l = 1;                                      % Array with DV names
for i = 1:Nodes
    for j = 1:Nodes                     % of the x_{ij}^k variables
        for k = 1:K
            NameDV (l,:)  = ['X_' num2str(i,'%02d') ',' num2str(j,'%02d') '_' num2str(k,'%02d')];
            l = l + 1;
        end
    end
end


% cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
cplex.addCols(obj, [], lb, ub, ctype, NameDV);


%%  Constraints
% 1. Demand Verification (#pax <= demand from i to j)
for i = 1:Nodes
    for k = 1:K
        C1 = zeros(1,DV);
        for j = 1:Nodes
            C1(Xindex(i,j,k)) =  1;
            C1(Xindex(j,i,k)) = -1;
        end
        if i == origin(k)       % i is element of origin of k
            cplex.addRows(demand(k), C1, demand(k), ...
                sprintf('Quantity_Constraint_%d_%d',i,k));
        elseif i == destination(k)  % i is element of destination of k
            cplex.addRows(-demand(k), C1, -demand(k), ...
                sprintf('Quantity_Constraint_%d_%d',i,k));
        else
            cplex.addRows(0, C1, 0, ...
                sprintf('Quantity_Constraint_%d_%d',i,k));
        end    
    end
end

% 2. Capacity constraint
for i = 1:Nodes
    for j = 1:Nodes
        C2 = zeros(1,DV);
        for k = 1:K
        C2(Xindex(i,j,k)) = 1;
        end
        cplex.addRows(0, C2, capa(i,j),sprintf('Capacity_Constraint_%d_%d',i,j));
    end
end

%%  Execute model
cplex.Param.mip.limits.nodes.Cur    = 1e+8;        %max number of nodes to be visited (kind of max iterations)
cplex.Param.timelimit.Cur           = 120;         %max time in seconds
%   Run CPLEX
cplex.solve();
cplex.writeModel([model '.lp']);

sol.cost = cplex.Solution.objval;
solution_DV = cplex.Solution.x;

%%  Postprocessing
%Store direct results
%    solution_x_ijk = round(reshape(solution_DV(1:1:Nodes*Nodes*K),Nodes,Nodes,K));
% 
% for k = 1:K
%     xlswrite(filsol,solution_x_ij,k,'C3:V22')
% end


%end

%Function to return index of decision variables
function out = Xindex(m,n,p)
    Nodes = 16;
    K = 40;
    out =  (m - 1)*Nodes*K + (n-1)*K + p;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
end



