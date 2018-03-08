%%  Initialization
% Claudia Raducanu and Luka Van de Sype

addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka


%%  Determine input
%   Select input file and sheet

inputA      =   [pwd '/Input_AE4424_Ass1P1.xlsx'];
%filsol      =   'Solutions.xlsx';

%% Inputs
% These are the ouputs of the function matrixsetup;
%     Nodes = set of airports
%     K     = set of commodities
%     cost  = Nodes*Nodes matrix of arc costs
%     capa  = Nodes*Nodes matrix of arc capacities
%     quant = 1*K row vector of demand
%     origin = 1*K row vector of each commodity origin
%     destination = 1*K row vector of each commodity destination

[Nodes, K, cost , capa , origin , destination, demand, costgraph] = matrixsetup(inputA,1);

%% Parameters
%Set of shortest path for each commodity k
Pshort = {}; % P{3,1}(2) is the second step in the path of k=3
P = zeros(K,1);
 for k = 1:K
     Pshort{k,1} = shortestpath(costgraph,origin(k),destination(k));
     numP(k) = numel(Pshort{k,1}); % number of airports used in each path
     P(k) = 1;
 end

% cost of transporting one unit of k on path p
costp = zeros(k,1) + 1000 ; % set default value: cost=1000
for k = 1:K
    for n = 1:(numP(k)-1)
        cost_part = zeros(n-1,1);
        cost_part(n) = cost(Pshort{k,1}(n),Pshort{k,1}(n+1)); 
    end
    costp(Pshort{k,1}) = sum(cost_part);
end

% %% Objective of OF
% % cost of transporting total of k on path p
% costtotal = costp.*demand';
% 
% 
% %% Not yet
% 
% %%  Initiate CPLEX model
% %   Create model 
% model                   =   'MCF_Model';  % name of model
% cplex                   =   Cplex(model); % define the new model
% cplex.Model.sense       =   'minimize';
% %   Decision variables
% DV                      =  K*A;  % Number of Decision Var (xijk)
% 
% 
% %%  Objective Function
%         
% obj                     =   [costtotal; M] ;                              %Cost and the big M method
% lb                      =   zeros(DV, 1);                                 %Lower bounds
% ub                      =   inf(DV, 1);                                   %Upper bounds
% ctype                   =   char(ones(1, (DV)) * ('I'));                  %Variable types 'C'=continuous; 'I'=integer; 'B'=binary
% 
% 
% l = 1;                          % Array with DV names  
% for k = 1:K
%     for p = 1:p(k)
%         NameDV (l,:)  = ['F_' num2str(p,'%02d')];
%         l = l + 1;
%     end
% end
% for a = 1:A                     % should be for every arc
%         NameDV (l,:)  = ['S_' num2str(p,'%02d')];
%         l = l + 1;
% end
% 
% % cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
% cplex.addCols(obj, [], lb, ub, ctype, NameDV);
% 
% 
% %%  Constraints
% % 0. Slack and Delta
% for a = 1:A
%     slack(a) = max(0,
% end
% 
% % 1. Demand Verification (#pax <= demand from i to j)
% for k = 1:K
%     C1 = zeros(1,DV);
%     for p = 1:P(k)
%         C1(Findex(k,p)) = 1    
%     end
%     cplex.addRows(1, C1, 1,sprintf('Demand_Verification_%d_%d',k,p));
% end
% 
% % 2. Capacity constraint
% for a = 1:A
%     C2 = zeros(1,DV);
%     for k = 1:K
%         for p = 1:P(k)
%             C2(Findex(k,p)) = demand(k)*delta(a,k,p);
%         end
%     end
%     if C2 - capa(a) >= 0
%         slack(a) = C2 - capa(a);
%     else
%         slack(a) = 0;
%     end
%     cplex.addRows(0, C2-slack(a), capa(a),sprintf('Capacity_Constraint_%d_%d',i,j));
% end
% 
% %%  Execute model
% cplex.Param.mip.limits.nodes.Cur    = 1e+8;        %max number of nodes to be visited (kind of max iterations)
% cplex.Param.timelimit.Cur           = 120;         %max time in seconds
% %   Run CPLEX
% cplex.solve();
% cplex.writeModel([model '.lp']);
% 
% sol.cost = cplex.Solution.objval;
% solution_DV = cplex.Solution.x;
% 
% %%  Postprocessing
% %Store direct results
% %    solution_x_ijk = round(reshape(solution_DV(1:1:Nodes*Nodes*K),Nodes,Nodes,K));
% % 
% % for k = 1:K
% %     xlswrite(filsol,solution_x_ij,k,'C3:V22')
% % end
% 
% 
% %end
% 
% %Function to return index of decision variables
% function out = Findex(p)
%     Nodes = 16;
%     K = 40;
%     out =  (m - 1)*Nodes*K + (n-1)*K + p;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
% end



