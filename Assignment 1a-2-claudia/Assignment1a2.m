%%  Initialization
% Claudia Raducanu and Luka Van de Sype
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
%   Select input file and sheet

input      =   [pwd '/Input_AE4424_Ass1P1.xlsx'];
%filsol      =   'Solutions.xlsx';

%% Inputs

[Nodes, K, cost, capacity, origin, destination, demand, airports] = ...
    matrixsetup(input);

%% Initialization 
%Parameters
    
    %Set of shortest path for commodity k
    
    sp.dist = cell(40,1);
    sp.path = cell(40,1);
    sp.pred = cell(40,1);
    
    for k = 1:K
        [sp.dist{k,1}, sp.path{k,1}, sp.pred{k,1}] = graphshortestpath(sparse(cost),origin(k),destination(k),...
                        'Directed', false);
    end
    
    
    
%     for k = 1:K
%          p_k_init{k,1} = graphshortestpath(costgraph,origin(k),destination(k));
%          numP(k) = numel(p_k_init{k,1}); % number of airports used in each path
%     end

%P = numel(p_k_init);

% cost of transporting one unit of k on path p
% costp = zeros(k,1) + 1000 ; % set default value: cost=1000
% for k = 1:K
%     for n = 1:(numP(k)-1)
%         cost_part = zeros(n-1,1);
%         cost_part(n) = cost(p_k_init{k,1}(n),p_k_init{k,1}(n+1)); 
%     end
%     costp(p_k_init{k,1}) = sum(cost_part);
% end

%% Objective of OF
% cost of transporting total of k on path p
% costtotal = costp.*demand';




%% Not yet

% % %%  Initiate CPLEX model
% % %   Create model 
% % model                   =   'MCF_Model';  % name of model
% % cplex                   =   Cplex(model); % define the new model
% % cplex.Model.sense       =   'minimize';
% % %   Decision variables
% % DV                      =  K*K;  % Number of Decision Var (xijk)
% % 
% % 
% % %%  Objective Function
% %         
% % 
% % obj                     =   costtotal ;
% % lb                      =   zeros(DV, 1);                                 %Lower bounds
% % ub                      =   inf(DV, 1);                                   %Upper bounds
% % ctype                   =   char(ones(1, (DV)) * ('I'));                  %Variable types 'C'=continuous; 'I'=integer; 'B'=binary
% % 
% % 
% % l = 1;                          % Array with DV names  
% % for p = 1:P
% %     NameDV (l,:)  = ['F_' num2str(p,'%02d')];
% %     l = l + 1;
% % end
% % 
% % 
% % % cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
% % cplex.addCols(obj, [], lb, ub, ctype, NameDV);
% % 
% % 
% % %%  Constraints
% % % 1. Demand Verification (#pax <= demand from i to j)
% % for i = 1:Nodes
% %     for k = 1:K
% %         C1 = zeros(1,DV);
% %         for j = 1:Nodes
% %             C1(Xindex(i,j,k)) =  1;
% %             C1(Xindex(j,i,k)) = -1;
% %         end
% %         if i == origin(k)       % i is element of origin of k
% %             cplex.addRows(demand(k), C1, demand(k), ...
% %                 sprintf('Quantity_Constraint_%d_%d',i,k));
% %         elseif i == destination(k)  % i is element of destination of k
% %             cplex.addRows(-demand(k), C1, -demand(k), ...
% %                 sprintf('Quantity_Constraint_%d_%d',i,k));
% %         else
% %             cplex.addRows(0, C1, 0, ...
% %                 sprintf('Quantity_Constraint_%d_%d',i,k));
% %         end    
% %     end
% % end
% % 
% % % 2. Capacity constraint
% % for i = 1:Nodes
% %     for j = 1:Nodes
% %         C2 = zeros(1,DV);
% %         for k = 1:K
% %         C2(Xindex(i,j,k)) = 1;
% %         end
% %         cplex.addRows(0, C2, capa(i,j),sprintf('Capacity_Constraint_%d_%d',i,j));
% %     end
% % end
% % 
% % %%  Execute model
% % cplex.Param.mip.limits.nodes.Cur    = 1e+8;        %max number of nodes to be visited (kind of max iterations)
% % cplex.Param.timelimit.Cur           = 120;         %max time in seconds
% % %   Run CPLEX
% % cplex.solve();
% % cplex.writeModel([model '.lp']);
% % 
% % sol.cost = cplex.Solution.objval;
% % solution_DV = cplex.Solution.x;
% % 
% % %%  Postprocessing
% % %Store direct results
% % %    solution_x_ijk = round(reshape(solution_DV(1:1:Nodes*Nodes*K),Nodes,Nodes,K));
% % % 
% % % for k = 1:K
% % %     xlswrite(filsol,solution_x_ij,k,'C3:V22')
% % % end
% % 
% % 
% % %end
% % 
% % %Function to return index of decision variables
% % function out = Xindex(m,n,p)
% %     Nodes = 16;
% %     K = 40;
% %     out =  (m - 1)*Nodes*K + (n-1)*K + p;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
% % end
% % 
% % 
% % 