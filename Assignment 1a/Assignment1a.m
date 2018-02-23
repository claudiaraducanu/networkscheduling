%%  Initialization
%addpath('D:Documents\cplex\matlab\x64_win64'); %Luka
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Satrio
%addpath('Users\Luka\Documents\IBM\ILOG\CPLEX_Studio1271'); %Bryan
clc
clearvars
close all
warning('off','MATLAB:lang:badlyScopedReturnValue')
warning('off','MATLAB:xlswrite:NoCOMServer')


%%  Determine input
%   Select input file and sheet
input        =   [pwd '/Input_AE4424_Ass1P1.xlsx'];
filsol      =   'Solutions.xlsx';

cost        = xlsread(input,1,'D2:D31');
capa        = xlsread(input,1,'E2:E31');
quant       = xlsread(input,2,'D2:D41');

origin      = xlsread(input,1,'D2:D41'); % origin of comodity k
destination = xlsread(input,1,'D2:D41'); % destination of comodity k

%%  Define formulations
%   The Sets
Nodes = 16;     % set of airports     
K =40;          % set of commodities

%%  Initiate CPLEX model
%   Create model 
model                   =   'MCF_Model';  % name of model
cplex                   =   Cplex(model); % define the new model
cplex.Model.sense       =   'minimize';
%   Decision variables
DV                    =  Nodes*Nodes*K;  % Number of Decision Var (xijk)


%%  Objective Function
 cost_OF       =   reshape(cost, Nodes*Nodes*K, 1);

obj                     =   [cost_OF];
lb                      =   zeros(DV, 1);                                 %Lower bounds
ub                      =   inf(DV, 1);                                   %Upper bounds
ctype                   =   char(ones(1, (DV)) * ('I'));                  %Variable types 'C'=continuous; 'I'=integer; 'B'=binary


l = 1;                                      % Array with DV names
for i = 1:Nodes
    for j = 1:Nodes                     % of the x_{ij}^k variables
        for k = 1:K
            NameDV (l,:)  = ['X_' num2str(i,'%02d') ',' num2str(j,'%02d') '_' num2str(0,'%02d')];
            l = l + 1;
        end
    end
end


% cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
cplex.addCols(obj, [], lb, ub, ctype, NameDV);


%%  Constraints
% 1. Demand Verification (#pax <= demand from i to j)
a  = zeros(K);
for i = 1:Nodes
    for k = 1:K
        for j = 1:Nodes
            C11 = zeros(1,DV);
            C12 = zeros(1,DV);
            C11(Xindex(i,j,k)) = 1;
            C12(Xindex(j,i,k)) = 1;
        end
        if i = origin(k)       % i is element of origin of k
            a(k) = quant(k);
        elseif i = destination(k)  % i is element of destination of k
            a(k) = -quant(k);
        else
            a(k) = 0;
        end
        C1 = C11-C12    
        cplex.addRows(a(k), C1, a(k),sprintf('Direct_Demand_Constraint_%d_%d_%d',i,j,k));
    end
end

% 2. Capacity constraint
for i = 1:Nodes
    for j = 1:Nodes
        for k = 1:K
        C2 = zeros(1,DV);
        C2(Xindex) = 1
        end
        cplex.addRows(0), C2, capa(i,j),sprintf('Direct_Demand_Constraint_%d_%d_%d',i,j,k));
    end
end




%%  Execute model
cplex.Param.mip.limits.nodes.Cur    = 1e+8;        %max number of nodes to be visited (kind of max iterations)
cplex.Param.timelimit.Cur           = 120;         %max time in seconds
%   Run CPLEX
cplex.solve();
cplex.writeModel([model '.lp']);

sol.profit = cplex.Solution.objval;
solution_DV = cplex.Solution.x;

%%  Postprocessing
%Store direct results
   solution_x_ijk = round(reshape(solution_DV(1:1:Nodes*Nodes*K),Nodes,Nodes,K));

for k = 1:K
    xlswrite(filsol,solution_x_ij',k,'C3:V22')
end


%end

%Function to return index of decision variables
function out = Xindex(m,n,p)
    Nodes = 16;
    K = 40;
    out =  (m - 1)*Nodes*K + (n-1)*K + p;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
end



