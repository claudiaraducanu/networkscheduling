%%  Initialization
% Claudia Raducanu and Luka Van de Sype
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
%   Select input file and sheet

input      =  'Input_AE4424_Ass1Verification.xlsx';
%filsol      =   'Solutions.xlsx';

%% Inputs

    [Nodes, K, cost, capacity, origin, destination, demand,s,t] = ...
    read_arc_v(input);

    A = size(s,1);
    capacity = nonzeros(capacity);
    
%% Initialization 
%Parameters
    
    %Set of shortest path for commodity k
    sp.dist = cell(K,1);
    sp.path = cell(K,1);
    sp.pred = cell(K,1);
    
    for k = 1:K
        [sp.dist{k,1}, sp.path{k,1}, sp.pred{k,1}] = graphshortestpath(sparse(cost),origin(k),destination(k),...
                        'Directed', true);
    end
    
    sp.arcs   = mat2cell([s t],ones(size(s,1),1));
   
   % If arc belongs to path 
    
    d_ij_p     = zeros(K,A);
    
     for k=1:K
         for a = 1:A
            pathlength = size(sp.path{k,1},2);
            if pathlength > 2 
                ainp = repelem(sp.path{k,1},2);
                ainp = ainp(2:end-1);
                ainp = transpose(reshape(ainp,2,pathlength-1));
                for ac = 1:(pathlength-1) 
                    if sum(ainp(ac,:) == sp.arcs{a,1}) == 2
                        d_ij_p(k,a)  = 1;
                    end
                end
            else
                if sum(sp.path{k,1} == sp.arcs{a,1}) == 2
                     d_ij_p(k,a)  = 1;
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
    DV                      =  K+A;  % Number of Decision Var (xijk)

%%  Objective Function
        
    c_p                     =   reshape(transpose(demand).*cell2mat(sp.dist),K,1);
    M                       =   ones(A,1)*1000;
    obj                     =   [c_p; M] ;
    lb                      =   zeros(DV, 1);                                 %Lower bounds
    ub                      =   inf(DV, 1);                                   %Upper bounds
    %ctype                   =   char(ones(1, (DV)) * ('I'));                  %Variable types 'C'=continuous; 'I'=integer; 'B'=binary

l = 1;
for k = 1:K
    %for p = 1:p(k)
        NameDV (l,:)  = ['F_' num2str(k,'%02d')];
        l = l + 1;
    %end
end 
for a = 1:A
    %for p = 1:p(k)
        NameDV (l,:)  = ['S_' num2str(a,'%02d')];
        l = l + 1;
    %end
end 
% cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
    RMP.addCols(obj, [], lb, ub);

%%  Constraints
% 1. Commodity constraint
    for k = 1:K
        C1 = zeros(1,DV);
        C1(Findex(k)) = 1;   
        RMP.addRows(1, C1, 1,sprintf('Commodity_%d',k));
    end
% 
% 2. Bundle constraint
    for a = 1:A
        C2 = zeros(1,DV);
        for k = 1:K
            C2(Findex(k)) = demand(k)*d_ij_p(k,a);
        end
        C2(Sindex(a,K)) = -1; 
        RMP.addRows(0, C2, capacity(a),sprintf('Bundle_%d',a));
    end

%%  Execute model
    RMP.Param.mip.limits.nodes.Cur    = 1e+8;        %max number of nodes to be visited (kind of max iterations)
    RMP.Param.timelimit.Cur           = 120;         %max time in seconds
% %   Run CPLEX
    
    RMP.solve();
    RMP.writeModel([model '.lp']);
    
    pi_ij = RMP.Solution.dual; % dual variables of bundle constraints
          =  RMP.Solution.x;
    
    
    
%Function to return index of decision variables
function out = Findex(k)
    out =  k;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
end
function out = Sindex(a,K)
    out = K + a;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
end



