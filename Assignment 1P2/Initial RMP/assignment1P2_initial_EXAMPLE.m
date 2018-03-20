%%  Initialization
% Claudia Raducanu and Luka Van de Sype
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka

clearvars
clear all

%%  Determine input
%   Select input file and sheet

input      =  'Input_Example.xlsx';

%% Inputs
[P, L, fare_p, demand , capacity, pathflights, flightnrs] = matrixsetup1P2_initial_EXAMPLE(input) ; 

% The binary value: Delta
delta = cell(P,1); % For each Path all Flights are checked: delta{p,1}(i)
for p = 1:P
    delta{p,1} = zeros(L,1);
    for i = 1:L
        if pathflights(p,1) == flightnrs(i) || pathflights(p,2) == flightnrs(i)
            delta{p,1}(i) = 1;            
        end
    end
end

% The daily unconstrained demand on flight(i): Q
Q = zeros(L,1);
for i = 1:L
    a = zeros(P,1);
    for p = 1:P
        a(p) = demand(p)*delta{p,1}(i);
    end
    Q(i) = sum(a);
end

  %%  Initiate CPLEX model
        %   Create model 
        model                 =   'Initial';  % name of model
        RMP                   =    Cplex(model); % define the new model
        RMP.Model.sense       =   'minimize';

        %   Decision variables
        DV                      =  P;  % Number of Decision Var (xijk)
   %% Objective function

        obj                     =   [fare_p] ;
        lb                      =   zeros(DV, 1);                                 %Lower bounds
        ub                      =   inf(DV, 1);                                   %Upper bounds              

        RMP.addCols(obj, [], lb, ub);
        
   %%  Constraints
    % 1. Capacity constraint
        for i = 1:L
            C1 = zeros(1,DV);
            for p = 1:P
                if delta{p,1}(i) ~= 0
                    C1(Tindex(p)) = 1;
                end
            end
            RMP.addRows(Q(i)-capacity(i), C1, inf, sprintf('Capacity_%d',i));
        end
       %%  Execute model
    RMP.Param.mip.limits.nodes.Cur    = 1e+8;        %max number of nodes to be visited (kind of max iterations)
    RMP.Param.timelimit.Cur           = 120;         %max time in seconds
 
    %   Run CPLEX
    RMP.solve();
    RMP.writeModel([model '.lp']);
    
        %%  Function to return index of decision variables
function out = Tindex(p) 
        out = p;
end
    
    
    
    
    
