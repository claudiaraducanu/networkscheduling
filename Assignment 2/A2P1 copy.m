%%  Initialization
% Claudia Raducanu and Luka Van de Sype
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64'); %Luka
% clc
 clearvars
% clear all

%%  Determine input
%%
input      =  'Assignment2.xlsx';

% Inputs
[P, R, L, fare, fare_r, demand, col, delta, Q, farecost, Bpr, ... 
          recap_p, recap_r, recaprate] = setup2P1(input);
[AC,B,timespace] = read_schedule(input);

% More variables
K = 4;      % set of fleet types
Lf = 208;   % set of aircraft flights (not going by bus)
Gk = size(timespace(1).ga,1) + size(timespace(2).ga,1) ... % total amount 
     + size(timespace(3).ga,1) + size(timespace(4).ga,1);  % of ground arcs (1364)
 
% Loop stoppers
not_opt_col = 1;
not_opt_row = 1;
colz = [];
rowz = [];

%% The costs
%%
% Cost of an aircraft being on a ground arc ( = 0)
y_cost = zeros(Gk,1);

% Cost of assigning ac type k to flight i (per K all L are given)
f_cost = [timespace(1).fl.Cost; timespace(2).fl.Cost; ...
          timespace(3).fl.Cost; timespace(4).fl.Cost];

% Cost of spilled passengers 
t_cost = zeros(numel(col(:,1)),1);
for p = 1:numel(col(:,1))    
    t_cost(p) = farecost(col(p,1),col(p,2));  % select only costs you need
end


 %%  Initiate CPLEX model
 %%
        %   Create model 
        model                 =   'IFAM';        % name of model
        RMP                   =    Cplex(model); % define the new model
        RMP.Model.sense       =   'minimize';

        %   Add the columns
        DV                      = Gk ...              % DV for FAM (y_ak)
                                  + K*Lf  ...         % DV for FAM (f_ik)
                                  + numel(col(:,1));  % DV for PMF part

        obj                     =   [y_cost; f_cost; t_cost] ; % y_ak, f_ik, t_pr
        lb                      =   zeros(DV, 1);   % Lower bounds
        ub                      =   inf(DV, 1);     % Upper bounds             
        
        RMP.addCols(obj, [], lb, ub);
        
   %%  Constraints
   %%
    % 1. Capacity constraint (for actual flight (without 'busflights')
        for i = 1:Lf
            C11 = zeros(1,DV); % Spilled pax from p
            C12 = zeros(1,DV); % Recaptured pax from r to p
            C13 = zeros(1,DV); % Capacity for each flight leg
            for pr = 1:DV
                if delta{col(pr,1),1}(i) ~= 0 
                    C11(pr)= 1;
                end
                if delta{col(pr,2),1}(i) ~= 0 && ismember(col(pr,2), col(:,1))
                    C12(pr) = Bpr(col(pr,1),col(pr,2)); 
                end
            end
            for k = 1:K
                C13(k) = AC.Seats(k);
            end
            C1 = C13 + C11 - C12;
            RMP.addRows(Q(i), C1, inf, sprintf('Capacity_%d',i));
        end
    % 2. Capacity constraint (for only 'busflights')
        for i = Lf+1:L
            C21 = zeros(1,DV); % Spilled pax from p
            C22 = zeros(1,DV); % Recaptured pax from r to p
            for pr = 1:DV
                if delta{col(pr,1),1}(i) ~= 0 
                    C21(pr)= 1;
                end
                if delta{col(pr,2),1}(i) ~= 0 && ismember(col(pr,2), col(:,1))
                    C22(pr) = Bpr(col(pr,1),col(pr,2)); 
                end
            end
            C2 = 4*54 + C21 - C22; % 4 busses with 54 seats for each flight
            RMP.addRows(Q(i), C2, inf, sprintf('Capacity_%d',i));
        end
        
    % 3. Each flight is operated by 1 aircraft type
        for i = 1:Lf
            C3 = zeros(1,DV);
            for k = 1:K
                C3(k) = 1;
            end
            RMP.addRows(1, C3, 1, sprintf('Flights_%d',i));
        end
        
    % 4. What goes in, goes out BY CLAUDIA
        for k = 1:K
            N_k = size(timespace(k).node,1);
            for n = 1:N_k
               [idx_ofl,idx_ifl, idx_on, idx_in] = fik_idx(timespace,k,n);
               for ii = 1:size(idx_ofl,1)
                   C4(Findex(k,idx_ofl(ii,1)) =  1;
               end
               for ii = 1:size(idx_ifl,1)
                   C4(Findex(k,idx_ifl(ii,1)) = -1;
               end            
               
            end
        end
            
    % 5. Fleet size is not exceeded
        for k = 1:K
           C5 = zeros(1,DV);
           for a = 1:numel(timespace(k).nga.Loc) %set of overnight ground arcs
               C5(a) = 1;
           end
           RMP.addRows(0, C5, AC.Units(k), sprintf('Fleetsize_%d',i));
        end

    %%  Execute model 
    %%
    sol = RMP.solve();
    OV = [];
    OV = [OV;RMP.Solution.objval];
    primal  = RMP.Solution.x;   
    dual    = RMP.Solution.dual;
    RMP.writeModel([model '.lp']);
    
    Alist = [RMP.Model.A];

%% SUPER BIG LOOP
%%
titer =0;
while not_opt_col == 1 || not_opt_row == 1

    titer = titer +1;
    disp('-------------------------------------------------');
    disp(['Iteration: ',num2str(titer)]);   
    disp('-------------------------------------------------');
    
    dual    = RMP.Solution.dual;
    
    
%% COLUMN GENERATION
%%
iter=0;
   while not_opt_col == 1
   %% Iterations
   %%
   iter = iter +1;
    disp('-------------------------------------------------');
    disp(['Column iteration: ',num2str(iter)]);   
    disp('-------------------------------------------------');
    
   %% Pricing Problem
   %%
   % Get dual variables 
    pi   = dual(1:L);                 % dual variables of capacity constraints (slack)
    sigma = zeros(P,1);
    if size(dual) > L
        for ro = 1:numel(rowz)
            sigma(rowz(ro)) = dual(L+ro);
        end
    end
    
    checkcol = 0; 
    R = P;
    
    for re = 1:numel(recap_p) % run only for those that can be recaptured
            pi_j = 0;
            pi_i = 0;
            
            for i = 1:L
                if delta{recap_p(re),1}(i) ~= 0     % search for the flight legs for each recap_p
                    pi_i = pi_i + pi(i);            % for those flight legs, add pi
                end
                if delta{recap_r(re),1}(i) ~= 0
                    pi_j = pi_j + pi(i);   
                end
            end
            
% recap_p is the path p that can be recaptured by path r in recap_r            
             if (fare(recap_p(re)) - pi_i) - (Bpr(recap_p(re),recap_r(re))*(fare(recap_r(re)) - pi_j )) - sigma(recap_p(re)) < 0 
                D = find(recap_r(re)==col(:,2));
                Dp = find(col(D,1)==recap_p(re));
                
                if isempty(Dp) == 1
                    checkcol = checkcol + 1; % to check if columns are added (loop stop)
                    
                    colz = [colz; [recap_p(re) recap_r(re)]]; % only the added columns
                    col = [col; [recap_p(re) recap_r(re)]];   % all columns
                    
                    t_cost = [t_cost;farecost(col(recap_p(re),1),col(recap_p(re),2))];
                    
                    % New objective function
                    A       = RMP.Model.A(:,recap_p(re));
                    
                    % add a constraint to every L if the leg is used in p.
                    deltasp = 0;
                    deltasp = find(delta{recap_p(re),1}(:)~=0); % find 
                    if isempty(deltasp) == 0        % if the list is not 0, add a constraint
                         A(deltasp)= 1;
                    end
                    
                    % substract rr to every L if the leg is used in r.
                    deltasr = 0;
                    deltasr = find(delta{recap_r(re),1}(:)~=0);
                    if isempty(deltasr) == 0        % if the list is not 0, add a constraint
                        A(deltasr)= A(deltasr) - Bpr(recap_p(re),recap_r(re)); 
                    end
                    
                    %Alist   = [Alist,A];
                    
                    obj     = farecost(recap_p(re),recap_r(re));
                    lb      = [0];
                    ub      = [inf];
                    %ctype   = []; 

                    RMP.addCols(obj,A,lb,ub);

                end
            end
    end
    
     
    RMP.solve();
    OV = [OV;RMP.Solution.objval];
    RMP.writeModel([model '.lp']);
    
    dual = RMP.Solution.dual;
    pi   = dual(1:L);                 % dual variables of capacity constraints (slack)
    sigma = zeros(P,1);
    if size(dual) > L
        for ro = 1:numel(rowz)
            sigma(rowz(ro))  = dual(L+ro);
        end
    end

    
    
    if checkcol == 0
        not_opt_col = 0;
    else 
        not_opt_col = 1;
        not_opt_row = 1;       
    end    
    
   end
     
 
%% ROW GENERATION
%%
     iteration = 0;

   while not_opt_row == 1
  %%  Initiate CPLEX model
  iteration = iteration +1;

    disp('-------------------------------------------------');
    disp(['Row iteration: ',num2str(iteration)]);   
    disp('-------------------------------------------------');
    size(rowz)
  
  DV                      = Gk + K*Lf ...       % DV for FAM part
                            + numel(col(:,1));  % DV for PMF part (incl new)
  primal = RMP.Solution.x;
     
        %% Separation Problem
    % 2. Demand constraint
    check = 0;
    for p = 1:P
        B = find(col(:,1)==p);
        
        if sum(primal(B)) > demand(p)
            D = find(p==rowz);
            
            if isempty(D) == 1
                check = check + 1;
                C = zeros(1,DV);
                C(B) = 1;
                RMP.addRows(-inf, C, demand(p), sprintf('Demand_%03d',p));
                rowz = [rowz;p];
            end
        end
    end 

    
           %%  Execute model 
    %   Run CPLEX
    RMP.solve();
    OV = [OV;RMP.Solution.objval];

    
%     %   Get dual variables
%     primal  = RMP.Solution.x;   
%     dual    = RMP.Solution.dual;
    
    
    if check == 0
        not_opt_row = 0;
    else
        not_opt_row = 1;
        not_opt_col = 1;
    end       
        
   end
    
end

%% Post Processing
%%
function out = Yindex(k, a)
    NGA = 34;
    out =  (k - 1) * NGA + a;  % Function given the variable index for each Y(k,a)
end

function out = Findex(k, i)
    LF = 208;
    out = 1364 + (k - 1) * LF + i;  % Function given the variable index for each Y(k,a)
end

function out = Tindex(i)
    LF = 208;
    out = 1364 + (k - 1) * LF + i;  % Function given the variable index for each Y(k,a)
end















