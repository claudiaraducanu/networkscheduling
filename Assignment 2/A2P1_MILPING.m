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
          recap_p, recap_r, recaprate,flightnrs] = setup2P1(input);
tic
[AC,B,timespace,count_time] = read_schedule(input);
time_start = toc; % after the time space network
find(timespace(1).fl.Departure > timespace(1).fl.Arrival);
% More variables
K  = size(AC,1);                % set of AIRCRAFT types
Lf = size(timespace(1).fl,1);   % set of aircraft flights (not going by bus)

GA  = zeros(1,size(timespace,2));
for k  = 1:size(timespace,2)
    GA(1,k) = size(timespace(k).gat,1);
end
Gk = sum(GA);




% Loop stoppers
not_opt_col = 1;
not_opt_row = 1;
colz = []; % list of added columns
rowz = []; % list of added rows
coliez = []; % list of the amount of added columns
rowiez = []; % list of the amount of added rows

%% The costs
%%
% Cost of an aircraft being on a ground arc ( = 0)
y_cost = zeros(Gk,1);

% Cost of assigning ac type k to flight i (per K all L are given)
f_cost =[]; 
for k = 1:K
    f_cost = [f_cost; timespace(k).fl.Cost;];
end

% Cost of spilled passengers 
t_cost = zeros(numel(col(:,1)),1);
for p = 1:numel(col(:,1))    
    t_cost(p) = farecost(col(p,1),col(p,2));  % select only costs you need
end

time_pre = toc;

 %%  Initiate CPLEX model
 %%
        %   Create model 
        model                  =   'IFAM';        % name of model
        IFAM                   =   Cplex(model); % define the new model
        IFAM.Model.sense       =   'minimize';
        
        % second model?
        model2                 =   'MILP';
        MILP                   =   Cplex(model2);
        MILP.Model.sense       =   'minimize';
        

        %   Add the columns
        DV                      = Gk ...              % DV for FAM (y_ak)
                                  + K*Lf  ...         % DV for FAM (f_ik)
                                  + numel(col(:,1));  % DV for PMF part

        obj                     =   [y_cost; f_cost; t_cost] ; % y_ak, f_ik, t_pr
        lb                      =   zeros(DV, 1);   % Lower bounds
        ub                      =   inf(DV, 1);     % Upper bounds     
        
        %% Naming DVs
        l = 1;   % Array with DV names  (OPTIONAL, BUT HELPS READING THE .lp FILE)
        for k = 1:k
            for n = 1:GA(k)                         % of the x_{ij} variables
                NameDV (l,:)  = ['Y_' num2str(k,'%02d') '_' num2str(n,'%03d')];
                l = l + 1;
            end
        end
        for k = 1:K 
            for i = 1:Lf                         
                NameDV (l,:)  = ['F_' num2str(k,'%02d') '_' num2str(i,'%03d')];
                l = l + 1;
            end
        end        
        for pr = 1:numel(col(:,1))
            NameDV (l,:)  = ['T_' num2str(pr,'%06d') ];
            l = l + 1;
        end
        
        ctype       = [];
        ctype2      = char(ones(1, (DV)) * ('I')); % for the MILP
        
        % Add the Columns
        
        IFAM.addCols(obj, [], lb, ub, ctype, NameDV);
        MILP.addCols(obj, [], lb, ub, ctype2, NameDV);
        
   %%  Constraints
   %%
    % 1. Capacity constraint (for actual flight (without 'busflights')
        for i = 1:Lf
            C11 = zeros(1,DV); % Spilled pax from p
            C12 = zeros(1,DV); % Recaptured pax from r to p
            C13 = zeros(1,DV); % Capacity for each flight leg
            
            for pr = 1:numel(col(:,1))
                % Activated the spilled passengers
                if delta{col(pr,1),1}(i) ~= 0 
                    C11(Tindex(pr))= 1;
                end
                % Activate the recaptured passengers in p (use Bpr)
                if delta{col(pr,2),1}(i) ~= 0 && ismember(col(pr,2), col(:,1))
                    C12(Tindex(pr)) = Bpr(col(pr,1),col(pr,2)); 
                end
            end
            % Calculate/Select the capacity for this constraint
            for k = 1:K
                C13(Findex(k,i)) = AC.Seats(k);
            end
            % Add all together and add row to cplex
            C1 = C13 + C11 - C12;
            IFAM.addRows(Q(i), C1, inf, sprintf('Capacity_%d',i));
            MILP.addRows(Q(i), C1, inf, sprintf('Capacity_%d',i));
        end
        
        
    % 2. Capacity constraint (for only 'busflights')
        for i = Lf+1:L
            C21 = zeros(1,DV); % Spilled pax from p
            C22 = zeros(1,DV); % Recaptured pax from r to p
            
            for pr = 1:numel(col(:,1))
                if delta{col(pr,1),1}(i) ~= 0 
                    C21(Tindex(pr))= 1;
                end
                if delta{col(pr,2),1}(i) ~= 0 && ismember(col(pr,2), col(:,1))
                    C22(Tindex(pr)) = Bpr(col(pr,1),col(pr,2)); 
                end
            end
            C2 = C21 - C22; % 4 busses with 54 seats for each flight
            IFAM.addRows((Q(i)-216), C2, inf, sprintf('Capacitie_%d',i));
            MILP.addRows((Q(i)-216), C2, inf, sprintf('Capacitie_%d',i));
        end
        
    % 3. Each flight is operated by 1 aircraft type
        for i = 1:Lf
            C3 = zeros(1,DV);
            for k = 1:K
                C3(Findex(k,i)) = 1;
            end
            IFAM.addRows(1, C3, 1, sprintf('Flights_%03d',i));
            MILP.addRows(1, C3, 1, sprintf('Flights_%03d',i));
        end
        
    % 4. What goes in, goes out BY CLAUDIA
         for k = 1:K
            N_k = size(timespace(k).node,1);
            for n = 1:N_k
               C4 = zeros(1,DV);
               [idx_ofl,idx_ifl, idx_on, idx_in] = fik_idx(timespace,k,n);
               for ii = 1:size(idx_ofl,1)
                   C4(Findex(k,idx_ofl(ii,1))) =  1;
               end
               for ii = 1:size(idx_ifl,1)
                   C4(Findex(k,idx_ifl(ii,1))) = -1;
               end            
                C4(Yindex(k, idx_on, timespace)) =  1;
                C4(Yindex(k, idx_in, timespace)) = -1;
                IFAM.addRows(0, C4, 0, sprintf('Inandout_%d_%0003d',k,n));
                MILP.addRows(0, C4, 0, sprintf('Inandout_%d_%0003d',k,n));
            end
        end
    % 5. Fleet size is not exceeded
        for k = 1:K
           C5 = zeros(1,DV);
           G  = numel(timespace(k).ga.Loc);
           NGA = numel(timespace(k).nga.Loc);
           for a = G+1:(G+NGA) %set of overnight ground arcs
               C5(Yindex(k,a,timespace)) = 1; 
           end
           [idx_fls,no] = cl_k(timespace,k);
           for i = 1:no
               C5(Findex(k,idx_fls(i))) = 1;
           end
           IFAM.addRows(0, C5, AC.Units(k), sprintf('Fleetsize_%d',k));
           MILP.addRows(0, C5, AC.Units(k), sprintf('Fleetsize_%d',k));
        end

    %%  Execute model 
    %%
    sol = IFAM.solve();
    IFAM.writeModel([model '.lp']);
    OV = [];
    OV = [OV;IFAM.Solution.objval];
    dual    = IFAM.Solution.dual;
    %RMP.writeModel([model '.lp']);
    constraints = size(IFAM.Model.A,1);
    Alist = [IFAM.Model.A];
 
 time_initial = toc;

%% SUPER BIG LOOP
%%
titer =0;
while not_opt_col == 1 || not_opt_row == 1

    titer = titer +1;
    disp('-------------------------------------------------');
    disp(['Iteration: ',num2str(titer)]);   
    disp('-------------------------------------------------');
    
    dual    = IFAM.Solution.dual;

    
    
%% COLUMN GENERATION
%%
iter=0;
   while not_opt_col == 1
   %% Iterations
   %%
   iter = iter +1;
   no_it(titer,1) = iter;
    disp('-------------------------------------------------');
    disp(['Column iteration: ',num2str(iter)]);   
    disp('-------------------------------------------------');
    
    
   %% Pricing Problem
   %%
   % Get dual variables 
    pi   = dual(1:L);                 % dual variables of capacity constraints (slack)
    sigma = zeros(P,1);
    if size(dual) > constraints
        for ro = 1:numel(rowz)
            sigma(rowz(ro)) = dual(constraints+ro);
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
                    A       = IFAM.Model.A(:,2378+recap_p(re));
                    
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
                    lb      = 0;
                    ub      = inf;
                    NameDV  = ['T_' num2str(recap_p(re)),'_' num2str(recap_r(re))];
                    
                    % ctype and adding Cols
                    ctype       = [];
                    ctype2     = char(ones(1, (1)) * ('I'));
                    IFAM.addCols(obj,A,lb,ub,ctype,NameDV);
                    MILP.addCols(obj,A,lb,ub,ctype2,NameDV);

                end
            end
    end
    
     
    IFAM.solve();
    OV = [OV;IFAM.Solution.objval];
    IFAM.writeModel([model '.lp']);
    
    dual = IFAM.Solution.dual;
    pi   = dual(1:L);                 % dual variables of capacity constraints (slack)
    sigma = zeros(P,1);
    if size(dual) > constraints
        for ro = 1:numel(rowz)
            sigma(rowz(ro)) = dual(constraints+ro);
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
  no_it(titer,2) = iteration;
    disp('-------------------------------------------------');
    disp(['Row iteration: ',num2str(iteration)]);   
    disp('-------------------------------------------------');
  
  DV                      = Gk + K*Lf ...       % DV for FAM part
                            + numel(col(:,1));  % DV for PMF part (incl new)
  primal = IFAM.Solution.x(Gk+Lf*4+1:end);
     
        %% Separation Problem
    % 2. Demand constraint
    check = 0;
    for p = 1:P
        Bc = find(col(:,1)==p);
        
        if sum(primal(Bc)) > demand(p)
            D = find(p==rowz);
            
            if isempty(D) == 1
                check = check + 1;
                C = zeros(1,DV);
                C(Tindex(Bc)) = 1;
                IFAM.addRows(-inf, C, demand(p), sprintf('Demand_%03d',p));
                MILP.addRows(-inf, C, demand(p), sprintf('Demand_%03d',p));
                rowz = [rowz;p];
            end
        end
    end 

    
           %%  Execute model 
    %   Run CPLEX
    IFAM.solve();
    OV = [OV;IFAM.Solution.objval];

    
%     %   Get dual variables
%     primal  = IFAM.Solution.x;   
%     dual    = IFAM.Solution.dual;
    
    
    if check == 0
        not_opt_row = 0;
    else
        not_opt_row = 1;
        not_opt_col = 1;
    end       
        
   end
    
    rowiez = [rowiez; size(rowz)];
    coliez = [coliez; size(colz)];
    disp(rowiez)
    disp(coliez)
end
time_loops = toc;

% SHOULD MILP BE HERE? HOW?
MILP.solve();
MILP.writeModel([model2 '.lp']);

OV = [OV;MILP.Solution.objval];

% NEW ROW, final checking
  primal = MILP.Solution.x(Gk+Lf*4+1:end);
     
        %% Adding rows after optimal MILP
for p = 1:P
    Bc = find(col(:,1)==p);

    if sum(primal(Bc)) > demand(p)
        D = find(p==rowz);

        if isempty(D) == 1
            C = zeros(1,DV);
            C(Tindex(Bc)) = 1;
            %IFAM.addRows(-inf, C, demand(p), sprintf('Demand_%03d',p));
            MILP.addRows(-inf, C, demand(p), sprintf('Demand_%03d',p));
            rowz = [rowz;p];
        end
    end
end 
rowiez = [rowiez; size(rowz)]; 

MILP.solve(); % Final solving
MILP.writeModel([model2 '.lp']);

% Make lists of solution, so it is easy to print
OV = [OV;MILP.Solution.objval];
Final = size(B,1)*4500 + MILP.Solution.objval ;
OV_Final = OV + size(B,1)*4500;
disp(Final)

time_final = toc;
%% Post Processing 
    for kk = 2:(size(GA,2)+1)
        ga_k(kk) = sum(GA(1:kk-1));%#ok<SAGROW>
    end
    
    for k = 1:K
        solution(k).groundarc          = MILP.Solution.x(ga_k(k)+1:ga_k(k+1)); %#ok<SAGROW>
        solution(k).overnight          = solution(k).groundarc(end-size(timespace(k).nga,1)+1:end);
        ac_verif(k)                    = sum(solution(k).overnight);
        [idx_fls,no] = cl_k(timespace,k);
    end

    n  = Gk;
    cp  = n + Lf;
    for k = 1:K
        solution(k).fl                 = timespace(k).fl(logical(MILP.Solution.x(n+1:cp)),:);
        n = n +Lf;
        cp = cp +Lf;
    end
    
time_end = toc;

    %% Export to .txt
    % Flights operated by B737-700
    writetable(solution(3).fl,'B737.txt','Delimiter',' ')  
    % Flights operated by B737-800
    writetable(solution(4).fl,'B738.txt','Delimiter',' ')  
    
%% Post Processing
%%
% Spilage computation
spilled_it = [];%int64( find(MILP.Solution.x(Gk+Lf*4+1:Gk+Lf*4+737)) );
spilled_am = [];
for spill = 1:79
    if int64(MILP.Solution.x(Gk+Lf*4+737+spill)) ~= 0
        spilled_am = [spilled_am; MILP.Solution.x(Gk+Lf*4+737+spill)];
        spilled_it = [spilled_it; colz(spill,:)];
    end
end
spilled_it = spilled_it-1; % as the itineraries in given info start from 0
spillage = sortrows([spilled_it, spilled_am]); % sort on itinerary p
spillage = [spillage(1:11,1),spillage(1:11,2),spillage(1:11,3),... % set it in different columns for the report
    spillage(12:22,1),spillage(12:22,2),spillage(12:22,3),...
    spillage(23:33,1),spillage(23:33,2),spillage(23:33,3)];
spillage = array2table(spillage,'VariableNames',...
   {'ItineraryP','ItineraryR','Pax','Itineray','Itinerayr','Pa','Itinera','Itinerayrr','Paxt'}); % make it a table
writetable(spillage,'spillage.txt','Delimiter',' ')  % write the table

spilled_tot = sum(MILP.Solution.x(Gk+Lf*4+737+1:Gk+Lf*4+737+79)); % total amount of spilled passengers


% Time computation
theend = toc;
times = [time_start;time_pre;time_initial;time_loops;time_final;time_end;theend];


%% Evolution of the Cost
% The Plot
% xlabels = {"Initial", "C1", "C2", "C3", "R1", "R2", "R3", "R4", "R5",...
%         "R6", "R7", "R8", "R9", "R10", "R11", "C1", "C2", "C3",...
%         "R1", "R2", "C1", "C2", "R1", "MILP", "MILP"};
% xlabels1 = ["Initial"; "C1"; "C2"; "C3"; "R1"; "R2"; "R3"; "R4"; "R5";...
%         "R6"; "R7"; "R8"; "R9"; "R10"; "R11"; "C1"; "C2"; "C3";...
%         "R1"; "R2"; "C1"; "C2"; "R1"; "MILP"; "MILP"];
% x_values = 1:25;
% plot(OV_Final)
% xlabel('Iterations')
% ylabel('Final cost [$]')
% line('XData', [1 1], 'YData', [2650000 2740000], 'LineStyle', '-', ...
%     'LineWidth', 1, 'Color','k')
% line('XData', [15 15], 'YData', [2650000 2740000], 'LineStyle', '-', ...
%     'LineWidth', 1, 'Color','k')
% line('XData', [20 20], 'YData', [2650000 2740000], 'LineStyle', '-', ...
%     'LineWidth', 1, 'Color','k')
% line('XData', [23 23], 'YData', [2650000 2740000], 'LineStyle', '-', ...
%     'LineWidth', 1, 'Color','k')
% set(gca,'Xtick',x_values,'XTickLabel',xlabels);

% The Table
%costtable = [xlabels1, OV_Final];


% %% B737 spilled flights
% 
% for pa = spilled_it
%     hoola = find(delta{pa+1,1}(:));
%     for fli = 1:57
%         if isequal(flightnrs(hoola),solution(3).fl(fli,1))
%             sp = [sp; solution(3).fl(fli,1)] % list of spillage flights for B737
%             it_sp = [it_sp; pa] % list of itineraries that B737 flies, with spillage
%         end
%     end
% end




%% Index Functions
%%
function out = Yindex(k, a, timespace)
    GA  = zeros(1,size(timespace,2));
    for kk  = 1:size(timespace,2)
        GA(kk) = size(timespace(kk).gat,1);
    end
    GA_k    = zeros(size(GA));   
    for kk = 2:size(GA,2)
        GA_k(kk) = sum(GA(1:kk-1));
    end
    out = GA_k(k) + a;
end

function out = Findex(k, i)
    LF = 208;
    out = 1546 + (k - 1) * LF + i;  % Function given the variable index for each F(k,i)
end

function out = Tindex(pr)
    out = 2378 + pr;  % Function given the variable index for each T(p,r)
end

















