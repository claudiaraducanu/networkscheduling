%function Triple_R()
%%  Initialization
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64')
    clc
 	clearvars
    close all

%% Assignment 1, Problem 1
% try        %  Initiate CPLEX model
%  Create model
    model                 =   'Assi_1_PR_2';
    RMP                   =   Cplex(model);
    RMP.Model.sense       =   'minimize';
    
    
%% Define base workspace (cost) Variables and reshape where necessary 
%%

%Flights
[num,fl_nr,all]=xlsread('Input_AE4424_Ass1P2.xlsx',1,'A2:A233');
[num2,fl_ori,all]=xlsread('Input_AE4424_Ass1P2.xlsx',1,'B2:B233');
[num3,fl_desti,all]=xlsread('Input_AE4424_Ass1P2.xlsx',1,'C2:C233');
[num4,text4,all]=xlsread('Input_AE4424_Ass1P2.xlsx',1,'D2:D233');
dep = datestr(num4, 'HH:MM:SS');
STD = datetime(dep);
[num5,text5,all]=xlsread('Input_AE4424_Ass1P2.xlsx',1,'E2:E233');
arr = datestr(num5, 'HH:MM:SS');
STA = datetime(arr);
[fl_cap,text6,all]=xlsread('Input_AE4424_Ass1P2.xlsx',1,'F2:F233');

%passenger itineraries
[pass_itin,text7,all]=xlsread('Input_AE4424_Ass1P2.xlsx',2,'A2:A738');
[num8,pass_ori,all]=xlsread('Input_AE4424_Ass1P2.xlsx',2,'B2:B738');
[num9,pass_desti,all]=xlsread('Input_AE4424_Ass1P2.xlsx',2,'C2:C738');
[pass_demand,text10,all]=xlsread('Input_AE4424_Ass1P2.xlsx',2,'D2:D738');
[pass_fare,text11,all]=xlsread('Input_AE4424_Ass1P2.xlsx',2,'E2:E738');
[~,pass_leg1,all]=xlsread('Input_AE4424_Ass1P2.xlsx',2,'F2:F738');
[~,pass_leg2,all]=xlsread('Input_AE4424_Ass1P2.xlsx',2,'G2:G739');

%recapture rates
[recap_from,text14,all]=xlsread('Input_AE4424_Ass1P2.xlsx',3,'A2:A300');
[recap_to,text15,all]=xlsread('Input_AE4424_Ass1P2.xlsx',3,'B2:B300');
[recap_rate,text16,all]=xlsread('Input_AE4424_Ass1P2.xlsx',3,'C2:C300');
[recap_fare_from,text17,all]=xlsread('Input_AE4424_Ass1P2.xlsx',3,'D2:D300');
[recap_fare_to,text18,all]=xlsread('Input_AE4424_Ass1P2.xlsx',3,'E2:E300');
% 
% %Flights
% [num,fl_nr,all]=xlsread('Input_Example',1,'A2:A7');
% [num2,fl_ori,all]=xlsread('Input_Example',1,'B2:B7');
% [num3,fl_desti,all]=xlsread('Input_Example',1,'C2:C7');
% [num4,text4,all]=xlsread('Input_Example',1,'D2:D7');
% dep = datestr(num4, 'HH:MM:SS');
% STD = datetime(dep);
% [num5,text5,all]=xlsread('Input_Example',1,'E2:E7');
% arr = datestr(num5, 'HH:MM:SS');
% STA = datetime(arr);
% [fl_cap,text6,all]=xlsread('Input_Example',1,'F2:F7');
% 
% %passenger itineraries
% [pass_itin,text7,all]=xlsread('Input_Example',2,'A2:A9');
% [num8,pass_ori,all]=xlsread('Input_Example',2,'B2:B9');
% [num9,pass_desti,all]=xlsread('Input_Example',2,'C2:C9');
% [pass_demand,text10,all]=xlsread('Input_Example',2,'D2:D9');
% [pass_fare,text11,all]=xlsread('Input_Example',2,'E2:E9');
% [~,pass_leg1,all]=xlsread('Input_Example',2,'F2:F9');
% [~,pass_leg2,all]=xlsread('Input_Example',2,'G2:G10');
% 
% %recapture rates
% [recap_from,text14,all]=xlsread('Input_Example',2,'A14');
% [recap_to,text15,all]=xlsread('Input_Example',2,'B14');
% [recap_rate,text16,all]=xlsread('Input_Example',2,'C14');
% [recap_fare_from,text17,all]=xlsread('Input_Example',2,'D14');
% [recap_fare_to,text18,all]=xlsread('Input_Example',2,'E14');




%% Define decision varaibles
%%

sol_list = [];
temp = 0;

demand_list = [];

for i = 1:size(fl_nr,1)
    temp = 0;
    for j = 1:size(pass_itin,1)
        if strcmp(fl_nr(i),pass_leg1(j)) == 1
            temp = temp + pass_demand(j);
        end
        if strcmp(fl_nr(i),pass_leg2(j)) == 1
            temp = temp + pass_demand(j);
        end
    end
    demand_list = [demand_list;temp];
end

lst = minus(demand_list,fl_cap);

    

deficit = (demand_list > fl_cap);
precies_goed = (demand_list == fl_cap);
excess = (demand_list < fl_cap);

opt_row = 0;
opt_col = 0;

zeven_drie_zeven = size(pass_itin,1);
twee_drie_twee = size(fl_nr,1);

van_naar= [[0:zeven_drie_zeven-1]',ones(zeven_drie_zeven,1)*999];



%Create actual DV's
        DV      =     size(pass_itin,1);                   %[i,j],k1; [i,j],k2;... 
        
    
%%  Objective Function

cost = pass_fare;

    obj                     =   cost';                                     % DV coefficients in the OF
    lb                      =   zeros(DV, 1);                              % Lower bounds
    ub                      =   inf(DV, 1);                                % Upper bounds
%     ctype                   =   char(ones(1, (DV)) * ('C'));             % Variable types 'C'=continuous; 'I'=integer; 'B'=binary


%% Setting up the Master template of the problem

% cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
% RMP.addCols(obj', [], lb, ub,[] ,NameDV);
RMP.addCols(obj', [], lb, ub)
% RMP.addCols(obj.o, [], obj.lb, obj.ub);


%% Constraint modeling
%%

%Spill itineraries
for fl = 1:size(fl_nr)
    C1      =   zeros(1, DV);    %Setting coefficient matrix with zeros
    for i = 1:size(pass_itin)
        if strcmp(fl_nr(fl),pass_leg1(i))
            C1(i) = 1;
        end
        if strcmp(fl_nr(fl),pass_leg2(i))
            C1(i) = 1;
        end
    end
    RMP.addRows(lst(fl), C1, inf,sprintf('Spill_FL_%03d',fl));
end



%% Solving the model
%%  Execute model
%     cplex.Param.mip.limits.nodes.Cur            = 1e+8;         %max number of nodes to be visited (kind of max iterations)
%     cplex.Param.timelimit.Cur                   =  180;         %max time in seconds


%   Run CPLEX
    opl = RMP.solve();
    sol_list = [sol_list;RMP.Solution.objval];
%     RMP.writeModel([model '.lp']);
%     RMP.writeSolution([model '.sol']);
    


%% Clear Vars for second run
% clearvars -except 


%  cplex.Solution.objval    
%     Balanced_sol(:,:,ahhh) = cplex.Solution.x
% [NameDV,char(ones((DV),1) * ('=')),num2str(Fuel_sol,'%02d')]
% k = waitforbuttonpress
% close all
% [NameDV,char(ones((DV),1) * ('=')),num2str(cplex.Solution.x,'%02d')]
% found = find(RMP.Solution.x);

% [NameDV(found,:),char(ones((size(found,1)),1) * ('=')),num2str(RMP.Solution.x(found),'%02d')]

% RMP.Solution.objval

% dual = RMP.Solution.dual;


[find(RMP.Solution.x),RMP.Solution.x(find(RMP.Solution.x))];
dual = RMP.Solution.dual;

added_cols = [];
added_rows = [];

in_col = 0;
in_row = 0;
in_grand = 0;




%% The grand loop
%%


while opt_col + opt_row ~= 2

    %% Column generation
    %%
    dual = RMP.Solution.dual;

    prev = size(RMP.Model.A,2)-zeven_drie_zeven;

    while opt_col ~= 1

        test_list = [];
        t_p = [];
        % determening the reduced prices of the different paths
        for i = 1:size(recap_rate,1)


            %FROM PRICE PROCESSING -------------

            price_from = recap_fare_from(i);
            from = recap_from(i); % The P's that can be recaptured

            first_from = pass_leg1(from+1); % First leg that can be recaptured
            second_from = pass_leg2(from+1); % Second leg that can be recaptured

            for j = 1:size(fl_nr,1) % for all L
                if strcmp(first_from,fl_nr(j)) == 1 % if
                    first_nr_from = j;
                end
                if strcmp(second_from,fl_nr(j)) == 1
                    second_nr_from = j;
                elseif sum(strcmp(second_from,fl_nr)) == 0
                    second_nr_from = 0;
                end
            end


            if second_nr_from == 0
                sum_pi_from = (dual(first_nr_from));
            else
                sum_pi_from = (dual(first_nr_from)+dual(second_nr_from));
            end



            %TO PRICE PROCESSING -------------
            price_to = recap_fare_to(i);
            to = recap_to(i);

            first_to = pass_leg1(to+1);
            second_to = pass_leg2(to+1);


            for j = 1:size(fl_nr,1)
                if strcmp(first_to,fl_nr(j)) == 1
                    first_nr_to = j;
                end
                if strcmp(second_to,fl_nr(j)) == 1
                    second_nr_to = j;
                elseif sum(strcmp(second_to,fl_nr)) == 0
                    second_nr_to = 0;
                end
            end



            if second_nr_to == 0
                sum_pi_to = (dual(first_nr_to));
            else
                sum_pi_to = (dual(first_nr_to)+dual(second_nr_to));
            end

            %REDUCED COST --------------------------

            rr = recap_rate(i);
            
            if size(added_rows,1) == 0
                reduced_cost = (price_from -  sum_pi_from) - rr*(price_to - sum_pi_to);
            else
                if sum(from == added_cols) > 0
                    entry = find(from == added_cols) + twee_drie_twee;
                    reduced_cost = (price_from -  sum_pi_from) - rr*(price_to - sum_pi_to) - dual(entry);
                else
                    reduced_cost = (price_from -  sum_pi_from) - rr*(price_to - sum_pi_to);
                end
            end
            t_p = [t_p;reduced_cost];



            % Adding the columns if necessary -------------------------------------
            if reduced_cost < 0 && sum(i ~= added_cols) == size(added_cols,2) %add column only if not yet in there and reduced cost < 0
                added_cols = [added_cols,i];

                obj = price_from - rr * price_to; %fare of new path
                A = zeros(size(RMP.Model.A,1),1);

                %add indices for the paths considered in the RMP
                A(first_nr_from) = 1;
                if second_nr_from ~= 0
                    A(second_nr_from) = 1;
                end

                A(first_nr_to) = -rr;
                if second_nr_to ~= 0
                    A(second_nr_to) = -rr;
                end
                
                %add the correct index in the A matrix for the added
                %seperation problem rows
                
                
                if sum(from == added_rows) > 0
                    indices = find(from == added_rows);
                    A(232+indices) = 1;
                end
                
                size(A)
                
                
                lb = [0];
                ub = [inf];
                ctype = [];
%                 name = ['From_' num2str(from,'%03d') '_To_' num2str(to,'%03d')];


                RMP.addCols(obj,A,lb,ub);
                van_naar = [van_naar; from,to];


                test_list = [test_list,A];
                opt_row = 0;
            end
        end

        RMP.solve();
        sol_list = [sol_list;RMP.Solution.objval];
        dual = RMP.Solution.dual;

        found = find(t_p<0);

        together = [added_cols,found'];
        together = unique(together);

        current = size(added_cols,2);
        increase = current - prev;

        if size(together,2) == size(added_cols,2) && increase == 0
            opt_col = 1;
            in_col = in_col + 1;
        end

        prev = current;
    end


    %% Row generation
    %%
    
    % Adding the relevant rows to the problem

    vorige = size(RMP.Model.A,1)-twee_drie_twee;
    
    while opt_row ~= 1
        oplossing = RMP.Solution.x;
        cols = size(RMP.Model.A,2);

        for p = 0:(size(pass_itin)-1)

            onze_p = (van_naar(:,1) == p);

            de_tp = times(oplossing,onze_p);
            links = sum(de_tp);
            rechts = pass_demand(p+1);

            if links > rechts && sum(p ~= added_rows) == size(added_rows,2)
                C1 = zeros(1,cols);
                C1(find(onze_p)) = 1;

                RMP.addRows(-inf, C1, rechts,sprintf('Seperation_problem_path%03d',p));
                added_rows = [added_rows,p];

                opt_col = 0;
            end

        end


        RMP.solve()
        sol_list = [sol_list;RMP.Solution.objval];


%         samen = [added_rows,found'];
%         samen = unique(samen);

        huidige = size(added_rows,2);
        toename = huidige - vorige;

        if toename == 0
            opt_row = 1;
            in_row = in_row + 1;
        end

        vorige = huidige;

    end
 

end



%% Post Processing in order to verify data
%%
close all





pax = [];
pax_l = [];
pax_j = [];
pax_tot = 0;

uitkomst = RMP.Solution.x;
for i = 1:zeven_drie_zeven
    itin = i - 1;
    
    leaving = 0;
    joining = 0;
    
    
    for j = 1:size(van_naar,1)
        if  van_naar(j,1) == itin
            leaving = leaving + uitkomst(j);
        end
        if  van_naar(j,2) == itin
            
            for k = 1:size(recap_rate,1)
                if van_naar(j,1) == recap_from(k) && van_naar(j,2) == recap_to(k)
                    recap = recap_rate(k);
                end
            end            
            joining = joining + uitkomst(j)*recap;
        end
        
    end
    
    pax_tot = pax_tot + leaving - joining;
    total_demand = pass_demand(i) - leaving + joining;
    pax = [pax;total_demand];
    pax_l = [pax_l;leaving];
    pax_j = [pax_j;joining];
end



temp = 0;

new_demand = [];

for i = 1:size(fl_nr,1)
    temp = 0;
    for j = 1:size(pass_itin,1)
        if strcmp(fl_nr(i),pass_leg1(j)) == 1
            temp = temp + pax(j);
        end
        if strcmp(fl_nr(i),pass_leg2(j)) == 1
            temp = temp + pax(j);
        end
    end
    new_demand = [new_demand;temp];
end
[new_demand,fl_cap,minus(fl_cap,new_demand)]
RMP.Solution.objval



%% path to the solution 
%%

plot([1:size(sol_list,1)],sol_list)



%% What are the best fligths to increase
%%

minmin = minus(fl_cap,new_demand);
opop = [minmin,dual([1:232])];
vluchten = [1 2 8 25 41 80 103 144 163 214 219:232];
xc = [vluchten',opop(vluchten,:)]
BBB = sortrows(xc,3)