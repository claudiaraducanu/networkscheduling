%% Assignment 2


tic
%function Robin&Bert()
%% Initialization
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');   
clc
clearvars
close all

%Create CPLEX Model
model                 =   'Assignment_2';
RMP                   =  Cplex(model);
RMP.Model.sense       =  'minimize';
    
    
%% Pre-processing step 
%%
types  = 5;

[num,fl_nr,all]=xlsread('Assignment2','Flight','A2:A233');
[num2,org,all]=xlsread('Assignment2','Flight','B2:B233');
[num3,dest,all]=xlsread('Assignment2','Flight','C2:C233');
[num4,text4,all]=xlsread('Assignment2','Flight','D2:D233');
dep = datestr(num4, 'HH:MM:SS');
STD = datetime(dep);
[num5,text5,all]=xlsread('Assignment2','Flight','E2:E233');
arr = datestr(num5, 'HH:MM:SS');
STA = datetime(arr);
[OC330,text33,all]=xlsread('Assignment2','Flight','F2:F233');
[OC340,text34,all]=xlsread('Assignment2','Flight','G2:G233');
[OC737,text37,all]=xlsread('Assignment2','Flight','H2:H233');
[OC738,text38,all]=xlsread('Assignment2','Flight','I2:I233');
[OCbus,text39,all]=xlsread('Assignment2','Flight','J2:J233');

% adding in nan for cells that are skipped in xlsread and creating grand
% list

OC330 = [nan(2,1);OC330;nan(14,1)];
OC340 = [nan(2,1);OC340;nan(14,1)];
OC737 = [nan(4,1);OC737;nan(16,1)];
OC738 = [nan(4,1);OC738;nan(16,1)];

All_the_lists = [OC330,OC340,OC737,OC738,OCbus];


% back to the file reading

[pass_itin,text7,all]=xlsread('Assignment2','Itinerary','A2:A738');
[num8,pass_ori,all]=xlsread('Assignment2','Itinerary','B2:B738');
[num9,pass_desti,all]=xlsread('Assignment2','Itinerary','C2:C738');
[pass_demand,text10,all]=xlsread('Assignment2','Itinerary','D2:D738');
[pass_fare,text11,all]=xlsread('Assignment2','Itinerary','E2:E738');
[num12,pass_leg1,all]=xlsread('Assignment2','Itinerary','F2:F738');
[num13,pass_leg2,all]=xlsread('Assignment2','Itinerary','G2:G738');

% leg2 = [nan(1,1);leg2;nan(6,1)];

[recap_from,text14,all]=xlsread('Assignment2','Recapture Rate','A2:A300');
[recap_to,text15,all]=xlsread('Assignment2','Recapture Rate','B2:B300');
[recap_rate,text16,all]=xlsread('Assignment2','Recapture Rate','C2:C300');
[recap_fare_from,text17,all]=xlsread('Assignment2','Recapture Rate','D2:D300');
[recap_fare_to,text18,all]=xlsread('Assignment2','Recapture Rate','E2:E300');

TAT = [50,60,30,35];
cap = [10,4,8,31];
Countertje = 0;
seats = [270,290,128,170,216];


for k = 1:4
    
    Airports = unique(org);
    ArrivalCube = zeros(232,2,34);
    DepartureCube = zeros(232,2,34);

    for i = 1:length(fl_nr)

        FlightNumber = i;
        Origin = org(i);
        Destination = dest(i);
        ArrivalTime = num5(i)+(TAT(k)/1440); 
        DepartureTime = num4(i);

        if ArrivalTime >= 1
            ArrivalTime = ArrivalTime-1;
        end
        
        if sum(num2str(All_the_lists(i,k))) ~= 253;
            Countertje = Countertje +1;
            for j = 1:length(Airports) 
            if strcmp(Airports(j),Origin) == 1
                DepartureCube(i,1,j) = FlightNumber;
                DepartureCube(i,2,j) = DepartureTime;
            end

            if strcmp(Airports(j),Destination) == 1
                ArrivalCube(i,1,j) = FlightNumber;
                ArrivalCube(i,2,j) = ArrivalTime;
            end  
            end
        end
        
    end

    %figure()
    datetick('x','HH:MM')
    hold on
    for j = 1:size(Airports)
        x_start = 0;
        x_end = 1;
        y = [j,j];
        %plot([x_start,x_end],y,'Color',[0,0.7,0.9],'LineWidth',0.2)
    end

    for j = 1:length(Airports)
        indices = find(ArrivalCube(:,1,j));
        x_jes = ArrivalCube(indices,2,j);
        y_jes = ones(size(x_jes,1),1)*j;
        scatter(x_jes,y_jes, '*')

        indices = find(DepartureCube(:,1,j));
        x_jes = DepartureCube(indices,2,j);
        y_jes = ones(size(x_jes,1),1)*j;
        scatter(x_jes,y_jes, 'o')
    end
    
    for i = 1:length(fl_nr)
        if sum(num2str(All_the_lists(i,k))) ~= 253;
            X1 = num4(i);
            X2 = num5(i)+(TAT(k)/1440);

            if X2 >= 1
                X2 = X2-1;
            end

            Origin = org(i);
            Destination = dest(i);

            for j = 1:length(Airports)
                if strcmp(Airports(j),Origin) == 1
                    Y1 = j;
                end     
                if strcmp(Airports(j),Destination) == 1
                   Y2 = j;
                end   
            end
            %plot([X1,X2],[Y1,Y2])
        end
    end
    
    GroundArcs = [];
    Teller = 1;
    Grootte = [];
    SamenModifiedGrootte = [];

    for j = 1:length(Airports)
        ArrivalTest = unique(sortrows(ArrivalCube(:,2,j)));
        DepartureTest = unique(sortrows(DepartureCube(:,2,j)));
        Samen = sortrows([ArrivalTest;DepartureTest]);
        SamenModified = Samen([3:length(Samen)]);

        Grootte = [Grootte;length(SamenModified)];
        SamenModifiedGrootte = [SamenModifiedGrootte;size(SamenModified,1)];
        
        if size(SamenModified) > 0;
            
            if size(SamenModified,1) == 1;
                GroundArcs = [GroundArcs;Teller, j,SamenModified(1),SamenModified(1)];
                Teller = Teller +1;
            end
            if size(SamenModified,1) ~= 1;
                for z = 1:(length(SamenModified)-1)
                    GroundArcs = [GroundArcs;Teller, j,SamenModified(z),SamenModified(z+1)];
                    Teller = Teller +1;
                end  
                GroundArcs = [GroundArcs;Teller, j,SamenModified(z+1),SamenModified(1)];
                Teller = Teller +1;
            end
        end
    end
    
    temp = GroundArcs(1,:);
    for q = 2:size(GroundArcs,1)
        arc_l = round(GroundArcs(q,3),7);
        arc_r = round(GroundArcs(q,4),7);
        trigger = false;
        trigger = or(arc_l ~=  arc_r,sum(GroundArcs(q,2)== GroundArcs(:,2)) == 1);
        if trigger
            temp = [temp;GroundArcs(q,:)];
        end
    end
    
    temp(:,1) = [1:size(temp,1)]';       
    GroundArcs = temp;
        
    
    
    if k == 1;
        GroundArcsA330 = GroundArcs;
    end
    if k == 2;
        GroundArcsA340 = GroundArcs;
    end
    if k == 3;
        GroundArcsB737 = GroundArcs;
    end
    if k == 4;
        GroundArcsB738 = GroundArcs;
    end
end


%% Flight leg assignment decission variables
%%

gr_arcs = size(GroundArcsA330,1)+size(GroundArcsA340,1)+size(GroundArcsB737,1)+size(GroundArcsB738,1);
cost = [];
big_m = 1e7;

for k = 1:types
    for i = 1:size(fl_nr,1)
        
        if All_the_lists(i,k) > 0
            temp_cost = All_the_lists(i,k);
        else
            temp_cost = big_m;
        end
        cost = [cost;temp_cost];
        
    end
end

cost = [cost;zeros(gr_arcs,1)];

%% Defining the objective function
%%

DV = size(cost,1);

    obj                     =   cost;                                     % DV coefficients in the OF
    lb                      =   zeros(DV, 1);                              % Lower bounds
    ub                      =   inf(DV, 1);                                % Upper bounds
%     ctype                   =   char(ones(1, (DV)) * ('C'));             % Variable types 'C'=continuous; 'I'=integer; 'B'=binary



%% Setting up the Master template of the problem
%%

% cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
% RMP.addCols(obj', [], lb, ub,[] ,NameDV);
RMP.addCols(obj, [], lb, ub)


%% Constraint modeling
%%

%Operate all flights
for fl = 1:size(fl_nr,1)
    C1      =   zeros(1, DV);    %Setting coefficient matrix with zeros
    
    for k = 1:types
        index = fl + (k-1)*232;
        C1(index) = 1;
    end
    
    RMP.addRows(1, C1, 1,sprintf('Operate_fl%03d',fl));
end


%Deny operation of flight if NA
for k = 1:types
    for i = 1:size(fl_nr,1)
        C2      =   zeros(1, DV);    %Setting coefficient matrix with zeros
        
        if All_the_lists(i,k) > 0
            nothing = 0;
        else
            index = i + (k-1)*232;
            C2(index) = 1;
            RMP.addRows(0, C2, 0,sprintf('Cannot_operate%03d',fl,'AC_type%02d',k));
        end
        
        
    end
end


%Conservation of aircraft
for k = 1:4
    if k == 1
        grond_arcs = GroundArcsA330;
        counter = size(fl_nr,1)*types;
    end
    if k == 2
        grond_arcs = GroundArcsA340;
        counter = size(fl_nr,1)*types + size(GroundArcsA330,1);
    end
    if k == 3
        grond_arcs = GroundArcsB737;
        counter = size(fl_nr,1)*types + size(GroundArcsA330,1) + size(GroundArcsA340,1);
    end
    if k == 4
        grond_arcs = GroundArcsB738;
        counter = size(fl_nr,1)*types + size(GroundArcsA330,1) + size(GroundArcsA340,1) + size(GroundArcsB737,1);
    end
    


    dep_time = num4;
    arr_time = num5+(TAT(k)/1440);
    arr_time(find(arr_time>1)) = arr_time(find(arr_time>1)) - 1;
    append = [];

    for airp = 1:size(Airports,1)
        for arcs = 1:sum(grond_arcs(:,2) == airp,1)
            C3      =   zeros(1, DV);    %Setting coefficient matrix with zeros

            ap = Airports(airp);



            %in and outflow of arcs

            found = find(grond_arcs(:,2) == airp);
            if arcs == 1
                in_index = found(end);
                out_index = found(arcs);
            else
                in_index = found(arcs-1);
                out_index = found(arcs);
            end

            time = round(grond_arcs(in_index,4),7);

            C3(counter + in_index) = 1;
            C3(counter + out_index) = -1;

            %in and outflow of flights
            found_dep = find(round(dep_time,7) == time);
            found_arr = find(round(arr_time,7) == time);

            departing = [];
            arriving = [];

            for i = 1:size(found_dep)
                if strcmp(org(found_dep(i)),ap) == 1
                    departing = [departing,found_dep(i)];
                end
            end

            for i = 1:size(found_arr)
                if strcmp(dest(found_arr(i)),ap) == 1
                    arriving = [arriving,found_arr(i)];
                end
            end

            if size(departing) > 0
               C3((k-1)*232 + departing) = -1;
            end
            if size(arriving) > 0
                C3((k-1)*232 + arriving) = 1;
            end


            %append = [append;airp,arcs,C3];
            RMP.addRows(0, C3, 0,sprintf('AC_type%02d',k,'Airport%02d',airp,'Arc%02d',arcs));
        end
    end
end


%sum of flights, chosen at midnight
for k = 1:4
    C4      =   zeros(1, DV);    %Setting coefficient matrix with zeros
    
    if k == 1
        grond_arcs = GroundArcsA330;
        base = size(fl_nr,1)*types;
    end
    if k == 2
        grond_arcs = GroundArcsA340;
        base = size(fl_nr,1)*types + size(GroundArcsA330,1);
    end
    if k == 3
        grond_arcs = GroundArcsB737;
        base = size(fl_nr,1)*types + size(GroundArcsA330,1) + size(GroundArcsA340,1);
    end
    if k == 4
        grond_arcs = GroundArcsB738;
        base = size(fl_nr,1)*types + size(GroundArcsA330,1) + size(GroundArcsA340,1) + size(GroundArcsB737,1);
    end

    fl_past_mid = [];
    ga_past_mid = [];

    for i = 1:size(fl_nr,1)
        ArrivalTime = num5(i)+(TAT(k)/1440); 
        DepartureTime = num4(i);
        if (ArrivalTime - DepartureTime) < 0
            fl_past_mid = [fl_past_mid,i];
        end
        if ArrivalTime > 1
            fl_past_mid = [fl_past_mid,i];
        end
    end
    
    fl_past_mid = unique(fl_past_mid);

    for i = 1:size(Airports,1)
        found = max(find(grond_arcs(:,2) == i));
        ga_past_mid = [ga_past_mid,found];
    end
    
    C4( (k-1)*232 + fl_past_mid) = 1;
    C4( base + ga_past_mid) = 1;
    
    
    RMP.addRows(0, C4, cap(k),sprintf('AC_type_%02d',k));
end

base_FAM_cols = size(RMP.Model.A,2);
base_FAM_rows = size(RMP.Model.A,1);

%% PART TWO, ADDING THE PMF --------------------------------------------------------------------------------------------------------
%%

%% objective function
%%

%addition of the 737 ficticious itineraries
maat = size(pass_fare,1);
    obj                     =   pass_fare;                                 % DV coefficients in the OF
    lb                      =   zeros(maat, 1);                              % Lower bounds
    ub                      =   inf(maat, 1);                                % Upper bounds
%     ctype                   =   char(ones(1, (DV)) * ('C'));             % Variable types 'C'=continuous; 'I'=integer; 'B'=binary

RMP.addCols(obj, [], lb, ub)




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


zeven_drie_zeven = size(pass_itin,1);
twee_drie_twee = size(fl_nr,1);

van_naar= [[0:zeven_drie_zeven-1]',ones(zeven_drie_zeven,1)*999];




%% Constraint modeling
%%



%Initial setup of capacity check in each leg
%   not considering recaptures
cols = size(RMP.Model.A,2);

for fl = 1:size(fl_nr,1)
    C5      =   zeros(1, cols);    %Setting coefficient matrix with zeros
    
    %flight capacity of AC operating
    for k = 1:types
        index = fl + (k-1)*232;
        C5(index) = seats(k);
    end
    
    %passengers spilled to ficticious itin   
    for i = 1:size(pass_itin)
        if strcmp(fl_nr(fl),pass_leg1(i))
            index = base_FAM_cols + i;
            C5(index) = 1;
        end
        if strcmp(fl_nr(fl),pass_leg2(i))
            index = base_FAM_cols + i;
            C5(index) = 1;
        end
    end
        
    RMP.addRows(demand_list(fl), C5, inf,sprintf('Capacity_fl_%03d',fl));
end


%% Initial solving of the RMP
%%

opl = RMP.solve();
dual = RMP.Solution.dual;


base_cols_post_RMP = size(RMP.Model.A,2);
base_rows_post_RMP = size(RMP.Model.A,1);


opt_row = 0;
opt_col = 0;

%lists for storing additional information
added_cols = [];
added_rows = [];
%sol_list = [row loops, col loops, total loops, solution]
sol_list = [0,0,0,RMP.Solution.objval];
    
in_col = 0;
in_row = 0;
in_grand = 0;

%% Row and Column generation
%%

while opt_col + opt_row ~= 2

    %% Column generation
    %%
    dual = RMP.Solution.dual;

    prev = size(RMP.Model.A,2)-base_cols_post_RMP;

    while opt_col ~= 1

        test_list = [];
        t_p = [];
        % determening the reduced prices of the different paths
        for i = 1:size(recap_rate,1)
    

            %FROM PRICE PROCESSING -------------

            price_from = recap_fare_from(i);
            from = recap_from(i);

            first_from = pass_leg1(from+1);
            second_from = pass_leg2(from+1);

            for j = 1:size(fl_nr,1)
                if strcmp(first_from,fl_nr(j)) == 1
                    first_nr_from = j;
                end
                if strcmp(second_from,fl_nr(j)) == 1
                    second_nr_from = j;
                elseif sum(strcmp(second_from,fl_nr)) == 0
                    second_nr_from = 0;
                end
            end


            if second_nr_from == 0
                sum_pi_from = (dual(base_FAM_rows+first_nr_from));
            else
                sum_pi_from = (dual(base_FAM_rows+first_nr_from)+dual(base_FAM_rows+second_nr_from));
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
                sum_pi_to = (dual(base_FAM_rows+first_nr_to));
            else
                sum_pi_to = (dual(base_FAM_rows+first_nr_to)+dual(base_FAM_rows+second_nr_to));
            end

            %REDUCED COST --------------------------

            rr = recap_rate(i);
            
            if size(added_rows,1) == 0
                reduced_cost = (price_from -  sum_pi_from) - rr*(price_to - sum_pi_to);
            else
                if sum(from == added_cols) > 0
                    entry = find(from == added_cols) + base_FAM_rows + twee_drie_twee;
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
                A(base_FAM_rows + first_nr_from) = 1;
                if second_nr_from ~= 0
                    A(base_FAM_rows + second_nr_from) = 1;
                end

                A(base_FAM_rows + first_nr_to) = -rr;
                if second_nr_to ~= 0
                    A(base_FAM_rows + second_nr_to) = -rr;
                end
                
                %add the correct index in the A matrix for the added
                %seperation problem rows
                
                
                if sum(from == added_rows) > 0
                    indices = find(from == added_rows);
                    A(base_FAM_rows + twee_drie_twee + indices) = 1;
                end
                

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
        sol_list = [sol_list;in_col,in_row,in_grand,RMP.Solution.objval];
        prev = current;
    end

    %% Row generation
    %%
    
    vorige = size(RMP.Model.A,1)-base_rows_post_RMP;
    
    while opt_row ~= 1
        cols = size(RMP.Model.A,2);
        oplossing = RMP.Solution.x([2647:cols],:);
        
        for p = 0:(size(pass_itin)-1)

            onze_p = (van_naar(:,1) == p);

            de_tp = times(oplossing,onze_p);
            links = sum(de_tp);
            rechts = pass_demand(p+1);

            if links > rechts && sum(p ~= added_rows) == size(added_rows,2)
                C1 = zeros(1,cols);
                C1(base_cols_post_RMP - 737 + find(onze_p)) = 1;

                RMP.addRows(-inf, C1, rechts,sprintf('Seperation_problem_path%03d',p));
                added_rows = [added_rows,p];

                opt_col = 0;
            end

        end


        RMP.solve()
        


%         samen = [added_rows,found'];
%         samen = unique(samen);

        huidige = size(added_rows,2);
        toename = huidige - vorige;

        if toename == 0
            opt_row = 1;
            in_row = in_row + 1;
        end

        vorige = huidige;
        sol_list = [sol_list;in_col,in_row,in_grand,RMP.Solution.objval];

    end

    in_grand = in_grand + 1;
    sol_list = [sol_list;in_col,in_row,in_grand,RMP.Solution.objval];
end

tijd = toc;


%% post processing
%%

allocation_tussen = [RMP.Solution.x(1:232),RMP.Solution.x(233:464),RMP.Solution.x(465:696),RMP.Solution.x(697:928),RMP.Solution.x(929:1160)];
tijd_tussen = toc;
% plot(1:size(sol_list,1),sol_list(:,4))



%% Now ensuring f_i = {0,1}
%%

%creating a new model named MILP with 
%all the parameters of RMP in addition to the specification of the ctype




%Create a second CPLEX Model with different ctype
model                 =   'Assignment_2_MILP';
MILP                  =  Cplex(model);
MILP.Model.sense      =  'minimize';
    

    obj                     =   RMP.Model.obj;                             % DV coefficients in the OF
    lb                      =   RMP.Model.lb;                              % Lower bounds
    ub                      =   RMP.Model.ub;                              % Upper bounds
    the_size = size(RMP.Model.A,2);
    ctype                   =   [char(ones(1, (1160)) * ('B')),char(ones(1, (the_size - 1160)) * ('C'))];             % Variable types 'C'=continuous; 'I'=integer; 'B'=binary

    A                       =   full(RMP.Model.A);
    lhs                     =   RMP.Model.lhs;
    rhs                     =   RMP.Model.rhs;

%setting up the column side of the problem
MILP.addCols(obj,[],lb,ub,ctype)
%adding in the rows
MILP.addRows(lhs, A,rhs);

%% solving the MILP problem
MILP.solve()



%% final row check
%%


cols = size(MILP.Model.A,2);
oplossing = MILP.Solution.x([2647:cols],:);
new_rows = [];
for p = 0:(size(pass_itin)-1)

    onze_p = (van_naar(:,1) == p);

    de_tp = times(oplossing,onze_p);
    links = sum(de_tp);
    rechts = pass_demand(p+1);

    if links > rechts && sum(p ~= added_rows) == size(added_rows,2)
        C1 = zeros(1,cols);
        C1(base_cols_post_RMP - 737 + find(onze_p)) = 1;

        MILP.addRows(-inf, C1, rechts,sprintf('Seperation_problem_path_milp%03d',p));
        new_rows = [new_rows,p];
    end
end

MILP.solve()

    
%% Post processing
%%

allocation_end = [MILP.Solution.x(1:232),MILP.Solution.x(233:464),MILP.Solution.x(465:696),MILP.Solution.x(697:928),MILP.Solution.x(929:1160)];
tijd_end = toc;
