clear all
clear vars

filename     = 'Assignment2.xlsx';

%% Input 
% import chedule of airline
[~,~,schedule] = xlsread(filename,1,'A2:I233');
[~,varnames]   = xlsread(filename,1,'A1:I1'); 
% if flight cannot be operated by aircraft type assign a verify high cost

schedule(strcmp(schedule, 'NA')) = {1000000};

schedule       = cell2table(schedule,...
                 'VariableNames',varnames);

             
% represent departure and arrival times in a day in minutes 
% from 0 (0:00) to 1439 (23:59)  

[h,m]                 = hms(datetime(schedule.Departure,'ConvertFrom','excel'));
schedule.Departure    = h*60+m;            % departure time in minutes
[h,m]                 = hms(datetime(schedule.Arrival,'ConvertFrom','excel'));
schedule.Arrival      = h*60+m;            % arrival time in minutes
             
% sort flights based on the origin airport and then by the departure time             
schedule = sortrows(schedule,[2,4]);

% determine the locations ( airports) in the time space networks

[a,~,X]  = unique(schedule.ORG); % airport names
A        = size(a,1);            % number of airports 
             
% separate the flights between the hub airports as they are operate by buses

indx             = find(strcmp(schedule.ORG, 'AEP') & strcmp(schedule.DEST, 'EZE')); 
indx             = [indx ;find(strcmp(schedule.ORG, 'EZE') & strcmp(schedule.DEST, 'AEP'))];
B                = schedule(indx,:); % Bus schedule
B                = B(:,1:5);         % Remove costs for each aircraft types as they are operated by buses
B.Cost           = 4500*ones(size(B,1),1); % add cost of operating a bus trip
schedule(indx,:) = []; % remove bus trips from airline schedule


% locatios ( airports) in time-space network
L        = size(schedule.ORG,1); % number of flight arcs
 
% aircraft data per type
[~,varnames] = xlsread(filename,4,'A1:D1');
[~,~,AC]   = xlsread(filename,4,'A2:D5');
AC = cell2table(AC,'VariableNames',varnames);

%% L: Flight arcs for each aircraft type k 
FL.A330         = schedule(:,1:6); 
FL.A330.Arrival = FL.A330.Arrival + AC.TAT(1);

FL.A340         = [schedule(:,1:5) schedule(:,7)]; 
FL.A340.Arrival = FL.A340.Arrival + AC.TAT(2);

FL.B737         = [schedule(:,1:5) schedule(:,8)]; 
FL.B737.Arrival = FL.B737.Arrival + AC.TAT(3);

FL.B738         = [schedule(:,1:5) schedule(:,9)]; 
FL.B738.Arrival = FL.B738.Arrival + AC.TAT(4);

%% N_k: Nodes for each aircraft type k
% Generate a vector of 1:1439 for each airport that represents the timeline
% at which the nodes can be located. 

N_k = zeros(1,A); 

 for j = 1:A
    idx          = [strcmp(FL.A330.ORG, a(j)) strcmp(FL.A330.DEST, a(j))] ;    % index of nodes departing and arriving at airport j
    node_j_k     = sortrows([table2cell(FL.A330(idx(:,1),[2 4])) ; table2cell( FL.A330(idx(:,2),[3 5]))],2);
    G_k          = 
 end 



% flight arc determination for each aircraft

% create table of flight arcs. 

% nodes.depart  = sortrows([ org num2cell(dep_time_min) ],[1,2]); % node 
% nodes.arr     = sortrows([ dest num2cell(arr_time_min)],[1,2]); % node



