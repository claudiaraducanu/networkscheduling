%function [AC,B,timespace] = read_schedule(filename)
clearvars
clear all
%read_schedule - Read the flight schedule information provided and
%converts it into a set of nodes and arcs that would represent the
%time-space network for each aircraft type the airline owns. 
%
% Inputs:
%    filename - .xls file containing the schedule of the airline 
%   
% Outputs:
%    AC         - [4,4] Table of aircraft specific information.
%    B          - Flights operated by buses
%    timespace  - Structure array with 4 fields`:
%                   - fl: flight arcs in the timespace network for aircraft
%                   type k
%                   - node: nodes in timepsace netork for aircraft type k
%                   - ga: ground arcs in the timespace network for aircraft
%                   type k
%                   - nga: ground arcs in the timespace network for 
%                   aircraft type k
% Author: Group 4
% April 2018;
% Assignment 2, Network Scheduling

%% ------------- BEGIN CODE -----------------------------------------------
filename = 'Assignment2.xlsx';
%% Input 
     %import schedule of airline
    schedule = readtable(filename);   
    
    %% Use this if we decide to include flights that cannot be operated by 
    % as specific aircraft type: if flight cannot be operated by aircraft 
    % type assign a very high cost
    %schedule(strcmp(schedule, 'NA')) = {1000000};
    
%%
    % represent departure and arrival times in a day in minutes 
    % from 0 (0:00) to 1439 (23:59)  

    [h,m]                 = hms(datetime(schedule.Departure,'ConvertFrom','excel'));
    schedule.Departure    = h*60+m;            % departure time in minutes
    [h,m]                 = hms(datetime(schedule.Arrival,'ConvertFrom','excel'));
    schedule.Arrival      = h*60+m;            % arrival time in minutes

    % sort flights based on the origin airport and then by the departure time             
    schedule = sortrows(schedule,[2,4]);

%% Separate the flights between the hub airports as they are operate by buses

    indx             = find(strcmp(schedule.ORG, 'AEP') & ...
                            strcmp(schedule.DEST, 'EZE')); 
    indx             = [indx ;find(strcmp(schedule.ORG, 'EZE') & ...
                            strcmp(schedule.DEST, 'AEP'))];
    B                = schedule(indx,:); % Bus schedule
    B                = B(:,1:5);         % Remove costs for each aircraft types as they are operated by buses
    B.Cost           = 4500*ones(size(B,1),1); % add cost of operating a bus trip
    schedule(indx,:) = [];               % remove bus trips from airline schedule


%% Aircraft type information 

    [~,varnames] = xlsread(filename,4,'A1:D1');
    [~,~,AC]   = xlsread(filename,4,'A2:D5');
    AC = cell2table(AC,'VariableNames',varnames);

    K  = size(AC,1); % The total number of aircraft types

%% L: Flight arcs for each aircraft type k 

    L        = size(schedule.ORG,1); % number of flight arcs
    a        = cell(K,1); % airports in time-space network for aircraft type k
    
    for k = 1:K
        timespace(k).fl         = [schedule(:,1:5) schedule(:,5+k)]; 
        % the flight arc ends at the turn around time
        timespace(k).fl.Arrival = timespace(k).fl.Arrival + AC.TAT(k); 
        timespace(k).fl.Properties.VariableNames{end} = 'Cost'; 
        
        %% Use this part of the code only if we decide to eliminate flights 
        % that cannot be flown by a specific aircraft type from the time
        % space network of that flight
        timespace(k).fl((strcmp(timespace(k).fl.Cost,'NA')),:) = [];
        a{k,1}(:)    = unique([timespace(k).fl.ORG ; timespace(k).fl.DEST ]);
    end

%% Find airports in time-space network    
    % If original approach determine the locations ( airports) in the time 
    % space networks
    
    %a  = unique(schedule.ORG); % airport names
    %A         = size(a,1);            % number of airports        
    

%% Time space network for each aircraft type k
  
    for k = 1:K
        timespace(k).node       = {};
        timespace(k).ga         = {};
        timespace(k).nga        = {};
        for j = 1:size(a{k,1},2)
            % index of nodes departing and arriving at airport j 
            % (there may be duplicate nodes at this point because multiple 
            % flight arcs can end at one node)        
            idx            = [strcmp(timespace(k).fl.ORG, a{k,1}(j)) ...
                            strcmp(timespace(k).fl.DEST, a{k,1}(j))] ; 
            
            % nodes at one airport j ordered in chronological time 
            %( duplicate nodes removed)
            node_j_time    = unique(cell2mat(sortrows([...
                    table2cell(timespace(k).fl(idx(:,1),4)); ...
                 table2cell( timespace(k).fl(idx(:,2),5))],1)),'sorted'); 

            % generate cell array of nodes at airport j with 
            % airport IATA code
            node_j_k       = cell(size(node_j_time,1),2);     
            node_j_k(:,1)  = {a{k,1}(j)};                            
            node_j_k(:,2)  = num2cell(node_j_time);

            % cell array that stores all the nodes for an aircraft type
            timespace(k).node      = [timespace(k).node; node_j_k]; 

            % ground arcs are just a pair of sequential nodes at an airport
            for n = 1:(size(node_j_k,1)-1)    
                timespace(k).ga = [timespace(k).ga ; ...
                    node_j_k(n,:) node_j_k(n+1,2)];
            end

            % night arcs are formed by the last entry in list of nodes 
            % (node_i_k) and the first entry
            timespace(k).nga   = [timespace(k).nga ; ...
                node_j_k(end,:) node_j_k(1,2)];

        end
        timespace(k).nga = cell2table(timespace(k).nga, ...
            'VariableNames',{'Loc','Departure','Arrival'});
        timespace(k).ga  = cell2table(timespace(k).ga,...
            'VariableNames',{'Loc','Departure','Arrival'});
    end

    disp('Number of ground arcs for A330: ') 
    disp(size(timespace(1).ga,1))
    disp('Number of overnight arcs for A330: ') 
    disp(size(timespace(1).nga,1))

%% ------------- END OF CODE ----------------------------------------------

