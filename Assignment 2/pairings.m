clearvars
clear all
%% Input: Flights operated by B737-700
filename                   = 'B737.txt';
flight.fl                  = readtable(filename);

% For daily problems each flight arc is replicated as many times as the 
% maximum number of calendar days allowed in a pairing but pairings are 
% generated only from flights operating on the first day.
%% Duplicate flight arcs 1x

day2 = flight.fl;
day2.Departure  = day2.Departure + 24*60;
day2.Arrival    = day2.Arrival   + 24*60;
flight. fl      = sortrows([ flight.fl ; day2],[2 ,1]); % concatinate 2 tables 

%% Determine node(airport,time)
% 
% % determine airports and number of airports in flight network
% a         = unique(timespace.fl.ORG);        % airport names
% A         = size(a,1);                       % number of airports 
% 
% timespace.node = {};
% for j = 1:A
%     % index of nodes in table of flights departing and arriving at airport j    
%     % (there may be duplicate nodes at this point because multiple 
%     % flight arcs can end at one node)     
%     idx            = [strcmp(timespace.fl.ORG, a(j)) ...
%                     strcmp(timespace.fl.DEST, a(j))] ; 
%     % nodes at one airport j ordered in chronological time 
%     %( duplicate nodes removed)                
%     node_j_time    = unique(sortrows([...
%             timespace.fl.Departure(idx(:,1)); ...
%             timespace.fl.Arrival(idx(:,2))],1),'sorted');                 
%     % generate cell array of nodes at airport j with 
%     % airport IATA code
% 
%     node_j_k       = cell(size(node_j_time,1),2);     
%     node_j_k(:,1)  = {a(j)};                            
%     node_j_k(:,2)  = num2cell(node_j_time);        
%     % cell array that stores all the nodes for an aircraft type
%     timespace.node      = [timespace.node; node_j_k]; 
% end
% timespace.node = cell2table(timespace.node, ...
%             'VariableNames',{'airport','time'});
% %duties.source = timespace.node(strcmp(timespace.node.airport, 'AEP'),:); 

%% Determine ground arcs. 
            