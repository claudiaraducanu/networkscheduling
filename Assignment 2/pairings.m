clearvars
clear all
%% Input: Flights operated by B737-700
filename1                   = 'B737.txt';
%filename2                   = 'B738.txt';


flight.fl7                  = readtable(filename1); 
flight.fl8                  = readtable(filename2);


flight.fl7.Arrival          = flight.fl7.Arrival - 30; 
flight.fl8.Arrival          = flight.fl8.Arrival - 35;
flight.fl                   = sortrows( [flight.fl7; flight.fl8],[2,4]);
% For daily problems each flight arc is replicated as many times as the 
% maximum number of calendar days allowed in a pairing but pairings are 
% generated only from flights operating on the first day.
%% Duplicate flight arcs 1x

% day2 = flight.fl;
% day2.Departure  = day2.Departure + 24*60;
% day2.Arrival    = day2.Arrival   + 24*60;
% flight. fl      = sortrows([ flight.fl ; day2],[2 ,4]); % concatinate 2 tables 
%% Determine node(airport,time)

% determine airports and number of airports in flight network
a         = unique(flight.fl.ORG);          % airport names
A         = size(a,1);                       % number of airports 
 
%% Determine sit connections
source = flight.fl(strcmp(flight.fl.ORG, 'AEP'),:); 
pairs  = cell(size(source,1),1);
for s = 1:size(source,1)
    S = source(s,:);
    req_1 = double(strcmp(S.DEST,flight.fl.ORG));
    req_3 = double(S.Arrival+45   <= flight.fl.Departure );
    req_4 = double(S.Arrival+3*60  >= flight.fl.Departure );
    next  = flight.fl(req_1 + req_3 + req_4 == 3,:);
    if size(next,1) == 1 
        pairs{s,1} = [ S.Flight next.Flight];
    else
        ss       = cell(size(next,1),1);
        ss(:)    = S.Flight;
        pairs{s,1} = [ss next.Flight];
    end
end


%% Determine ground arcs. 
            