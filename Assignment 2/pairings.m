clearvars
clear all
%% Input: Flights operated by B737-700
filename1                   = 'B737.txt';
filename2                   = 'B738.txt';


flight.fl7                  = readtable(filename1); 
flight.fl8                  = readtable(filename2);

%% Requirement on minimum sit connection between flights
flight.fl7.Arrival          = flight.fl7.Arrival - 30 ; 
flight.fl8.Arrival          = flight.fl8.Arrival - 35 ;

flight.fl                   = sortrows( [flight.fl7; flight.fl8],[2,4]);
flight.fl.Arrival(flight.fl.Arrival < 0) = flight.fl.Arrival(flight.fl.Arrival < 0)+ 1440; 
%% Determine node(airport,time)
% determine airports and number of airports in flight network
a         = unique(flight.fl.ORG);           % airport names
A         = size(a,1);                       % number of airports 
 
%% Determine sit connections 
pairs  = cell(size(flight.fl.ORG,1),1);
for s = 1:size(flight.fl.ORG,1)
    S = flight.fl(s,:);        
    req_1 = strcmp(S.DEST,flight.fl.ORG);
    % Requirement on maximum idle time between flights
    max_start = S.Arrival+180;
    min_start = S.Arrival+45;
    if max_start >= 1440
        max_start = max_start - 1440;
    end
    if min_start >= 1440
        min_start = min_start - 1440;
    end  
    req_3 = double(min_start  <= flight.fl.Departure );
    req_4 = double(max_start  >= flight.fl.Departure );
    % combine requirements 1 and 4 to determine next possible flights
    next  = flight.fl(req_1 + req_3 +req_4 == 3,:);
    if size(next,1) == 1 
        pairs{s,1} = [ S.Flight next.Flight];
    else
        ss       = cell(size(next,1),1);
        ss(:)    = S.Flight;
        pairs{s,1} = [ss next.Flight];
    end
    if size(next,1) == 0 && strcmp(S.DEST,'AEP')==1
        pairs{s,1} = 'remove';
    end
end


%% Determine ground arcs. 
            