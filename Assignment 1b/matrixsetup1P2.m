% %function [nodes, fare, demand , capacity, origin , dest, recapture, s, ...
% %           t] = matrixsetup(filename)
%        %% Read file for network
%    
%     [~, ditinerary.origin]  = xlsread(filename,2,'B2:B738'); % names of departure airports
%     [~, ditinerary.dest]    = xlsread(filename,2,'C2:C738'); % names of arrival airports
%     [~, ditinerary.flightnr]= xlsread(filename,2,'F2:F738'); % flight number Leg 1
%     ditinerary.fare         = xlsread(filename,2,'E2:E738'); % fare associated to all arcs
%     ditinerary.demand       = xlsread(filename,2,'D2:D738'); % demand associated to all arcs
%     
%     [~, dflight.origin]  = xlsread(filename,2,'B2:B738'); % names of departure airports
%     [~, dflight.dest]    = xlsread(filename,2,'C2:C738'); % names of arrival airports
%     [~, dflight.flightnr]= xlsread(filename,2,'F2:F738'); % flight number 
%     dflight.capacity     = xlsread(filename,1,'F2:F233'); % capacity associated to existing arc
%     
%     dpath.recapture     = xlsread(filename,3,'C2:C300'); % recapture rate associated to some arcs
%        %% Trying to figure it out
%     % itenaries
%     itinerary.data.origin     = ditinerary.origin;
%     itinerary.data.dest       = ditinerary.dest;
%     itinerary.data.flightnr   = ditinerary.flightnr; %How to do it, 2 flight legs?
%     itinerary.data.fare       = ditinerary.fare;
%     itinerary.data.demand     = ditinerary.demand;
%     
%     flight.data.origin      = dflight.origin;
%     flight.data.dest        = dflight.dest;
%     flight.data.flightnr    = dflight.flightnr;
    
   
    
    
    
    
    
    
    
    
    
    

