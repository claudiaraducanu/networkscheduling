        %% EXAMPLE Input
function [P, L, fare, recaprate, recap_p, recap_r, ... 
          demand , capacity, pathflights, flightnrs] = matrixsetup1P2_EXAMPLE(filename)  

 %  filename = 'Input_Example.xlsx';
   
   % All dependend on P (itineraries)
    itinerary.flightnr          = xlsread(filename,2,'F2:G9');   % flight numbers both legs
    itinerary.paths             = xlsread(filename,2,'A2:A9'); % itinerary number
    itinerary.fare              = xlsread(filename,2,'E2:E9'); % fare associated to itineraries
    itinerary.demand            = xlsread(filename,2,'D2:D9'); % demand associated to itineraries
    
   % All dependend on L (flights)
    flight.flightnr             = xlsread(filename,1,'A2:A7'); % flight number 
    flight.capacity             = xlsread(filename,1,'F2:F7'); % capacity associated to existing arc
    
   % All about recapturing
    recapture.rate  = xlsread(filename,2,'C14'); % recapture rate associated to some arcs
%     recapture.p     = xlsread(filename,2,'A14'); % recaptured from path p
%     recapture.r     = xlsread(filename,2,'B14'); % recaptured to path r
    recapture.pfare = xlsread(filename,2,'D14'); % fare of original p
    recapture.rfare = xlsread(filename,2,'E14'); % fare of recapture r

        %% EXAMPLE Output INITIALISATION!!!!!!!!!!!!!!!!!!!!
   % What do we need
    P = numel(itinerary.paths); 
    L = numel(flight.capacity);  % ugly way of doing things
    
    fare   = itinerary.fare;    % Initially all can be recaptured to p=0, with fare_0 = 0
%     fare_p       = recapture.pfare ;
%     fare_r       = recapture.rfare ;
    recap_p      = recapture.p ;
    recap_r      = recapture.r ;
    recaprate    = recapture.rate ;
    pathflights  = itinerary.flightnr ;
    flightnrs    = flight.flightnr;
    demand       = itinerary.demand;
    capacity     = flight.capacity;
end
%     

% Can p be recaptured?
% recapture = ones(P,2)*1000
% for p = 1:P
%     if p
% end

%% TO DO
% Make matrix for bpr (0 if not)
% Only need fare_p_all




    

