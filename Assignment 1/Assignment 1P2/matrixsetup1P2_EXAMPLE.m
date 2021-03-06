        %% EXAMPLE Input
function [P, R, L, fare, fare_r, demand, capacity, col, delta, Q, costfull, Bpr, ... 
         recap_p, recap_r, recaprate] = matrixsetup1P2_EXAMPLE(filename)  

% filename = 'Input_Example.xlsx';
   
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
    recapture.p     = xlsread(filename,2,'A14'); % recaptured from path p
    recapture.r     = xlsread(filename,2,'B14'); % recaptured to path r


        %% EXAMPLE Output INITIALISATION!!!!!!!!!!!!!!!!!!!!
   % What do we need
    P = numel(itinerary.paths); 
    L = numel(flight.capacity);  % ugly way of doing things
    
    fare   = itinerary.fare;    % Initially all can be recaptured to p=0, with fare_0 = 0
    recap_p      = recapture.p ;
    recap_r      = recapture.r ;
    recaprate    = recapture.rate ;
    pathflights  = itinerary.flightnr ;
    flightnrs    = flight.flightnr;
    demand       = itinerary.demand;
    capacity     = flight.capacity;
    
        %% Calculations
Bpr = xlsread('Bpr_EXAMPLE',1,'A1:I8') ;

Rcol = 9;
Pcol = 1:P;
col = zeros(P,2);
for r = Rcol
    for p = Pcol
        col(p,:) = [p r];
    end
end

% The binary value: Delta
delta = cell(P,1); % For each Path all Flights are checked: delta{p,1}(i)
for p = 1:P
    delta{p,1} = zeros(L,1);
    for i = 1:L
        if pathflights(p,1) == flightnrs(i) || pathflights(p,2) == flightnrs(i)
            delta{p,1}(i) = 1;            
        end
    end
end

delta{9,1} = zeros(L,1);

% The daily unconstrained demand on flight(i): Q
Q = zeros(L,1);
for i = 1:L
    a = zeros(P,1);
    for p = 1:P
        a(p) = demand(p)*delta{p,1}(i);
    end
    Q(i) = sum(a);
end

% adding the ficticious fare and element
fare_r = [fare;0] ; % adding the ficticious recapture fare
R = P +1;


costfull = zeros(P,R);      % full cost matrix
for p = 1:P
    for r = 1:R
        costfull(p,r) = fare(p) - Bpr(p,r)*fare_r(r);
    end
end
end






