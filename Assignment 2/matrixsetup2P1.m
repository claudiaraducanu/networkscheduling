%function [nodes, fare, demand , capacity, origin , dest, recapture, s, ...
%           t] = matrixsetup(filename)
       %% Read file for network

    % Not used
%     [~, ditinerary.origin]  = xlsread(filename,2,'B2:B738'); % names of departure airports
%     [~, ditinerary.dest]    = xlsread(filename,2,'C2:C738'); % names of arrival airports   

%     [~, dflight.origin]  = xlsread(filename,2,'B2:B738'); % names of departure airports
%     [~, dflight.dest]    = xlsread(filename,2,'C2:C738'); % names of arrival airports    
   
    
        %% EXAMPLE Input

%    clearvars
%    clear all
    filename = 'Assignment2.xlsx';

        %% EXAMPLE Input
%function [P, R, L, fare, fare_r, demand, capacity, col, delta, Q, costfull, Bpr, ... 
%          recap_p, recap_r, recaprate] = matrixsetup2P1(filename)  

   
   % All dependend on P (itineraries)
    [~, ditinerary.flightnr]    = xlsread(filename,2,'F2:G738');   % flight numbers both legs
    itinerary.flightnr          = ditinerary.flightnr; 
    itinerary.paths             = xlsread(filename,2,'H2:H738'); % itinerary number
    itinerary.fare              = xlsread(filename,2,'E2:E738'); % fare associated to itineraries
    itinerary.demand            = xlsread(filename,2,'D2:D738'); % demand associated to itineraries
    
   % All dependend on L (flights)
    [~, dflight.flightnr]       = xlsread(filename,1,'A2:A233'); % flight number 
    flight.flightnr             = dflight.flightnr;
    %flight.capacity             = xlsread(filename,1,'F2:F233'); % capacity associated to existing arc
    
   % All about recapturing
    recapture.rate = xlsread(filename,3,'C2:C300'); % recapture rate associated to some arcs
    recapture.p   = xlsread(filename,3,'G2:G300'); % recaptured from path p
    recapture.r   = xlsread(filename,3,'H2:H300'); % recaptured to path r


        %% EXAMPLE Output INITIALISATION!!!!!!!!!!!!!!!!!!!!
   % What do we need
    P = numel(itinerary.paths); 
    L = numel(flight.flightnr);  % ugly way of doing things
    
    fare         = itinerary.fare;    % Initially all can be recaptured to p=0, with fare_0 = 0
    recap_p      = recapture.p ;
    recap_r      = recapture.r ;
    recaprate    = recapture.rate ;
    pathflights  = itinerary.flightnr ;
    flightnrs    = flight.flightnr;
    demand       = itinerary.demand;
    %capacity     = flight.capacity;
    
    %% CALCULATIONS
    Bpr = xlsread('Bpr',1,'A1:ABJ738') ;

    Rcol = 738;
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
            if isequal(pathflights(p,1), flightnrs(i)) || isequal(pathflights(p,2), flightnrs(i))
                delta{p,1}(i) = 1;            
            end
        end
    end

    delta{738,1} = zeros(L,1);

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

%end
% %% Write Bpr file
% ptor = zeros(P);
% for p = 1:P
%     for r = 1:P
%         for laa = 1:numel(recap_p)
%             if p == recap_p(laa) && r == recap_r(laa)
%                 ptor(p,r) = recaprate(laa);
%             end
%         end
%     end
% end
% 
% xlswrite('Bpr.xlsx',ptor)







    

