function [P, R, L, fare, fare_r, demand, col, delta, Q, farecost, Bpr, ... 
         recap_p, recap_r, recaprate] = setup2P1(filename)  

% This file contains the function to read most of the input file and some
% initial calculation. 
%
% Inputs:
%    filename - .xls file containing the information of the airline 
%   
% Outputs:
%    P          - number of passenger itineraries
%    R          - number of passenger itineraries + the ficitious
%    L          - number of flights in the schedule
%    fare       - fare for each itinerary
%    fare_r     - fare for the recaptured passengers from p to r
%    demand     - demand for each itinerary
%    col        - 2 column matrix, containing all columns considered in the
%                 problem
%    delta      - cell array, with P cells, containing L elements. Each
%                 element is 1 if flight is in that itinerary and is 0 
%                 otherwise.
%    Q          - the daily unconstrained demand on flight i
%    costfull   - PxR matrix, containing all costs for each O&D combination
%    Bpr        - PxR matrix, containing the recapture rates for the p's 
%                 that can be recaptured by r (written in file: Bpr.xls)
%    recap_p    - list of all itineraries that can be recaptured
%    recap_r    - list of all itineraries that can recapture
%    recaprate  - list of the recapturerates (same order as recap_p/r)
% Author: Group 4
% April 2018;
% Assignment 2, Network Scheduling

%     clearvars
% %    clear all
%     filename = 'Assignment2.xlsx';
       %% Read file for network   
   % All dependend on P (itineraries)
    [~, ditinerary.flightnr]    = xlsread(filename,2,'F2:G738');   % flight numbers both legs
    itinerary.flightnr          = ditinerary.flightnr; 
    itinerary.paths             = xlsread(filename,2,'H2:H738'); % itinerary number
    itinerary.fare              = xlsread(filename,2,'E2:E738'); % fare associated to itineraries
    itinerary.demand            = xlsread(filename,2,'D2:D738'); % demand associated to itineraries
    
   % All dependend on L (flights)
    [~, dflight.flightnr]       = xlsread(filename,5,'A2:A233'); % flight number 
    flight.flightnr             = dflight.flightnr;
    
   % All about recapturing
    recapture.rate = xlsread(filename,3,'C2:C300'); % recapture rate associated to some arcs
    recapture.p   = xlsread(filename,3,'G2:G300'); % recaptured from path p
    recapture.r   = xlsread(filename,3,'H2:H300'); % recaptured to path r


        %% Input
   % What do we need
    P = numel(itinerary.paths); 
    L = numel(flight.flightnr); 
    
    fare         = itinerary.fare;    % Initially all can be recaptured to p=0, with fare_0 = 0
    recap_p      = recapture.p ;
    recap_r      = recapture.r ;
    recaprate    = recapture.rate ;
    pathflights  = itinerary.flightnr ;
    flightnrs    = flight.flightnr;
    demand       = itinerary.demand;
    
    %% CALCULATIONS
    Bpr = xlsread('Bpr',1,'A1:ABJ738') ; % Bpr is printed in a file (See code below)

    % Making the initial list with all columns used. First column is
    % initial path p (first column), the recaptured path r (second column)
    Rcol = 738;
    Pcol = 1:P;
    col = zeros(P,2);
    for r = Rcol
        for p = Pcol
            col(p,:) = [p r]; % Add all the initial columns.
        end
    end


    % The binary value: Delta
    delta = cell(P,1); % For each Path p it is checked if the Flights i is used: delta{p,1}(i)
    for p = 1:P
        delta{p,1} = zeros(L,1);
        for i = 1:L
            if isequal(pathflights(p,1), flightnrs(i)) || ...
                    isequal(pathflights(p,2), flightnrs(i))
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


    farecost = zeros(P,R);  % matrix containing costs of all p,r combinations
    for p = 1:P
        for r = 1:R
            farecost(p,r) = fare(p) - Bpr(p,r)*fare_r(r);
        end
    end

end


% %% Write Bpr file
%  % Bpr is a P by R matrix containing the recapture rates for the p's that
%  % can be recaptured by r
% 
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







    

