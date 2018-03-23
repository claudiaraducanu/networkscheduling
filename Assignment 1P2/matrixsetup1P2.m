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
%    filename = 'Input_AE4424_Ass1P2.xlsx';

        %% EXAMPLE Input
function [P, L, fare, demand, capacity, pathflights, flightnrs] ... 
          = matrixsetup1P2_EXAMPLE(filename)  

   
   % All dependend on P (itineraries)
    [~, ditinerary.flightnr]    = xlsread(filename,2,'F2:G738');   % flight numbers both legs
    itinerary.flightnr          = ditinerary.flightnr; 
    itinerary.paths             = xlsread(filename,2,'H2:H738'); % itinerary number
    itinerary.fare              = xlsread(filename,2,'E2:E738'); % fare associated to itineraries
    itinerary.demand            = xlsread(filename,2,'D2:D738'); % demand associated to itineraries
    
   % All dependend on L (flights)
    [~, dflight.flightnr]       = xlsread(filename,1,'A2:A233'); % flight number 
    flight.flightnr             = dflight.flightnr;
    flight.capacity             = xlsread(filename,1,'F2:F233'); % capacity associated to existing arc
    
   % All about recapturing
%     recapture.rate     = xlsread(filename,3,'C2:C300'); % recapture rate associated to some arcs
%     recapture.p   = xlsread(filename,3,'G2:G300'); % recaptured from path p
%     recapture.r   = xlsread(filename,3,'H2:H300'); % recaptured to path r


        %% EXAMPLE Output INITIALISATION!!!!!!!!!!!!!!!!!!!!
   % What do we need
    P = numel(itinerary.paths); 
    L = numel(flight.capacity);  % ugly way of doing things
    
    fare   = itinerary.fare;    % Initially all can be recaptured to p=0, with fare_0 = 0
%     recap_p      = recapture.p ;
%     recap_r      = recapture.r ;
%     recaprate    = recapture.rate ;
    pathflights  = itinerary.flightnr ;
    flightnrs    = flight.flightnr;
    demand       = itinerary.demand;
    capacity     = flight.capacity;
end

%% Write Bpr file
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







    

