%function [nodes, fare, demand , capacity, origin , dest, recapture, s, ...
%           t] = matrixsetup(filename)
       %% Read file for network
   
    [~, dnetwork.originP]  = xlsread(filename,2,'B2:B738'); % names of departure airports
    [~, dnetwork.destP]    = xlsread(filename,2,'C2:C738'); % names of arrival airports
    dnetwork.fare          = xlsread(filename,2,'E2:E738'); % fare associated to all arcs
    dnetwork.demand        = xlsread(filename,2,'D2:D738'); % demand associated to all arcs
    dnetwork.capacity      = xlsread(filename,1,'F2:F233'); % capacity associated to existing arc
    dnetwork.recapture     = xlsread(filename,3,'C2:C300'); % recapture rate associated to some arcs
