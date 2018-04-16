%function [idx_ofl,idx_ifl] = fik_idx(timespace,k,n)
%fik_idx - Identifies the inbound and outbound flights at each node in the 
%the timspace network
%converts it into a set of nodes and arcs that would represent the
%time-space network for each aircraft type the airline owns. 
%
% Inputs:
%    filename - .xls file containing the schedule of the airline 
%    
% Outputs:
%    idx_ofl  - idx of flights inbound to the node
%    idx_ifl  - idx of flights outbound from the node
% Functions: 
%    read_schedule
% Author: Group 4
% April 2018;
% Assignment 2, Network Scheduling
%% ------------- BEGIN CODE -----------------------------------------------
    node_name = timespace(k).node{n,1};
    node_time = timespace(k).node{n,2};

    % Outbound flights from node n in N^k
    idx_ofl = double(strcmp(node_name,timespace(k).fl.ORG)) + ...
            double(timespace(k).fl.Departure == node_time);
    idx_ofl = find(idx_ofl == 2);
    %O_kn   = timespace(k).fl.Flight(idx_ofl);

    % Inbound flights from node n in N^k
    idx_ifl = double(strcmp(node_name,timespace(k).fl.DEST)) + ...
            double(timespace(k).fl.Arrival == node_time);
    idx_ifl = find(idx_ifl == 2);
    %I_kn   = timespace(k).fl.Flight(idx_ifl);
%end
%% ------------- END OF CODE ----------------------------------------------