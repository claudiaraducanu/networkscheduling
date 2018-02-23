%% Initialize
clearvars
clear all

%% Read file for network
filename     = 'Input_AE4424_Ass1P1.xlsx';
[~, dnetwork.origin]   = xlsread(filename,1,'B2:B31'); % names of departure airports
[~, dnetwork.dest]     = xlsread(filename,1,'C2:C31'); % names of arrival airports
dnetwork.acost         = xlsread(filename,1,'D2:D31'); % cost associated to existing arc
dnetwork.acap          = xlsread(filename,1,'E2:E31'); % cost associated to existing arc
 
%% Turn graph into undirected graph
% Since the fight routes are operated daily in both directions.

network.data.origin        = [ dnetwork.origin; dnetwork.dest];
network.data.dest          = [ dnetwork.dest;   dnetwork.origin];
network.data.cost          = [ dnetwork.acost; dnetwork.acost];
network.data.cap           = [ dnetwork.acap; dnetwork.acap];

%% Determine O-D pair cost

network.gcost       = digraph(network.data.origin,network.data.dest,...
                        network.data.cost);
network.gcost.Nodes.Number = (1:16)'; 
nodes               = numnodes(network.gcost);
[s,t]               = findedge(network.gcost);
network.mcost       = full(sparse(s,t,network.gcost.Edges.Weight,...
                    nodes,nodes));
indx                = find(network.mcost == 0);
network.mcost(indx) = 1000;

%% Determine O-D pair capacity

network.gcap        = digraph(network.data.origin,network.data.dest,...
                        network.data.cap);
[s,t]               = findedge(network.gcap);
network.mcap        = full(sparse(s,t,network.gcap.Edges.Weight,...
                    nodes,nodes));

%% Read file for commodities


                

