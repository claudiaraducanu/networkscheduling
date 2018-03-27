function [nodes, commodities, arcs , origin , dest, demand, s1, ...
           t1,capacity] = read_arc(filename)
    %% Read file for network
%     clear all
%     clearvars 
%     
%     filename = 'Input_AE4424_Ass1P1.xlsx';
%     
    [~, net.data.origin]   = xlsread(filename,1,'B2:B31'); % names of departure airports
    [~, net.data.dest]     = xlsread(filename,1,'C2:C31'); % names of arrival airports
    network.data.cost         = xlsread(filename,1,'D2:D31'); % cost associated to existing arc
    network.data.cap          = xlsread(filename,1,'E2:E31'); % cost associated to existing arc
   
  
    network.data.origin   = [ net.data.origin ; net.data.dest];  % names of departure airports
    network.data.dest     = [ net.data.dest ; net.data.origin]; % names of arrival airports
    network.data.cost     = [ network.data.cost; network.data.cost]; % cost associated to existing arc
    network.data.cap      = [ network.data.cap ; network.data.cap]; % cost associated to existing arc
   
    
    %% Determine O-D pair cost

    network.gcost       = digraph(network.data.origin,network.data.dest,...
                            network.data.cost);
    network.gcost.Nodes.Number = (1:size(network.gcost.Nodes.Name,1))'; 
    nodes               = numnodes(network.gcost);
    
    %% 
    [s1,t1]               = findedge(network.gcost);
    %cost                  = network.gcost.Edges.Weight(findedge(network.gcost,s1,t1));
   
    %% Determine O-D pair capacity

    network.gcap        = digraph(network.data.origin,network.data.dest,...
                            network.data.cap);
    capacity            = network.gcap.Edges.Weight(findedge(network.gcap,s1,t1));
    %network.gcost.Edges.Capacity = capacity;
    
    %% Airport names and numbers
    airports =  network.gcost.Nodes;               
    arcs     = network.gcost;
    
    %% Read file for commodities

    [~, cargo.origin]   = xlsread(filename,2,'B2:B41'); % names of origin airport for commodities
    [~, cargo.dest]     = xlsread(filename,2,'C2:C41'); % names of destination airports for commodities
    demand              = xlsread(filename,2,'D2:D41'); % demand for each commodity

    demand              = transpose(demand);

    commodities         = size(cargo.origin,1);
    
    origin = zeros(1,commodities);
    dest = zeros(1,commodities);

    for i= 1:commodities
        for j = 1:nodes
            if strcmp(airports.Name(j),cargo.origin(i)) == 1
                origin(i) = airports.Number(j);
            end
            if strcmp(airports.Name(j),cargo.dest(i)) == 1
                dest(i) = airports.Number(j);
            end
        end
    end

end

                

