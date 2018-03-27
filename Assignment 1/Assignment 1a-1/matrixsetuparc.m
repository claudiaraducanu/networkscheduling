function [nodes, commodities, arcs, origin , dest, demand, ntw] = matrixsetuparc(filename,undirected)

    %% Read file for network
   
    [~, dnetwork.origin]   = xlsread(filename,1,'B2:B31'); % names of departure airports
    [~, dnetwork.dest]     = xlsread(filename,1,'C2:C31'); % names of arrival airports
    dnetwork.acost         = xlsread(filename,1,'D2:D31'); % cost associated to existing arc
    dnetwork.acap          = xlsread(filename,1,'E2:E31'); % cost associated to existing arc
    %% Turn graph into undirected graph
    % Since the fight routes are operated daily in both directions.
    
    if undirected == 1
    
        network.data.origin        = [ dnetwork.origin; dnetwork.dest];   % destinations are also origins
        network.data.dest          = [ dnetwork.dest;   dnetwork.origin]; % origins are also destinations
        network.data.cost          = [ dnetwork.acost; dnetwork.acost];  
        network.data.cap           = [ dnetwork.acap; dnetwork.acap];
    else 
        network.data.origin        = dnetwork.origin;  
        network.data.dest          = dnetwork.dest;  
        network.data.cost          = dnetwork.acost;
        network.data.cap           = dnetwork.acap; 
        
    end
    %% Determine O-D pair cost

    network.gcost       = digraph(network.data.origin,network.data.dest,...
                            network.data.cost);
    network.gcost.Nodes.Number = (1:size(network.gcost.Nodes.Name,1))'; 
    nodes               = numnodes(network.gcost);
    [s,t]               = findedge(network.gcost);
    
    %% Determine O-D pair capacity

    network.gcap        = digraph(network.data.origin,network.data.dest,...
                            network.data.cap);
    capacity            = network.gcap.Edges.Weight(findedge(network.gcap,s,t));
    network.gcost.Edges.capacity = capacity;
      
    arcs                = numedges(network.gcost);
    
    ntw                 = network.gcost;
    
    %% Airport names and numbers
    airports =  network.gcost.Nodes;               

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

                

