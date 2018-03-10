function [nodes, commodities, cost , capacity, origin , dest, demand, s, ...
           t] = read_arc_v(filename)
    %% Read file for network
    
    [~, network.data.origin]   = xlsread(filename,1,'B2:B8'); % names of departure airports
    [~, network.data.dest]     = xlsread(filename,1,'C2:C8'); % names of arrival airports
    network.data.cost         = xlsread(filename,1,'D2:D8'); % cost associated to existing arc
    network.data.cap          = xlsread(filename,1,'E2:E8'); % cost associated to existing arc
   
    %% Determine O-D pair cost

    network.gcost       = digraph(network.data.origin,network.data.dest,...
                            network.data.cost);
    network.gcost.Nodes.Number = (1:size(network.gcost.Nodes.Name,1))'; 
    nodes               = numnodes(network.gcost);
    
    %% 
    [s,t]               = findedge(network.gcost);
    cost                = sparse(s,t,network.gcost.Edges.Weight,...
                        nodes,nodes);

    %% Determine O-D pair capacity

    network.gcap        = digraph(network.data.origin,network.data.dest,...
                            network.data.cap);
    [s,t]               = findedge(network.gcap);
    capacity            = sparse(s,t,network.gcap.Edges.Weight,...
                        nodes,nodes);

    %% Airport names and numbers
    airports =  network.gcost.Nodes;               

    %% Read file for commodities

    [~, cargo.origin]   = xlsread(filename,2,'B2:B5'); % names of origin airport for commodities
    [~, cargo.dest]     = xlsread(filename,2,'C2:C5'); % names of destination airports for commodities
    demand              = xlsread(filename,2,'D2:D5'); % demand for each commodity


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

                

