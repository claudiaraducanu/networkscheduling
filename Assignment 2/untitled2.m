clear all
clearvars

filename              = 'Assignment2.xlsx';
[~,~,timespace]       = read_schedule(filename);
a = 1;

    GA  = zeros(1,size(timespace,2));
    for k  = 1:size(timespace,2)
        GA(k) = size(timespace(k).gat,1);
    end
    GA_k    = zeros(size(GA));   
    for k = 2:size(GA,2)
        GA_k(k) = sum(GA(1:k-1));
    end
    out = GA_k(k) + a; 