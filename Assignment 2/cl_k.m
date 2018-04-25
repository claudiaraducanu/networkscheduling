function [idx_fls,no] = cl_k(timespace,k)
    idx_fls = find(timespace(k).fl.Departure > timespace(k).fl.Arrival); 
    no      = size(idx_fls,1);
end