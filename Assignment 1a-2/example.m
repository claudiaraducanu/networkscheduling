    


    param.delta = cell(K,1);
    
    for k=1:K
        param.delta{k,1} = zeros(P_k(k),A);
        for p = 1:P_k(k)
            for a = 1:A
                pathlength = size(sp.path{k,1}{p,1},2);
                if pathlength > 2 
                    ainp = repelem(sp.path{k,1}{p,1},2);
                    ainp = ainp(2:end-1);
                    ainp = transpose(reshape(ainp,2,pathlength-1));
                    for ac = 1:(pathlength-1) 
                        if sum(ainp(ac,:) == sp.arcs{a,1}) == 2
                            param.delta{k,1}(p,a)  = 1;
                        end
                    end
                else
                    if sum(sp.path{k,1}{p,1} == sp.arcs{a,1}) == 2
                     param.delta{k,1}(p,a)  = 1;
                    end
                end
            end
         end
    end 