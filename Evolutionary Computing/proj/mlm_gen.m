function y = mlm_gen(pop, HP, desc, genNO) %#ok<*AGROW>
% MLM_GEN
% @brief        Generates new popultion based on a description
% @param  pop  	Current population
% @param  HP    The current population's H-group
% @param  desc  The current populaiton's description
% @return       The new population with calculated fitness values based on the given description

    % method: replace non-H-group with new individuals
    y = zeros(size(pop, 1) - size(HP, 1), size(pop, 2));
    for i=1:size(y, 1)
        y(i,1:(size(pop, 2) - 1)) = desc.instantize(i, size(y, 1), genNO);
        y(i,end) = fitfunc(y(i,:));
    end
    % select from generated population, and with 80% prob. initiate a 
    % one-point XOVER, custom mutation
    k = [];
    while(size(k, 1) < size(y, 1))
        [o, b] = mlm_tournament_selection(y);
        if(randn < 0.8)
            cp = randperm(size(o, 2) - 2, 1);
            k1 = [o(1:cp), b(cp+1:end-1), 0];
            k2 = [o(cp+1:end-1), b(1:cp), 0];
        else
            k1 = o; k2 = b;
        end
        for i=1:2
            if(randn < 0.1)
                l = desc.instantize(size(k, 1), size(y, 1), genNO);
                assert(isequal(size(k1), size(k2)), 'The size didn\''t match');
                p = randperm(size(k1, 2) - 1, 1);
                if(i==1), k1(p) = l(p); 
                else k2(p) = l(p); end
            end
        end
        k(end+1, :) = k1; 
        k(end, end) = fitfunc(k(end,1:end-1));
        k(end+1, :) = k2;
        k(end, end) = fitfunc(k(end,1:end-1));
    end
    if size(k, 1) > size(y, 1), k = k(1:size(y), :); end
    assert(isequal(size(k), size(y)), 'The size didn\''t match');
    y = k;
end

function [x1, x2] = mlm_tournament_selection(pop)
    X = pop(randperm(length(pop), ceil(2*sqrt(size(pop, 1)))), :);
    X = sortrows(X, size(X, 2));
    x1 = X(1,:); x2 = X(2, :);
end