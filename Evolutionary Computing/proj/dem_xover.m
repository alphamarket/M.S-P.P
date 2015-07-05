function [c1, c2] = dem_xover(pop, gao)
% DEM_XOVER
% @brief        Xovers two individuals from population with respect to a GA option
% @param  pop  	The population to make the children from
% @param  gao   The genetic option
% @return       Two children

    % randomly choose 2 individuals
    n = size(pop, 1);
    l = size(pop, 2);
    % select parents
    parents = pop(randperm(n,2),:);
    % choose a crossover point
    cp = randperm(l, 1);
    % should not perform xover?
    if rand < (1 - gao.xover_prob), c1 = parents(1, :); c2 = parents(2, :); return; end
    % crossover
    c1 = [parents(1, 1:cp), parents(2, cp+1:end)];
    c2 = [parents(2, 1:cp), parents(1, cp+1:end)];
end