function fitness = fitfunc(pop)
% FITFUNC
% @brief        Calculates the fitness value
% @param  pop  	The population or individual to calc. the fitness for
% @return       The fitness value(s)
    fitness = zeros(size(pop, 1), 1);
    for i=1:size(pop, 1)
        fitness(i) = ackley(pop(i, :), 20, 0.2, 2 * pi);
    end
end