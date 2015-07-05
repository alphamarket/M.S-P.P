function m = dem_mutation(individuals, gao, new_gen_func)
% DEM_MUTATION
% @brief                Mutates a collection of individuals
% @param  individuals   The individuals to mutate
% @param  gao           The genetic option
% @param  new_gen_func  The new gene generator
% @return               The mutated individuals
    m = zeros(size(individuals));
    for i=1:size(individuals, 1)
        m(i, :) = individuals(i, :); 
        if rand < (1 - gao.mutation_prob), continue; end
        m(i, randperm(size(m(i, :), 2), 1)) = new_gen_func();
    end
end