function d = gen_genes(size, length, margin)
% GEN_DATA
% @brief                Generates new individual(s)
% @param  size          The size of population
% @param  length        The lenght of each chromosomes
% @param  margin   		The maximum range of each genes
% @return               The generated individual(s)
    d = -margin + rand(size, length) * 2 * margin;
end