function pop = init(pop_size, lchrom, genes_margin)
% INIT
% @brief                Initializes a population
% @param  pop_size      The population's size
% @param  lchrom        The lenght of chromosomes
% @param  genes_margin  The maximum range of each genes
% @return               The generated popualtion with calculated fitness
    pop = gen_genes(pop_size, lchrom, genes_margin); 
    pop(:, end) = fitfunc(pop);    
end