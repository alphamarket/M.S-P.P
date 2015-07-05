function pop = lem_darwinian_mode(pop, dar_probe, dar_threshold, ga_option, genes_margin, elit_ratio_func, save_history_func)
% LEM_DARWINIAN_MODE
% @brief                        The executor of darwinian evolutionary mode in LEM algorithm
% @param  pop                  	The input population
% @param  dar_probe             The darwinian probe
% @param  dar_threshold         The darwinian threshold
% @param  ga_option             The ordinary GA options
% @param  genes_margin          The range of valid genes in to generate new genes
% @param  elit_ratio_func       The function to fetch elit selection's ratio in every generation
% @param  save_history_func     The function for saving history in every generation
% @return                       The result pop from darwinian evolution
    disp('Darwinian Evolutionary mode');
    fprintf('%-12sgen#%-10smax(fit.)%-10savg(fit.)%-10smin(HP)\n', ' ', ' ', ' ', ' ');
    for genNO=1:dar_probe
        % sort pop and fetch elits
        [elits, pop] = dem_elit_selection(pop, elit_ratio_func(genNO));
        % print stat
        dem_stat(pop, genNO - 1);
        % save pop history
        save_history_func(pop);
        % check dar-threshold
        if elits(1, end) <= dar_threshold, break; end
        % create new population
        new_pop = elits;
        % forever
        while true
            % fetch children
            [c1, c2]    = dem_xover(pop,  ga_option);
            % mutate them
            cc          = dem_mutation([c1; c2], ga_option, @() gen_genes(1, 1, genes_margin));
            % calc. fitness of them
            cc(:, end)  = fitfunc(cc);
            % put into new population
            new_pop    = [new_pop; cc]; %#ok<AGROW>
            % continue while size of new population is less than current population
            if size(new_pop, 1) >= size(pop, 1), break; end
        end
        % fail-safe for population size overflow
        if size(new_pop, 1) > size(pop, 1), new_pop = new_pop(1:size(pop), :); end
        assert(isequal(size(pop), size(new_pop)), 'The size didn\''t match');
        % swap populations
        pop = new_pop;
    end
    % if dar_probe < 1
    if isempty(genNO), genNO = 1; end
    % sort pop and fetch elits
    [~, pop] = dem_elit_selection(pop, elit_ratio_func(genNO));
    % print stat
    dem_stat(pop, genNO);
    % save pop history
    save_history_func(pop);
end