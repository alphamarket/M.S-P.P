function [pop, desc] = lem_learn_mode(pop, learn_probe, learn_threshold, elit_ratio_func, save_history_func)
% LEM_LEARN_MODE
% @brief                        The executor of machine learning mode in LEM algorithm
% @param  pop                  	The input population
% @param  learn_probe           The learning probe
% @param  learn_threshold       The learning threshold
% @param  elit_ratio_func       The function to fetch elit selection's ratio in every generation
% @param  save_history_func     The function for saving history in every generation
% @return                       The result pop from machine learning evolution
    disp('Machine Learning mode');
    fprintf('%-12sgen#%-10smax(LP)%-10smin(LP)%-10smax(HP)%-10smin(HP)\n', ' ', ' ', ' ', ' ', ' ');
    for genNO=1:learn_probe
        % sort pop
        pop = sortrows(pop, size(pop, 2));
        % fetch {H|L}-group
        [HP, LP] = mlm_fetchXgroups(pop, elit_ratio_func(genNO));
        % print stat
        mlm_stat(pop, genNO - 1, HP, LP);
        % save pop history
        save_history_func(pop);
        % check learn-threshold
        if HP(1,end) <= learn_threshold, break; end
        % fetch description => do elitism + new population
        pop((size(HP, 1)+1):end,:) = mlm_gen(pop, HP, mlm_getDesc(pop, HP, LP, 'nn'), genNO);
    end
    % if learn_probe < 1
    if isempty(genNO), genNO = 1; end
    % sort pop
    pop = sortrows(pop, size(pop, 2));
    % save pop history
    save_history_func(pop);
    % fetch {H|L}-group
    [HP, LP] = mlm_fetchXgroups(pop, elit_ratio_func(genNO));
    % print stat
    mlm_stat(pop, genNO, HP, LP);
    % return the final description
    desc = mlm_getDesc(pop, HP, LP, 'nn');
end