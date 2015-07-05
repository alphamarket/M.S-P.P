function [pop, h, desc] = lem()
% LEM
% @brief    Executes the LEM algorithm
% @return   Returns the final population and a history of populations in evolutions

    %
    % intialize parameters
    %
    global hist
    global elit_ratio
    format shortEng
    % general params.
    elit_ratio          = 0.3;
    genes_margin        = 40;
    pop_size            = 100;
    lchrom              = 2 + 1;
    % machine learning mode's param
    learn_probe         = 100;
    learn_threshold     = 01e-15;
    % darwinian evolutionary mode's params
    dar_probe           = 000;
    dar_threshold       = learn_threshold / 10;
    % genetic options
    gao.xover_prob      = 0.8;
    gao.mutation_prob   = 0.1;
    % make population size even
    if(mod(pop_size, 2) == 1) pop_size = pop_size + 1; end
    % init population
    pop = init(pop_size, lchrom, genes_margin);
    % history container
    hist = {};
    % invoke machine learning mode
    [pop, desc]...
        = lem_learn_mode       (...
            pop,...
            learn_probe,...
            learn_threshold,...
            @elit_ratio_func,...
            @save_hist);
    disp('Machine Learning mode has done...');
    % invoke darwinian evolutionary mode
    pop = lem_darwinian_mode   (...
            pop,...
            dar_probe,...
            dar_threshold,...
            gao,...
            genes_margin,...
            @(~) elit_ratio,...
            @save_hist);
    disp('Darwinian mode has done...');
    % resort and print final outputs
    pop = sortrows(pop, size(pop, 2));
    fprintf('Best found solution has fitness of: %f\n', pop(1, end));
    % return the history
    h = hist;
end

function ratio = elit_ratio_func(genNO)
% ELIT_RATIO_FUNC
% @brief        Calculates the elit selection ratio dynamically
% @param genNO  The current genetation# which elits are selecting from
% @return       The elit selection ratio 
    global elit_ratio
    ratio = elit_ratio / sqrt(genNO);
end

function save_hist(pop)
% SAVE_HIST
% @brief        Saves history of populations
% @param pop   The population to be saved
    global hist;
    hist{end+1} = pop;
end