function [elits, pop] = dem_elit_selection(pop, elit_ratio)
% DEM_ELIT_SELECTION
% @brief                Selects elits
% @param  pop          	The population to select elits from
% @param  elit_ratio    The ratio of elits to be selected
% @return               The elits and sorted population
    pop = sortrows(pop, size(pop, 2));
    [elits, ~] = mlm_fetchXgroups(pop, elit_ratio);
end