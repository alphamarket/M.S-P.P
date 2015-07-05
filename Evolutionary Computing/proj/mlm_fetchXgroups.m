function [HP, LP] = mlm_fetchXgroups(pop, ratio)
% MLM_FETCHXGROUPS
% @brief            Fetched two H-groups and L-group from a popualtion
% @param  pop      	The population to fetch x-group from
% @param  ratio     The ratio of elits to be selected
% @return           The fetched H-group and L-group
    if ratio > 0.5, ratio = 0.5; end
    HP = pop(1:ceil(size(pop) *      ratio )    , :); 
    LP = pop(  ceil(size(pop) * (1 - ratio)):end, :);
end