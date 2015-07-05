function mlm_stat(pop, genNO, HP, LP) %#ok<INUSL>
% MLM_STAT
% @brief        Prints a statistical report for a population
% @param  pop  	The population
% @param  genNO The current generation#
% @param  HP    The H-group population
% @param  LP    The L-group population
    fprintf('    %10i%-10s%.5e%-6s%.5e%-7s%.5e%-7s%.5e\n', ...
        genNO, ' ', LP(end,end), ' ', LP(1,end), ' ', HP(end,end), ' ', HP(1,end));
end