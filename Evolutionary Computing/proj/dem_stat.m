function dem_stat(pop, genNO)
% DEM_STAT
% @brief        Prints a statistical report for a population
% @param  pop  	The population
% @param  genNO The current generation#
    a = pop(:,end);
    frmt_str = '    %10i%-10s%.5e%-10s%.5e%-7s%.5f\n';
    if isinf(max(a)), frmt_str = '    %10i%-10s%.5e%-18s%.5e%-7s%.5f\n'; end
    fprintf(frmt_str, genNO, ' ', max(a), ' ', mean(a(~isinf(a))), ' ', min(a));
end