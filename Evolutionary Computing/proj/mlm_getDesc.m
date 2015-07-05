function desc = mlm_getDesc(pop, HP, LP, clss)
% MLM_GETDESC
% @brief        Get a description based on passed H-group and L-group
% @param  pop  	The curret population which H-group and L-group has been selected from
% @param  HP    The H-group population
% @param  LP    The L-group population
% @param  clss  The class of description generator `nn => neural network`
% @return       The generated description
    switch clss
        case 'nn'
            desc = mlm_desc('nn', mlm_getDesc_NN(pop, HP, LP));
        otherwise
            error('description class `%s` is not defined!', clss); 
    end
end