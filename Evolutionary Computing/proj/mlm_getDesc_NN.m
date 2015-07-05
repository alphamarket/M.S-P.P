function net = mlm_getDesc_NN(pop, HP, LP) %#ok<INUSL>
% MLM_GETDESC_NN
% @brief        Get a neural network description
% @param  pop  	The curret population which H-group and L-group has been selected from
% @param  HP    The H-group population
% @param  LP    The L-group population
% @return       The fited neural network to x-groups

    % prepare groups to train the nn
    HP(1:size(HP), end) = 0; LP(1:size(LP), end) = 1;
    % the nn's input & target pop
    x = [HP(:,1:end-1); LP(:,1:end-1)]'; t = [HP(:,end); LP(:,end)]';
    % try to fetch best nn layar configuration
    cnet = [3]; %#ok<NBRAK>
%     qlay = floor(sqrt(sqrt(size(pop, 2)))); 
%     while(qlay > 1), cnet = [cnet cnet]; qlay = floor(sqrt(qlay)); end %#ok<AGROW>
    net = configure(feedforwardnet(cnet), t, x); net.trainParam.showWindow = false;
    % train the nn
    net = train(net, t, x);
end