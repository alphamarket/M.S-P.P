norm_trainset = normc(trainset(:, 1:end-1));
norm_testset  = normc(testset(:, 1:end-1));
dim_b3 = 1000;
fprintf('Reducing deminsions from %d -> %d....', size(trainset, 2) - 1, dim_b3);
eigvector = PCA(norm_trainset, dim_b3);
trainset_b3 = (norm_trainset * eigvector);
testset_b3  = (norm_testset * eigvector);
fprintf('[DONE]\nTraining MLP....');
x = trainset_b3'; t = trainset(:, end)';
% set hidden layer
cnet_b3 = 30;
% configure the network
net_b3 = configure(feedforwardnet(cnet_b3), x, t); 
net_b3.trainFcn = 'trainscg';
net_b3.trainParam.showWindow = 0;
% train the nn
net_b3 = train(net_b3, x, t);
fprintf('[DONE]\nTesting model....');
% predict the tests
pred_b3 = round(net_b3(testset_b3'));
pred_b3(pred_b3 > 2) = 2;
pred_b3(pred_b3 < 1) = 1;
fprintf('[DONE]\nConfusion Matrix:\n');
C = confusionmat(pred_b3, testset(:, end)');
disp(dataset({C, 'pred_1', 'pred_2'}, 'obsnames', {'real_1', 'real_2'}))

clear C x t norm_trainset norm_testset eigvector