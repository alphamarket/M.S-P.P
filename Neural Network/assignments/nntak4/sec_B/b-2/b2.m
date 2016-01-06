fprintf('Loading CNN data...');
% load the pre-trained CNN
net_gg = load('../imagenet-vgg-f.mat');
fprintf('[DONE]\nExtracting features....');
train_features = zeros(size(trainset, 1), 1000);
test_features = zeros(size(testset, 1), 1000);
% extract trainset's features
for i=1:size(trainset, 1)
    % preprocess for network configuration
    image = single(reshape(trainset(i, 1:end-1), net_gg.meta.normalization.imageSize(1:2)));
    % fake an RGB
    image = cat(3, image, image, image);
    % feed the network
    res = vl_simplenn(net_gg, image);
    % fetch the features
    train_features(i, :) = res(end-1).x;
end
% extract testset's features
for i=1:size(testset, 1)
    % preprocess for network configuration
    image = single(reshape(testset(i, 1:end-1), net_gg.meta.normalization.imageSize(1:2)));
    % fake an RGB
    image = cat(3, image, image, image);
    % feed the network
    res = vl_simplenn(net_gg, image);
    % fetch the features
    test_features(i, :) = res(end-1).x;
end
fprintf('[DONE]\nTraining MLP....');
x = train_features'; t = trainset(:, end)';
% set hidden layer
cnet_b2 = 30;
% configure the network
net_b2 = configure(feedforwardnet(cnet_b2), x, t); 
net_b2.trainFcn = 'trainscg';
net_b2.trainParam.showWindow = 0;
% train the nn
net_b2 = train(net_b2, x, t);
fprintf('[DONE]\nTesting model....');
% predict the tests
pred_b2 = round(net_b2(test_features'));
pred_b2(pred_b2 > 2) = 2;
pred_b2(pred_b2 < 1) = 1;
fprintf('[DONE]\nConfusion Matrix:\n');
C = confusionmat(pred_b2, testset(:, end)');
disp(dataset({C, 'pred_1', 'pred_2'}, 'obsnames', {'real_1', 'real_2'}))

clear C x t net_gg res image i
