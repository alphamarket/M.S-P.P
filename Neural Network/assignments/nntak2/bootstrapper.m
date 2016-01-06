function [net, tr, images, targets] = bootstrapper(section, cnet)
    clc
    if(nargin < 2), cnet = 100; end
    disp('$ loading datasets');
    images = load('data/training/images.mat');
    train_labels = load('data/training/labels.mat');
    test_images = load('data/test/images.mat');
    test_labels = load('data/test/labels.mat');
    images = [images.i, test_images.images];
    train_labels = [train_labels.j; test_labels.labels];
    targets = zeros(10, size(images, 1));
    for i=1:size(images, 2), targets(train_labels(i) + 1, i) = 1; end
    assert(size(images, 2) == size(targets, 2));
    disp('$ configuring the network');
    net = configure(feedforwardnet(cnet), images, targets);
    switch(section)
        case 1.1
            disp('$ batch training');
            net.trainFcn = 'traingd';
        case 1.2
            disp('$ online training');
            net.trainFcn = 'trainc';
        case 1.3
            disp('$ gradient decent with momentum');
            net.trainFcn = 'traingdm';
        case 3.11
        	net.trainFcn = 'trainscg';
        case 3.12
        	net.trainFcn = 'trainrp';
        case 3.13
    		% makes mem. overflow!!
        	net.trainFcn = 'trainbr';
    	case 3.14
    		% makes mem. overflow!!
        	net.trainFcn = 'trainlm';
        case 3.21
            net.trainFcn = 'traingd';
        	net.performFcn = 'crossentropy';
        case 3.22
            net.trainFcn = 'trainscg';
        	net.performFcn = 'crossentropy';
        case 3.23
            net.trainFcn = 'trainrp';
        	net.performFcn = 'crossentropy';
        otherwise
            error(strcat('Undefined section#', section))
    end
    net.divideFcn='divideind';
    [net.divideParam.trainInd,net.divideParam.valInd,net.divideParam.testInd] = divideind(size(images, 2), 1:5e+4, (5e+4) + 1:6e+4, (6e+4)+1:7e+4);
    net.trainParam.show = 1;
    net.trainParam.goal = 0.01;
    net.trainParam.epochs = 30;
    sprintf('Images.size = '), size(images), sprintf('Labels.size = '), size(targets)
    disp('$ starting to train');
    [net, tr] = train(net, images, targets);
    disp('$ training done');
end