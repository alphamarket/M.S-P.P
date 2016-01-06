function [net, info] = script
    clc, clear
    load_data
    imdb = getImdb(trainset_data, trainset_label, testset_data, testset_label);
    
    opts.train.batchSize = 1;
    opts.train.numEpochs = 30;
    opts.train.learningRate = 0.001 ;

    opts.expDir = fullfile('data','baseline');
    opts.networkType = 'simplenn' ;
    
    f=1/100 ;
    net.layers = {} ;
    net.layers{end+1} = struct('type', 'conv', ...
                               'weights', {{f*randn(10,10,1,20, 'single'), zeros(1, 20, 'single')}}, ...
                               'stride', 1, ...
                               'pad', 0) ;
    net.layers{end+1} = struct('type', 'pool', ...
                               'method', 'max', ...
                               'pool', [5 5], ...
                               'stride', 5, ...
                               'pad', 0) ;
    net.layers{end+1} = struct('type', 'conv', ...
                               'weights', {{f*randn(10,10,20,6, 'single'), zeros(1,6,'single')}}, ...
                               'stride', 1, ...
                               'pad', 0) ;
    net.layers{end+1} = struct('type', 'softmax') ;
    
    % Meta parameters
    net.meta.inputSize = [256 256 1] ;
    net.meta.trainOpts.learningRate = 0.001 ;
    net.meta.trainOpts.numEpochs = 20 ;
    net.meta.trainOpts.batchSize = 10 ;

    % Fill in defaul values
    net = vl_simplenn_tidy(net) ;

    [net, info] = cnn_train(net, imdb, getBatch, ...
      opts.train, ...
      'val', find(imdb.images.set == 3)) ;

function fn = getBatch
    fn = @(x,y) getSimpleNNBatch(x,y);
    
function [images, labels] = getSimpleNNBatch(imdb, batch)
    images = imdb.images.data(:,:,:,batch);
    labels = imdb.images.labels(:,batch);
