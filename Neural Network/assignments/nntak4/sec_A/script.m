clc, clear
images = cell(6, 1);
fprintf('Loading CNN data...');
% load the pre-trained CNN
net = load('../imagenet-vgg-f.mat');
fprintf('[DONE]\nLoading %d custom images...\n', size(images, 1));
figure
for i=1:size(images, 1)
    image_name = sprintf('images/%d.jpg', i);
    images{i} = imread(image_name);
    fprintf('Preprocessing the image %s\n', image_name);
    im_ = single(images{i}) ; % note: 0-255 range
    im_ = imresize(im_, net.meta.normalization.imageSize(1:2)) ;
    im_ = im_ - net.meta.normalization.averageImage;
    fprintf('Classifing the image %s\n', image_name);
    % run the CNN
    res = vl_simplenn(net, im_) ;
    % show the classification result
    scores = squeeze(gather(res(end).x)) ;
    [bestScore, best] = max(scores);
    subplot(2, 3, i); imagesc(images{i}) ;
    title(sprintf('%s (%d), score %.3f', net.meta.classes.description{best}, best, bestScore));
end