function imdb = getImdb(trainset_data, trainset_label, testset_data, testset_label)
x1 = reshape(trainset_data', 256, 256, size(trainset_data, 1));
x2 = reshape(testset_data', 256, 256, size(testset_data, 1));
y1 = trainset_label';
y2 = testset_label';

set = [ones(1,size(y1, 2)) 3*ones(1,size(y2, 2))];
data = single(reshape(cat(3, x1, x2),256,256,1,[]));
dataMean = mean(data(:,:,:,set == 1), 4);
data = bsxfun(@minus, data, dataMean) ;

imdb.images.data = data;
imdb.images.data_mean = dataMean;
imdb.images.labels = cat(2, y1, y2);
imdb.images.set = set;
imdb.meta.sets = {'train', 'val', 'test'};
imdb.meta.classes = arrayfun(@(x)sprintf('%d',x),1:6,'uniformoutput',false);