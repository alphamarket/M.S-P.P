fprintf('Loading images....');
if(exist('data.cache', 'file'))
    in = importdata('data.cache');
    data = in.data;
    label = in.label;
else
    data = [];
    label = [];
    meta = importdata('images/labels.dat');
    for j=1:length(meta)
        C = textscan(meta{j}, '%f %f %f %f %f %f %s');
        file = C{end}{1};
        if(~exist(sprintf('images/%s.tiff', file), 'file')), continue; end
        % read and convert to gray
        image = imread(sprintf('images/%s.tiff', file));
        % vectorize and store the image
        data(end+1, :) = reshape(image, [1 size(image, 1) * size(image, 2)]);
        % store the label
        label(end+1, :) = [C{1:end-1}];
    end
    save('data.cache', 'data', 'label');
end
fprintf('[DONE]\n');
shuffle = randperm(size(data, 1));
trainset_indices = shuffle(1:floor(0.8*length(shuffle)));
testset_indices  = shuffle(ceil(0.8*length(shuffle)):end);
trainset_data = data(trainset_indices, :);
testset_data  = data(testset_indices, :);
trainset_label = label(trainset_indices, :);
testset_label = label(testset_indices, :);

clear in shuffle trainset_indices testset_indices