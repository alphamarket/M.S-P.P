function imdb = getImdb(imDir)
    counter         = 0;
    imSize          = [28 28 1];
    imSet           = imList(imDir);
    imSeeds         = randperm(length(imSet));
    trainIndices    = 1:floor(length(imSeeds)*.8);
    testIndices     = floor(length(imSeeds)*.8)+1:length(imSeeds);
    trainSet        = imSet(imSeeds(trainIndices));
    validationSet   = imSet(imSeeds(testIndices));
    labels          = zeros(length(imSet), 1);
    set             = zeros(length(imSet), 1);
    Images          = zeros([imSize, length(imSet)], 'single');
    labelstr        = { 'ANG', 'DEP', 'HAP', 'SAD', 'SUP', 'NEU' };

    [Images, labels, set, counter]...
                    = readSet(trainSet, Images, labels, labelstr, set, 1, counter, imSize);
    [Images, labels, set, ~]...
                    = readSet(validationSet, Images, labels, labelstr, set, 3, counter, imSize);
    
    data            = Images;
    dataMean        = mean(data(:, : , :, set == 1), 4);
    data            = bsxfun(@minus, data, dataMean);
    classes         = cell(1, length(trainSet));
    
    for i = 1 : length(trainSet)
        classes{i}  = labelstr{labels(i)};
    end
    
    imdb.meta.sets  = {'train', 'val', 'test'};
    r               = randperm(length(set), length(set));

    imdb.images.set = set(r);
    imdb.images.data...
                    = data(:, :, :, r);
    imdb.images.data_mean...
                    = dataMean;
    imdb.images.labels...
                    = labels(r);
    
    imdb.meta.classes...
                    = classes;

function out = imList(baseDir)
    out = {};
    stack = push_cell({}, baseDir);
    while(~isempty(stack))
        [stack, cwd] = pop_cell(stack);
        files = dir(cwd);
        for index = 1:length(files)
            file = files(index);
            if(strcmp(file.name, '.') || strcmp(file.name, '..')), continue; end
            full_path = strcat(cwd, '/', file.name);
            if(file.isdir), stack = push_cell(stack, full_path);
            else
                [~, ~, ext] = fileparts(full_path);
                if(strmatch(ext, {'.jpg', ',jpeg', '.png', '.tiff'}, 'exact') > 0) %#ok<MATCH3>
                    out = push_cell(out, full_path); 
                end
            end
        end
    end

function cell = push_cell(cell, item)
    if(isempty(cell)), cell = {item}; return; end
    cell = [cell; item];

function [cell, item] = pop_cell(cell)
    if(isempty(cell)), error('cannot pop! cell is empty already!'); end
    item = cell{1, :};
    cell(1, :) = [];

function [Images, labels, set, counter] = readSet(imset, Images, labels, labelstr, set, set_value, counter, imSize)
    for jj = 1:size(imset, 1)
      counter = counter +1; 
      Images(:,:,:, counter) = imresize(single(imread(imset{jj})),imSize(1:2));
      for l=1:length(labelstr)
          if(strfind(imset{jj}, strcat('/', labelstr{l}, '/')))
              labels(counter) = l;
          end
      end
      set(counter) = set_value;
    end