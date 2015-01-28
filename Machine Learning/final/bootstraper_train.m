clc, clear, close;
% map categories to there numerical values
adds = {'data/mri/body', 'data/mri/brain', 'data/mri/knee', 'data/mri/shoulder', 'data/mri/spine'};
lbls = {1, 2, 3, 4, 5};
% fail-safe
if(size(adds) ~= size(lbls))
    error('Size mis-matched');
end
% load data
[rtrv_data, raw_data] = load_data(adds, lbls, max(cell2mat(lbls)));
% add Logistic Regression utilities to current path
addpath('lr', 'jsonlab');
% load options
opt = loadjson('opt.json');
% output the mixture degree value
fprintf('Mixing features with %i degree.\n', opt.permutation_degree);
% anonymous function for training
train = @(train_data) train_logistic(rtrv_data, train_data);
% define result variable
train_results = []; %#ok<NASGU>
% train by demand
switch(opt.select_train_method)
    case 'test_train'
        if(opt.load_from_saved_thetas)
            l = load(input('Please enter the source file''s path: ', 's'));
            train_results = l.train_results;
            opt.permutation_degree = l.permutation_degree;
        end
        % genetate randomized index vector
        indexes = randperm(size(raw_data, 1));
        % fetch train indexes
        train_indexes = indexes(1:floor(opt.test_train_method.train * size(raw_data, 1)));
        % fetch test indexes
        test_indexes  = indexes(size(train_indexes, 2) + 1:end);
        % generate permuted data from original data
        perm_data = [permute_map_features(raw_data(:, 1:end-1), opt.permutation_degree) raw_data(:, end)];
        % fetch training actual data
        train_data = perm_data(train_indexes, :);
        % fetch test actual data
        test_data  = perm_data(test_indexes,  :);
        % train with training data
        if(~opt.load_from_saved_thetas), train_results = train(train_data); else test_data = perm_data; end
        [general_error, errors_map] = test_logistic(train_results, test_data);
        fprintf('General error rate: %2.2f\nClass error rate: \n', general_error * 100);
        disp_error_map(errors_map);
        if(opt.load_from_saved_thetas)
            return
        end
    otherwise
        error('Invalid training method: %s', opt.select_train_method);
end
% ask for save results
while(1)
    switch(lower(input('Do you want to save results?[Y/n] ', 's')))
        case 'y'
            file = input('Please enter the destination file''s name: ', 's');
            permutation_degree = opt.permutation_degree;
            save(file, 'train_results', 'permutation_degree');
            disp('Saved.');
            break
        case 'n'
            return
    end
end