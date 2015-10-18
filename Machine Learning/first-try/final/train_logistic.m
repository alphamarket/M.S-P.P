function results = train_logistic(cat_data, prem_data)
    % occupy memory for result outcomes
    results = cell(size(cat_data));
    % occupy memory for end indexes of each categories
    end_indexs = [0; zeros(size(cat_data))];
    % for each category
    for d=1:size(cat_data, 1)
        % compute the end index of that category
        end_indexs(d + 1) = end_indexs(d) + size(cat_data{d}, 1);
    end
    % for each categoory
    for d=1:size(cat_data, 1)
        fprintf('Traing for label# %i\n', d);
        % train with respect to the category
        % store the training results
        results{d} = train_for_label(prem_data, d, end_indexs);
    end
end

function theta = train_for_label(prem_data, label, end_indexes)
    % fail safe for valid data
    if label > size(prem_data, 1) && label < 0
        error('label overflow');
    end
    % mute the labels except current label
    prem_data = mute_labels(prem_data, label, end_indexes);
    % fetch features
    X = prem_data(:, 1:end-1);
    % fetch the labels
    y = prem_data(:, end);
    % fetch the sizes
    [m, n] = size(X);
    % add the '1' threshold
    X = [ones(m, 1) X];
    % intial theta - all zero
    initial_theta = zeros(n + 1, 1);
    % define a cost function
    costFunction = @costFunctionReg;
    lambda = 1;
    % compute and display initial cost and gradient
    fprintf('Cost at initial theta (zeros): %f\n', costFunction(initial_theta, X, y, lambda));
    % set options for fminunc
    options = optimset('GradObj', 'on', 'MaxIter', 10000, 'Display', 'iter');
    disp('Start training....');
    % run fminunc to obtain the optimal theta
    % this function will return theta and the cost 
    [theta, cost] = fminunc(@(t)(costFunction(t, X, y, lambda)), initial_theta, options);
    % print the final cost value returned by fminunc
    fprintf('Cost at theta found by fminunc: %f\n', cost);
end

function data = mute_labels(data, except_label, end_indexes)
   % find the current's label's end index
   i = end_indexes(except_label:except_label+1);
   % all labels go zero
   data(:,end) = 0;
   % make current's label one
   data(i(1)+1:i(2), end) = 1;
end