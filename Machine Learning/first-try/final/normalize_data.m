function data = normalize_data(data)
    % do the normalization
    y = [];
    normalize_factor = 10;
    for i =1:size(data, 1)
        x = [0]; %#ok<NBRAK>
        for j=1:normalize_factor:size(data(i, :), 2) - (normalize_factor + 1)
            x = [x mean(data(i, j:j+normalize_factor))]; %#ok<AGROW>
        end
        y = [y; [x data(i, end)]]; %#ok<AGROW>
    end
    data = y;
end