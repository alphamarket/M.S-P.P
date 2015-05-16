function rules = gen_rules(l_max, x1, x2, func, mf_func, mf_opr, apply_cf)
    assert(l_max > 1, 'variable undetflow');
    if nargin < 7
        apply_cf = true;
    end
    if nargin < 6
        mf_opr = @(mv1, mv2) mv1 * mv2;
    end
    if nargin < 5
        mf_func = @mf_triangular;
    end
    rules = [];
    if(size(x1, 1) ~= 1) x1 = x1'; end
    if(size(x2, 1) ~= 1) x2 = x2'; end
    % setup the training set
    train_set = [x1', x2', value2class(func(x1, x2))'];
    % the implementation part of eq. (8, 9, 13, 14)
    for k=2:l_max
        for i = 1:k
            for j = 1:k
                bg1i = train_set(train_set(:, 3) == 1, :);
                bg2i = train_set(train_set(:, 3) == 2, :);
                bg1 = 0; bg2 = 0;
                for c=1:size(bg1i, 1)
                    bg1 = bg1 + mf_opr(mf_func(bg1i(c, 1), i, k), mf_func(bg1i(c, 2), j, k));
                end
                for c=1:size(bg2i, 1)
                    bg2 = bg2 + mf_opr(mf_func(bg2i(c, 1), i, k), mf_func(bg2i(c, 2), j, k));
                end
                cf = 1;
                if apply_cf
                    % eq. (12)
                    cf = abs(bg1 - bg2) / (bg1 + bg2);
                end
                % the numerical rule structure
                % format : 
                %   [K] [region{x1, x2}] [class] [confidence]
                rules(end+1, :) = [k i j 0 cf]; %#ok<AGROW>
                % eq. (10, 11)
                if(bg1 < bg2)
                    rules(end, 4) = 2;
                elseif(bg1 > bg2)
                    rules(end, 4) = 1;
                end
            end
        end
    end
end