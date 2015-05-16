function z = classifier(rules, xp, mf_func, mf_opr, apply_cf)
    if nargin < 5
        apply_cf = true;
    end
    if nargin < 4
        mf_opr = @(mv1, mv2) mv1 .* mv2;
    end
    if nargin < 3
        mf_func = @mf_triangular;
    end
    % fetch associated rules for each class
    g1 = rules(rules(:, 4) == 1, :);
    g2 = rules(rules(:, 4) == 2, :);
    ag1 = 0; ag2 = 0;
    % eq. (16, 17, 18, 19)
    for i=1:size(g1, 1)
        j = g1(i, :);
        x = mf_opr(mf_func(xp(:, 1), j(2), j(1)), mf_func(xp(:, 2), j(3), j(1)));
        if apply_cf
            x = x .* j(end);
        end
        ag1 = max(x, ag1);
    end
    for i=1:size(g2, 1)
        j = g2(i, :);
        x = mf_opr(mf_func(xp(:, 1), j(2), j(1)), mf_func(xp(:, 2), j(3), j(1)));
        if apply_cf
            x = x .* j(end);
        end
        ag2 = max(x, ag2);
    end
    assert(isequal(size(ag1), size(ag2)), 'size mis-matched');
    % prepare classification output results
    % format:
    %   [class] [confidence]
    z = zeros(size(ag1, 1), 2);
    for i = 1:size(z, 1)
        w = 0;
        if(ag1(i) > ag2(i)), w = 1;
        elseif(ag1(i) < ag2(i)), w = 2;
        end
        z(i, :) = [w abs(ag1(i) - ag2(i))];
    end
end