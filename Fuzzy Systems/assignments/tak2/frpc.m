function frpc(l_max)
    if nargin < 1
        l_max = 10;
    end
    func = @(x1, x2) -0.25*sin(2*pi*x1) + x2 - 0.5;
    patterns_size = [40, 100];
    % generate test patterns
    [t_x1, t_x2] = gen_data(100);
    xp = [t_x1', t_x2', value2class(func(t_x1, t_x2))'];
    % foreach pattern size
    for ps = 1:size(patterns_size, 2)
        fprintf('The number of patterns = %i\n\n', patterns_size(ps));
        %
        % rule generation part
        %
        % generate training data
        [x1, x2] = gen_data(patterns_size(ps));
        for iter=1:6
            % get current iteration's training/testing options
            [mf_func, mf_opr_train, mf_opr_test, apply_cf] = get_opt(iter);
            for l=2:l_max
                % generate rules
                rules = gen_rules(l, x1, x2, func, mf_func, mf_opr_train, apply_cf);
                [true_class, ~] = test_rule(rules, [x1', x2', value2class(func(x1, x2))'], mf_func, mf_opr_test, apply_cf);
                if(true_class == 1), break; end
            end
            % print_rules(rules);
            %
            % testing part
            % 
            [true_class, unclass] = test_rule(rules, xp, mf_func, mf_opr_test, apply_cf);
            fprintf('Correct = %.2f\t\tUnclass = %.2f\n', true_class , unclass);
        end
        fprintf('\n');
    end
end

function [true_class, unclass] = test_rule(rules, xp, mf_func, mf_opr_test, apply_cf)
    % try to classify the test patterns
    z = classifier(rules, xp(:, 1:2), mf_func, mf_opr_test, apply_cf);
    % print_classification_res(xp(:, 1:2), z);
    %
    % evaluating part
    % 
    assert(isequal(size(z, 1), size(xp, 1)), 'size mis-matched');
    true_class = 0; unclass = 0;
    for i=1:size(z, 1)
        if z(i, 1) == 0, unclass = unclass + 1;
        elseif z(i, 1) == xp(i, 3), true_class = true_class + 1; end
    end
    true_class = true_class / size(z, 1);
    unclass = unclass / size(z, 1);
end

function [mf_func, mf_opr_train, mf_opr_test, apply_cf] = get_opt(iter)
    mf_func = @mf_triangular;
    mf_opr_train = @(mv1, mv2) mv1 * mv2;
    mf_opr_test = @(mv1, mv2) mv1 .* mv2;
    apply_cf = true;
    switch iter
        case 1
            fprintf('Triangular\tProduct\t\tWith CF\t\t\t');
        case 2
            mf_func = @mf_trapezoid;
            fprintf('Trapezoid\tProduct\t\tWith CF\t\t\t');
        case 3
            apply_cf = false;
            fprintf('Triangular\tProduct\t\tWithout CF\t\t');
        case 4
            mf_func = @mf_trapezoid;
            apply_cf = false;
            fprintf('Trapezoid\tProduct\t\tWithout CF\t\t');
        case 5
            mf_opr_train = @(mv1, mv2) min(mv1, mv2);
            mf_opr_test = @(mv1, mv2) min(mv1, mv2);
            apply_cf = true;
            fprintf('Triangular\tMin\t\tWith CF\t\t\t');
        case 6
            mf_opr_train = @(mv1, mv2) min(mv1, mv2);
            mf_opr_test = @(mv1, mv2) min(mv1, mv2);
            apply_cf = false;
            fprintf('Triangular\tMin\t\tWithout CF\t\t');
        otherwise
            error('iter overflow');
    end
end