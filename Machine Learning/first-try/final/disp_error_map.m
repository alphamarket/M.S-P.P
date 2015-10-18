function disp_error_map(em)
    k = keys(em);
    v = values(em);
    for i=1:size(v, 2)
        x = v{i};
        fprintf('[%i] => { e: %f | c: %i }\n', k{i}, x(1), x(2));
    end
end