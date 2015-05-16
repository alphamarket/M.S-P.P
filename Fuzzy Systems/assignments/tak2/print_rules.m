function print_rules(rules)
    for r = 1:size(rules, 1)
        if(rules(r, 4) == 0)
            fprintf('[%i] No rule generated for region [A%i, A%i]\n', rules(r, 1), rules(r, 2), rules(r, 3));
            continue
        end
        fprintf('[%i] if x1 is A%i and x2 is A%i then x belongs to G%i with %.2f confidence.\n', rules(r, 1), rules(r, 2), rules(r, 3), rules(r, 4), rules(r, 5));
    end
end