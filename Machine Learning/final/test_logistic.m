function [error, error_map] = test_logistic(thetas, test_data)
    error = 0;
    dispstat('','init');
    error_map = containers.Map('KeyType','double', 'ValueType','any');
    misclass_map = containers.Map('KeyType','double', 'ValueType','any');
    pref_map = containers.Map('KeyType','double', 'ValueType','any');
    for t=1:size(test_data, 1)
        dispstat(sprintf('%%%i of testing is done....', floor(t*100 / size(test_data, 1))));
        i = test_data(t, 1:end-1);
        ac = test_data(t, end);
        if(~error_map.isKey(ac))
            error_map(ac) = [0 0];
            misclass_map(ac) = containers.Map('KeyType','double', 'ValueType','any');
            pref_map(ac) = [0 0]; % [ { TPR = TP/P } { FPR = FP/N }]
        end
        em = error_map(ac);
        error_map(ac) = [em(1) em(2) + 1];
        [class, prob, hist] = classify(i, thetas);
        switch true
            case class == ac        
                if(ac==5)
                    ;
                end
                pm = pref_map(ac);
                pref_map(ac) = [pm(1) + 1 pm(2)]; % inc TP by one
            case class ~= ac
                if(~pref_map.isKey(class))
                    pref_map(class) = [0 0]; % [ { TPR = TP/P } { FPR = FP/N }]
                end
                pm = pref_map(class);
                pref_map(class) = [pm(1) pm(2) + 1]; % inc FP by one
                error = error + 1;
                error_map(ac) = [(em(1) * size(test_data, 1) + 1) / size(test_data, 1) em(2) + 1];
                mm = misclass_map(ac);
                if(~mm.isKey(hist(2)))
                    mm(hist(2)) = 0;
                end
                mm(hist(2)) = mm(hist(2)) + 1;
                misclass_map(ac) = mm;
        end
    end
    for i=1:5
        if(~isKey(error_map, i)), continue; end
        em = error_map(i);
        pm = pref_map(i);
        N = 0;
        for j=1:5
            if(~isKey(error_map, j)), continue; end
            if i == j
                continue
            end
            tem = error_map(j);
            N = N + tem(2);
        end
        pref_map(i) = [pm em(2) N]; % [ TPR FPR P N ]
        pm = pref_map(i);
        fprintf('For classifier# %i: { TPR : %.5f | FPR : %.5f | ROC : %.4f }\n', i, pm(1) / pm(3), pm(2) / pm(4), (pm(1) / pm(3))/(pm(2) / pm(4)));  
    end
    error = error / size(test_data, 1);
    disp_misclass_map(misclass_map)
end

function disp_misclass_map(mm)
    k = keys(mm);
    v = values(mm);
    for i=1:size(v, 2)
        fprintf('----------- %i -----------\n', k{i});
        kk = keys(v{i});
        vv = values(v{i});
        x = [];
        for j=1:size(kk, 2)
            x = [x; kk{j} vv{j};]; %#ok<AGROW>
        end
        -sortrows(-x, 2)
    end
end