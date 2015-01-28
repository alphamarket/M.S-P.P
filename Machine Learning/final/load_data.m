function [rtrv_data, raw_data] = load_data(addresses, labels, category_count)
    % retrieved data
    rtrv_data = cell(category_count, 1);
    raw_data = [];
    update_required = 0;
    for cat=1:size(addresses, 2)
        disp(['Processing:    ' cell2mat(addresses(cat))]);
        [~, name] = fileparts(cell2mat(addresses(cat)));
        cache_file = sprintf('data/cached/%s.cxr.mat', name);
        lbl = labels{cat};
        if exist(cache_file, 'file')
            disp('    Loading from cached data file...');
            rtrv_data{lbl} = load(cache_file, 'd', '-ascii');
            rtrv_data{lbl} = rtrv_data{lbl};
        else
            d = normalize_data(load_data_from_images(cell2mat(addresses(cat)), lbl));
            save(cache_file, 'd', '-ascii');
            raw_data = [raw_data; d]; %#ok<AGROW>
            rtrv_data{lbl} = d;
            update_required = 1;
        end
    end
    if update_required
        save('data/cached/raw.cxr.mat', 'raw_data', '-ascii');
    else
        raw_data = load('data/cached/raw.cxr.mat', '-ascii');
    end
end