function features = permute_map_features(features, degree)
    if(degree > 0)
        s=size(features, 2);
        for i=1:s
            for j=1:s
                if i==j 
                    continue 
                end
                features = [features mapFeature(features(:, i), features(:, j), degree)]; %#ok<AGROW>
            end
        end
    end
    features = log(features + 1);
end