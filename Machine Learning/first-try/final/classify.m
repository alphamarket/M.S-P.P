function [class, prob, hist] = classify(X, thetas)
    class = -1;
    prob = -inf;
    X = [1 X];
    hist = [];
    for c=1:size(thetas, 1)
        p = sigmoid(X*thetas{c});
        if(prob < p)
            class = c;
            prob =  p;
        end
        hist(end + 1, :) = [c, p]; %#ok<AGROW>
    end
    hist = abs(sortrows(-hist, 2));
end

