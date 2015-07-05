function iplot(hist, k, w)
    close all;
    if nargin < 2, k = 1; end
    if nargin < 3, w = size(hist, 2); end
    z = ceil(sqrt(w - k + 1));
    l = [];
    for i=k:w
        y = hist{i};
        l = [l subplot(z, z, i - k + 1)]; %#ok<AGROW>
        plot(y(:, 1), y(:, 2), 'b.'), title(i);
    end
    linkaxes(l, 'xy');
end