function v = mf_trapezoid(x, i, k)
    assert(i <= k && i >= 1, 'value overflow');
    % eq. (4)
    aik = (i - 1) / (k-1);
    % eq. (5)
    bk  = 1 / (k-1);
    % eq. (6)
    v = max(min(2 - 2 * abs(x - aik) / bk, 1), 0);
end