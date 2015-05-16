function v = mf_triangular(x, i, k)
    assert(i <= k && i >= 1, 'value overflow');
    % eq. (4)
    aik = (i - 1) / (k-1);
    % eq. (5)
    bk  = 1 / (k-1);
    % eq. (3)
    v = max(1 - abs(x - aik) / bk, 0);
end