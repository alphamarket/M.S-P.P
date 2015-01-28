clc;
index = floor(rand * size(raw_data, 1));
fprintf('Testing subject# %i with known class# %i ', index, raw_data(index, end));
test_subject = raw_data(index, 1:end-1);
[class, prob, hist] = classify(test_subject, results);
fprintf(' : Classified as# %i with probability of %f\n', class, prob);
Possibilities = abs(sortrows(-hist, 2)) %#ok<NOPTS>