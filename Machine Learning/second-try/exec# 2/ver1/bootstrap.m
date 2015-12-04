clc, clear
% section 1
[x, y] = xlsread('data/tak2train.xlsx');
label = y(:, end);
label(1) = [];
marks = x(:, end);
x(:, end) = [];
x(:, 1) = [];
[x_test, ~] = xlsread('data/tak2test.xlsx');

% section 2
w_sec2 = pinv(x' * x) * x' * marks
test_grads_sec2 = x_test * w_sec2

% section 3
eta = 2e-5;
w_sec3 = rand(size(x_test, 2), 1);
mse = [];
for iter=1:100
	diffs = (marks - x * w_sec3);
	w_sec3 = w_sec3 + eta * (diffs' * x)';
	mse(end + 1) = sum(abs(diffs) .^ 2) / numel(x); 
end
plot(mse, '-.b*')
w_sec3
test_grads_sec3 = x_test * w_sec3

% section 4
marks(marks < 12) = 0;
marks(marks >= 12) = 1;

w_sec4 = pinv(x' * x) * x' * marks
test_grads_sec4 = x_test * w_sec4