function [d, fig1, fig2] = plot_details(file)
	d = load(file);
	close all, fig1 = figure; hold, grid on
	plot(d(:, 4), '-.g.')
	plot(d(:, 5), '-.r.')
	plot(d(:, 6), '-.b.')
	legend('Trainset', 'Evalset', 'Testset', 'location', 'southeast');
	xlabel('Epoch'), ylabel('Accuracy'), title('Accuracy over Train/Eval/Test sets')
	fig2 = figure; hold, grid on
	plot(d(:, 1), '-.g.')
	plot(d(:, 2), '-.r.')
	plot(d(:, 3), '-.b.')
	legend('Trainset', 'Evalset', 'Testset', 'location', 'northeast');
	xlabel('Epoch'), ylabel('MSE'), title('MSE over Train/Eval/Test sets')
	saveas(fig1, strcat(file, '_acc.png'));
	saveas(fig2, strcat(file, '_mse.png'));
end
