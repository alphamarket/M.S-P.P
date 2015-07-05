function tsp_plot_output() 
	close all;
	disp('loading data....');
	x=load('out.dat'); 
	y=load('log.out');
	dist = 0;
	for k=1:size(x,1)-1
	    dist = dist + norm(x(k+1,:) - x(k,:),2);
	end
	dist = dist + norm(x(end,:) - x(1,:),2);
	fprintf('COST IS: %.2f\n', dist);
	disp('ploting data....');
	x(end+1, :) = x(1, :);
	xx = x(:,1); xy = x(1:size(xx, 1),2);
	zxx = x(:,1); zxy = x(1:size(zxx, 1),2);
	plot(xx, xy), title('Found TSP solution');
	hold
	plot(zxx, zxy, 'r.');
	figure
	plot(y(:,2:4)), legend('max. fitness', 'avg. fitness', 'min. fitness'), title('The evolution history')
end