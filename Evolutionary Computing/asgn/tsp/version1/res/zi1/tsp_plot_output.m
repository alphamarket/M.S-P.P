% function tsp_plot_output() 
% 	close all;
% 	disp('loading data....');
% 	x=load('out.res'); 
% 	% y=load('log.out');
% 	dist = 0;
% 	for i=1:2:size(x, 1)-1
% 		dist = dist + norm(x(i,:) - x(i+1,:), 2);
% 	end
% 	dist = dist + norm(x(1,:) - x(end,:), 2);
% 	fprintf('%.2f\n', dist);
% 	% return;
% 	disp('ploting data....');
% 	xx = x(:,1); xy = x(1:size(xx, 1),2);
% 	zxx = x(:,1); zxy = x(1:size(zxx),2);
% 	% plot(xx, xy, '.'), title('Cities location');
% 	% figure
% 	plot(xx, xy), title('Found TSP solution');
% 	hold
% 	plot(zxx, zxy, 'r.');
% 	% figure
% 	% plot(y(:,2:4)), legend('max. fitness', 'avg. fitness', 'min. fitness'), title('The evelution history')
% end
function tsp_plot_output() 
	close all;
	disp('loading data....');
	x=load('out.res'); 
	dist = 0;
	for i=1:2:size(x, 1)-1
		dist = dist + norm(x(i,:) - x(i+1,:), 2);
	end
	dist = dist + norm(x(1,:) - x(end,:), 2);
	fprintf('COST IS: %.2f\n', dist);
	disp('ploting data....');
	xx = x(:,1); xy = x(1:size(xx, 1),2);
	zxx = x(:,1); zxy = x(1:size(zxx),2);
	plot(xx, xy), title('Found TSP solution');
	hold
	plot(zxx, zxy, 'r.');
end