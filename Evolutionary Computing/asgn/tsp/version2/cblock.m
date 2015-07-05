function cluster()
clc, close all;
x=load('cblock.out');
cmin=min(x(:,:,:));
cmax=max(x(:,:,:));
% plot(x(:,1), x(:,2), 'rx');
figure;
hold;
fprintf('Clusteres are in [%i, %i]\n', cmin(3), cmax(3));
colorlist=hsv(cmax(3) + 1);
for i=cmin(3):cmax(3)
    i;
    A = x(find(x(:,3) == i), 1:2);
    plot(A(:,1), A(:,2), '.', 'Color', colorlist(i,:))
end
end
