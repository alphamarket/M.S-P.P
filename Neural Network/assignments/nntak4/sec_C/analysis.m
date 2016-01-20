close all, clear, clc
load('data/exec/5in1/data.mat');
val_top1 = [];
val_top5 = [];
trn_top1 = [];
trn_top5 = [];
trn_obj  = [];
val_obj  = [];
for i=1:length(x)
    info = x{i}.info;
    val_top1 = [val_top1; info.val.error(1, :)];
    val_top5 = [val_top5; info.val.error(2, :)];
    trn_top1 = [trn_top1; info.train.error(1, :)];
    trn_top5 = [trn_top5; info.train.error(2, :)];
    val_obj  = [val_obj; info.val.objective]; 
    trn_obj  = [trn_obj; info.train.objective]; 
end
val_top1_mean = mean(val_top1);
val_top5_mean = mean(val_top5);
trn_top1_mean = mean(trn_top1);
trn_top5_mean = mean(trn_top5);
val_obj_mean  = mean(val_obj); 
trn_obj_mean  = mean(trn_obj);

figure, hold;
plot(trn_top1_mean, '.-', 'linewidth', 2);
plot(trn_top5_mean, '.-', 'linewidth', 2, 'color',[0 0.5 0]);
plot(val_top1_mean, '.--');
plot(val_top5_mean, '.--', 'color',[0 0.5 0]);
legend('train top-1 e', 'train top-5 e', 'val top-1 e', 'val top-5 e');
xlabel('training epoch'); ylabel('avg. error');
title('error');
grid on;

figure, hold;
plot(trn_obj_mean, '.-', 'linewidth', 2);
plot(val_obj_mean, '.--');
legend('train', 'val');
xlabel('training epoch'); ylabel('avg. energy');
title('objective');
grid on;

fprintf('\nAverage train error:\n\ttop-1: %f\n\ttop-5: %f\n', trn_top1_mean(end), trn_top5_mean(end));
fprintf('Average validation error:\n\ttop-1: %f\n\ttop-5: %f\n', val_top1_mean(end), val_top5_mean(end));
fprintf('Average objective:\n\ttrain: %f\n\tvalidation: %f\n', trn_obj_mean(end), val_obj_mean(end));