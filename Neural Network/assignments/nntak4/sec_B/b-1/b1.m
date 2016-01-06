fprintf('Building 1-NN model....');
knnm = fitcknn(trainset(:, 1:end-1), trainset(:, end), 'NumNeighbors', 1);
fprintf('[DONE]\nTesting the model...');
pred_b1 = predict(knnm, testset(:, 1:end-1));
fprintf('[DONE]\n');
C = confusionmat(pred_b1, testset(:, end));
fprintf('Confusion Matrix:\n');
disp(dataset({C, 'pred_1', 'pred_2'}, 'obsnames', {'real_1', 'real_2'}))

clear C