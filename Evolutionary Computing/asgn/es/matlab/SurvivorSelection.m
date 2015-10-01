% mu & lambda
% elit selection
NewPop = [NewPop; Pop(1, :)];
% calc. fitness
Cost = SetFitness(NewPop(:, 1:Dim));
% sort fitnesses
[Cost, inx] = sort(Cost);
% sort pop based on fitnesses order
NewPop = NewPop(inx, :);
% fetch top solutions
Pop = NewPop(1:PopSize, :);