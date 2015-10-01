function Cost = SetFitness(Pop)
    Cost = zeros(size(Pop, 1), 1);
    for i=1:size(Pop, 1)
        Cost(i) = ackley(Pop(i, :), 20, 0.2, 2 * pi);
    end
end