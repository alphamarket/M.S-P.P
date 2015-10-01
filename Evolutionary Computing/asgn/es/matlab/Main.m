Initialize;

BestCost = inf;

BestIter = 0;

for Iter = 1:MaxGen
    
    NewPop = [];
    
    for i = 1:ChildSize-1
        P1 = Pop(round(rand*(PopSize-1))+1,:);
        P2 = Pop(round(rand*(PopSize-1))+1,:);
        Child = XoverAndMut(P1, P2, CP, MP, Dim, MinSigma, TawPrim, TN, LowB, HighB);
        NewPop = [NewPop; Child]; %#ok<AGROW>
    end
    
    SurvivorSelection;
    
    disp([Iter Cost(1)]);
    
    if BestCost > Cost(1)
        BestCost = Cost(1);
        BestIter = Iter;
        BestSolution = Pop(1,1:Dim);
    end
end

BestIter

BestCost

BestSolution