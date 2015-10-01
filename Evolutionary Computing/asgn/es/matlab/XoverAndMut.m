function Child = XoverAndMut(P1, P2, CP, MP, Dim, MinSigma, TawPrim, TN, LowB, HighB)
    rndNum = rand;
    sp = P1;
    if(sp < 0.5), sp = p2; end
    if(rndNum < MP)
        Child = MutatFunc(sp, Dim, TawPrim, TN, MinSigma, LowB, HighB);
    elseif((1 - rndNum) < CP)
        Child = CrossOverFunc(P1, P2, Dim, LowB, HighB);
    else
      Child = sp;
    end
end

function Off = MutatFunc(Parent,Dim,TawPrim,TN,MinSigma,LB, HB)

    TPN = TawPrim * normrnd(0,1);
    
    SigmaN = Parent(Dim+1:2*Dim) * exp(TPN + TN);
    
    for i = 1:Dim, if(SigmaN(i) < MinSigma && SigmaN(i) > -MinSigma), SigmaN(i) = MinSigma; end, end
    
    Child = Parent(1:Dim) + SigmaN * normrnd(0,1);
    
    Child(Child > HB) = HB;
    
    Child(Child < LB) = LB;
    
    Off = [Child , SigmaN];
end


function Off = CrossOverFunc(Parent1,Parent2,Dim,LB, HB)
    
    Off(1:Dim) = (Parent1(1:Dim) * rand + Parent2(1:Dim) * rand);
    
    Off(Off > HB) = HB;
    
    Off(Off < LB) = LB;
    
    for i=(Dim + 1):(2 * Dim)
        if(rand < 0.5)
           Off(i) = Parent1(i);
        else
           Off(i) = Parent2(i);
        end  
    end
end

