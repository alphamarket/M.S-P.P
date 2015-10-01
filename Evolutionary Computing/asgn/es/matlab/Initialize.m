clc
clear all
format shortE
format loose

% Algorithm Parameters
PopSize = 20;
ChildSize = 200;
MaxGen = 500;
CP = 0.5;
MP = 0.5;
MinSigma =  1e-5;

% Problem Statement
Dim = 30;
LowB = -5;
HighB = 5;

% Initial Population
Pop(:, 1:Dim) = rand(PopSize, Dim) * (HighB - LowB) + LowB;
Pop(:, (Dim + 1):(2 * Dim)) = rand(PopSize, Dim);

Cost        = SetFitness(Pop(:, 1:Dim));

[Cost, inx] = sort(Cost);
Pop         = Pop(inx,:);

TawPrim     = 1 / sqrt(2 * Dim);
Taw         = 1 / sqrt(2 * sqrt(Dim));
TN          = Taw * normrnd(0,1);