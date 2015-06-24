library(e1071)
library(class)

cat("> Building `Naive Bayesian` model\n")

nb_model 	<- naiveBayes(trainset[,1:(ncol(trainset) - 1)], trainset[,ncol(trainset)])