#!/usr/bin/env Rscript
system('clear')

library(caret)

source("load_data.r")
# un-important featuers in dataset obtained from `<tree.r>::varImp()`
# exclude_columns <- c(
# 		'foreign_worker',
# 		'own_telephone',
# 		'num_dependents',
# 		'existing_credits',
# 		'residence_since',
# 		'job',
# 		'other_payment_plans',
# 		'housing',
# 		'cc_age'
# 	)

# dataset 	<- dataset[, -which(names(dataset) %in% exclude_columns)]
# trainset 	<- trainset[, -which(names(trainset) %in% exclude_columns)]
# testset 	<- testset[, -which(names(testset) %in% exclude_columns)]

source("nb.r")
source("svm.r")
source("tree.r")
source("assoc.r")

cat("\n-----------------------------------------------------------\n\n")
cat("> Testing the models\n")

cat("$ Making predictions from models\n")
# By making probabilistic prediction we may be able to use 
# Choquet Integral as combiner predictor
nb_prediction 	<- predict(nb_model, testset[,1:(ncol(testset) -1)], type='raw')

svm_prediction 	<- predict(svm_model, testset[,1:(ncol(testset) -1)], type='probabilities')

tr_prediction 	<- predict(tree_model, testset[,1:(ncol(testset) -1)], type='prob')

# Defining the \mu indicator
# Interpretation:
#		How a combination of primary predictors is sufficient to 
# 				make a judgement on instance's class!?
mu <- new.env()
# 1: NB Classifier
# 2: SVM Classifier
# 3: TREE Classifier
mu[[","]] 		<- 0
mu[["1,2,3"]] 	<- 1		# All 		Classifiers
mu[["1"]] 		<- 0.33		# Only NB 	Classifier
mu[["2"]] 		<- 0.33		# Only SVM 	Classifier
mu[["3"]] 		<- 0.33		# Only DT 	Classifier
mu[["1,2"]] 	<- 0.8		# NB & SVM 	Classifiers
mu[["2,3"]] 	<- 0.66		# SVM & DT 	Classifiers
mu[["1,3"]] 	<- 0.66		# NB & DT 	Classifiers

cat("$ Combining predictions using Fuzzy Choquet Integral\n")

# Predictions predicted by FCL(Fuzzy Choquet Integral)
fci_predictions <- c()

for(i in 1:nrow(testset)) {
	# join the predictions for current instance
	class 	<- rbind(nb_prediction[i,], svm_prediction[i,], tr_prediction[i,])
	# sort all predictions
	gs 		<- apply(class, 2, sort, index.return=TRUE)$good
	tau 	<- c(gs$ix)
	x_tau 	<- as.matrix(c(0, gs$x), dimnames = NULL, rownames = NULL)
	choquet_val <- 0.0;
	# FCL's Sum
	for(classifier in 1:length(tau)) {
		# signature of classifiers in A(i) in FCI
		collaborative_sig = paste(sort(tau[classifier:length(tau)]), collapse = ',')
		# sum current item in sum-iteration
		choquet_val <- choquet_val + 
			(x_tau[classifier + 1] - x_tau[classifier]) * mu[[collaborative_sig]];
	}
	# make a prediction based on FCI's value
	fci_predictions <- cbind(fci_predictions, ifelse(choquet_val > 0.5, "good", "bad"))
}

cat("\n-----------------------------------------------------------\n")
cat("        [Naive Bayesian model contingency table]\n\n")

nb_prediction 	<- predict(nb_model, testset[,1:(ncol(testset) -1)], type='class')

confusionMatrix(nb_prediction, testset[,ncol(testset)])

cat("\n-----------------------------------------------------------\n")
cat("              [SVM model contingency table]\n\n")

svm_prediction 	<- predict(svm_model, testset[,1:(ncol(testset) -1)])

confusionMatrix(svm_prediction, testset[,ncol(testset)])

cat("\n-----------------------------------------------------------\n")
cat("            [Decision tree contingency table]\n\n")

tr_prediction 	<- predict(tree_model, testset[,1:(ncol(testset) -1)], type='class')

confusionMatrix(tr_prediction, testset[,ncol(testset)])

cat("\n-----------------------------------------------------------\n")
cat("       [Overall classification contingency table]\n\n")

confusionMatrix(fci_predictions, testset[,ncol(testset)])