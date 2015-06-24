# splitdf function will return a list of training and testing sets
splitdf 	<- function(dataframe, training_portion = 0.7, seed=NULL) {
	if (!is.null(seed)) set.seed(seed)
	index 		<- 1:nrow(dataframe)
	tindex 		<- sample(index, trunc(length(index) * training_portion))
	trainset 	<- dataframe[tindex, ]
	testset 	<- dataframe[-tindex, ]
	list(trainset = trainset, testset = testset)
}

iprintf 	<- function(...) invisible(cat(sprintf(...)))
isize 		<- function(a) iprintf("[%i, %i]\n", nrow(a), ncol(a)) 

library(foreign)

cat("> Loading dataset\n")
if(!(file.exists("data/German credit fraud dataset/trainset.arff") && 
	 file.exists("data/German credit fraud dataset/testset.arff"))) {

	dataset 	<- read.arff('data/German credit fraud dataset/credit_fruad.arff')

	cat("> Shuffling dataset\n")
	for(iter in 1:100) 
		dataset <- dataset[sample(nrow(dataset)), ]

	cat("> Building train/test sets\n")
	splits 		<- splitdf(dataset, training_portion = 0.7, seed=808)
	trainset 	<- splits$trainset
	testset 	<- splits$testset

	write.arff(trainset, "data/German credit fraud dataset/trainset.arff")
	write.arff(testset, "data/German credit fraud dataset/testset.arff")
} else {
	cat("> Reading from cache sets!\n")
	trainset 	<- read.arff('data/German credit fraud dataset/trainset.arff')
	testset 	<- read.arff('data/German credit fraud dataset/testset.arff')
	dataset 	<- rbind(trainset, testset)
}