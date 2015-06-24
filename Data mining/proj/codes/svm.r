library(kernlab)

cat("> Building `SVM` model\n")

svm_model 	<- ksvm(class~., data = trainset, kernel = "rbfdot", prob.model=TRUE)