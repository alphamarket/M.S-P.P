suppressMessages(library(arules))

min_conf = 0.90
min_supp = 0.25

iprintf("> Building `Association` models with minimun(conf = %.2f, supp = %.2f)\n", min_conf, min_supp)

cat("\n------------------------- apriori -------------------------\n")

rules 	<- apriori(dataset[,grep('factor', sapply(dataset, class))], parameter = list(conf = min_conf, supp = min_supp))

cat("\n-----------------------------------------------------------\n\n")

cat("> Association rules extracted\n")

cat("> Pruning association rules by target variable(class)\n")

assoc_model <- subset(rules, subset = rhs %pin% "class")

xrules <-  as(assoc_model, "data.frame")

counter <- 0
pruned_rules <- list()
pruned_rules_counter = 0
for(i in xrules[,1]) {
	counter <- counter + 1
	x <- strsplit(i, " => ", TRUE)[[1]]
	if(grepl("class=(good|bad)",x[2])) {
		pruned_rules_counter <- pruned_rules_counter + 1
		pruned_rules[[pruned_rules_counter]] <- xrules[counter, ]
	}
}

iprintf("> Total %i associated rules with target varibale(class) have been extracted. (%%%.2f of total associated rules)\n", length(assoc_model), length(assoc_model) * 100 / length(rules))

cat("> Rules :\n")

for(rule in pruned_rules)
	iprintf("\t%s\tconf. = %.2f\tsupp. = %.2f\tlift = %.2f\n", gsub("(Factor w/ [0-9] levels )\"(.*)\"", "$1", rule$rules), rule$confidence, rule$support, rule$lift)