# Package
library("support.CEs")


# Example 1: Unlabeled Design
des1 <- rotation.design(attribute.names = list(
 Region = c("Reg_A", "Reg_B", "Reg_C"), 
 Eco = c("Conv.", "More", "Most"), 
 Price = c("1", "1.1", "1.2")), 
 nalternatives = 2, nblocks = 1, row.renames = FALSE, 
 randomize = TRUE, seed = 987)
des1

questionnaire(choice.experiment.design = des1)

data("syn.res1")
syn.res1[1:3,]

desmat1 <- make.design.matrix(choice.experiment.design = des1, 
 optout = TRUE, 
 categorical.attributes = c("Region", "Eco"), 
 continuous.attributes = c("Price"),
 unlabeled = TRUE)
desmat1[1:3,]

dataset1 <- make.dataset(respondent.dataset = syn.res1, 
 choice.indicators = 
 c("q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8", "q9"), 
 design.matrix = desmat1)
dataset1[1:3,]

clogout1 <- clogit(RES ~ ASC + Reg_B + Reg_C + More + Most + 
 More:F + Most:F + Price + strata(STR), data = dataset1)
clogout1
gofm(clogout1)

mwtp(output = clogout1, monetary.variables = c("Price"), 
 nonmonetary.variables = 
 c("Reg_B", "Reg_C", "More", "Most", "More:F", "Most:F"), 
 percentile.points = c(5, 95), seed = 987)


# Example 2: Labeled Design
des2 <- Lma.design(attribute.names = list(
 Eco = c("Conv.", "More", "Most"), 
 Price = c("1", "1.1", "1.2")), 
 nalternatives = 3, nblocks = 2, row.renames = FALSE, seed = 987)
des2

questionnaire(choice.experiment.design = des2)

data("syn.res2")
syn.res2[1:3,]

desmat2 <- make.design.matrix(choice.experiment.design = des2, 
 optout = TRUE, 
 categorical.attributes = c("Eco"), 
 continuous.attributes = c("Price"), 
 unlabeled = FALSE)
desmat2[1:4,]

dataset2 <- make.dataset(respondent.dataset = syn.res2, 
 choice.indicators = 
 c("q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8", "q9"), 
 design.matrix = desmat2)
dataset2[1:4,]

clogout2 <- clogit(RES ~ ASC1 + More1 + Most1 + Price1 + 
 ASC2 + More2 + Most2 + Price2 + ASC3 + More3 + Most3 + Price3 + 
 strata(STR), data = dataset2)
clogout2
gofm(clogout2)

mwtp(output = clogout2, 
 monetary.variables = c("Price1", "Price2", "Price3"), 
 nonmonetary.variables = list(
 c("More1", "Most1"), c("More2", "Most2"), c("More3", "Most3")), 
 percentile.points = c(5, 95), seed = 987)

