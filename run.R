source("EvolutionDif.R")

datos <- read.csv("test1.csv")

ejeED <- DE(datos, 66)
#cluster1 <- ejeED[[1]]
#cluster2 <- ejeED[[2]]
#cluster3 <- ejeED[[3]]
#cluster4 <- ejeED[[4]]
#jaccard <- vectorJaccard(cluster1, cluster2, cluster3, cluster4)
#jaccard

ejeED <- DE(datos, 66)
#cluster1 <- ejeED[[1]]
#cluster2 <- ejeED[[2]]
#cluster3 <- ejeED[[3]]
#cluster4 <- ejeED[[4]]
#jaccard <- vectorJaccard(cluster1, cluster2, cluster3, cluster4)
#jaccard
