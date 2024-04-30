##building co-occurence networks for mutations##

library(tidyverse)
library(cooccur)
library(visNetwork)


library(usethis) 
usethis::edit_r_environ()
##load in presence absence table from UpSet plots ## 

presabs145<-C5_upset6

View(presabs145)


#make first column row names

presabs145_2<-presabs145[,-1]
rownames(presabs145_2)<-presabs145[,1]

co<-print(cooccur(presabs145_2, spp_names=TRUE))

presabs145_2<-as.matrix(presabs145_2)

##attempt using igraph 

#generate co-occurrence matrix with sample x mutation not mutation x sample 

presabs145_3<-t(presabs145_2)

#generate co-occurrence matrix 

co_mat2<-t(presabs145_2)%*% presabs145_2

diag(co_mat2)<-0

write.csv(co_mat2, file="Colony5_Matrix.csv")
          
tree<-ape::read.tree(text="(S20024,(S20034,S20029),((S20023,S20018),(S10019,S20022)));")

plot(tree)
#assign dim names
dimnames(co_mat2)<-list(colnames(presabs145_2), colnames(presabs145_2))

tr<-nj(co_mat2)
plot(tr, "p")
write.tree(tr, file="tree_fixed.txt")

graph_from_adjacency_matrix(co_mat2, mode="upper", weighted=TRUE)
g2<-set_vertex_attr(g2, "v_weight", value=colSums(presabs145_3))

t<-sample_smallworld(co_mat2, children=5, mode="undirected")

plot(t, vertex.size=10)

tr<-make_tree(co_mat2, children=5, mode="undirected")

plot(tr, vertex.size=10)
#set diagnoal values to 0 
diag(co_mat2) <- 0

#assign dim names
dimnames(co_mat(list()))

#create graph from adjacency matrix 
# edge weights are equal to frequency of co-occurence

graph_from_adjacency_matrix(co_mat2, mode="upper", weighted=T, diag=FALSE )

g<-set_vertex_attr(g, "v_weight", value=colSums(presabs145_2))

plot(g, vertex.size=V(g)$v_weight*.1, edge.widgth=E(g)$weight*.5)


                                                            

