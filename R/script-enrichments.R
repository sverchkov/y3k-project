# Using Rolemodel for gene ontology
library("Rolemodel")
# The yeast database is org.Sc.sgd.db
library("org.Sc.sgd.db")

# Assume we have the attachments.by.action object

gene.set <- attachments.by.action[[1]]

# ISSUE: Rolemodel can't deal with the gene names I have
# TODO: Remap into something it can use

res <- rmTable( gene.set, lib = "org.Sc.sgd", n.upp = 20, n.low = 1, nupstart = 10, by = 1 )