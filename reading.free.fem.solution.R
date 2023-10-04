rm(list=ls())
source("read.mesh.R")
library(fdaPDE)

# reading free fem mesh
filename <-"data/square_32.mesh"
domain <- read.mesh(filename)
mesh <- create.mesh.2D(nodes=domain$nodes, triangles = domain$elements)
plot(mesh)

# exporting free fem mesh 
# write.mesh(mesh)

nrow(mesh$nodes)
# loading free fem solution (script/heatEquation.edp) at time T
coeff <- read.table("script/solution_T.txt", header = F)

FEMbasis <- create.FEM.basis(mesh)
plot(FEM(coeff,FEMbasis))
