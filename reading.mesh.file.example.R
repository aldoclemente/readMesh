rm(list=ls())
source("read.mesh.R")
library(fdaPDE)

# Square
filename <-"data/square.mesh"
domain <- read.mesh(filename)
mesh <- create.mesh.2D(nodes=domain$nodes, triangles = domain$elements)
plot(mesh)

# Surface
filename <-"data/sphere.surface.mesh"
domain <- read.mesh(filename)
mesh <- create.mesh.2.5D(nodes=domain$nodes, triangles = domain$elements)
plot(mesh)

# Cube
filename <-"data/cube.mesh"
domain <- read.mesh(filename)
mesh <- create.mesh.3D(nodes=domain$nodes, tetrahedrons = domain$elements)
plot(mesh)

# Sphere3D
filename <-"data/sphere3D.mesh"
domain <- read.mesh(filename)
mesh <- create.mesh.3D(nodes=domain$nodes, tetrahedrons = domain$elements)
plot(mesh)

