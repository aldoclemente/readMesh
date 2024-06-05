rm(list=ls())
if(system.file(package = "fdaPDE") == "") install.packages("fdaPDE")
source("read.mesh.R")
library(fdaPDE)

# Square
filename <-"data/square.mesh"
domain <- read.mesh(filename)
mesh <- create.mesh.2D(nodes=domain$nodes, triangles = domain$elements, 
                       segments = domain$faces)
plot(mesh)
range( domain$facesmarkers )
range( domain$nodesmarkers )

n = c(16,32,64,128)
for(i in 1:length(n)){
  filename <- paste0("data/square_",n[i],".mesh")
  domain <- read.mesh(filename)
  mesh <- create.mesh.2D(nodes=domain$nodes, triangles = domain$elements)
  
  folder <- paste0("data/unit_square_",n[i],"/")
  if(!dir.exists(folder)) dir.create(folder)
  write.mesh(mesh, folder=folder)
}
# Surface
filename <-"data/sphere.surface_0.125.mesh"
domain <- read.mesh(filename)

range(domain$facesmarkers) # ha senso
range(domain$nodesmarkers) # per imporre dirichlet ? 

mesh <- create.mesh.2.5D(nodes=domain$nodes, triangles = domain$elements)
plot(mesh)
dim(mesh$nodes)
range(domain$nodesmarkers)

h = 0.25 / 2^(0:3)
for(i in 1:length(h)){
  filename <- paste0("data/sphere.surface_",h[i],".mesh")
  domain <- read.mesh(filename)
  mesh <- create.mesh.2.5D(nodes=domain$nodes, triangles = domain$elements)
  
  folder <- paste0("data/sphere_surface_",gsub("[.]", "_",h)[i],"/")
  if(!dir.exists(folder)) dir.create(folder)
  write.mesh(mesh, folder=folder)
}

# Cube
filename <-"data/cube.mesh"
domain <- read.mesh(filename)
mesh <- create.mesh.3D(nodes=domain$nodes, tetrahedrons = domain$elements,
                       nodesattributes = domain$boundary)
plot(mesh)

# Sphere3D
filename <-"data/sphere3D.mesh"
domain <- read.mesh(filename)
mesh <- create.mesh.3D(nodes=domain$nodes, tetrahedrons = domain$elements,
                       nodesattributes = domain$boundary)
plot(mesh)

# 2d write freefem mesh -----------------------------------------------
n <- c(16,32)
i = 1
filename <- paste0("data/square_",n[i],".mesh")
domain <- read.mesh(filename)
range(domain$facesmarkers) # ha senso
range(domain$nodesmarkers) # per imporre dirichlet ? 

mesh <- create.mesh.2D(nodes=domain$nodes, triangles = domain$elements)
plot(mesh)

write_freefem(mesh)

write.table(mesh$nodes[1:10,], file="input_locations.txt", 
            row.names = FALSE, col.names = FALSE)

