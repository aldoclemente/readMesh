rm(list=ls())
if(system.file(package = "fdaPDE") == "") install.packages("fdaPDE")
source("read.mesh.R")
library(fdaPDE)

# Square
filename <-"data/square.mesh"
domain <- read.mesh(filename)
mesh <- create.mesh.2D(nodes=domain$nodes, triangles = domain$elements, 
                       segments = domain$faces)
plot(mesh, pch=".")
points(mesh$nodes[domain$nodesmarkers==1,], col="red", pch=16)
points(mesh$nodes[domain$nodesmarkers==2,], col="blue", pch=16)
points(mesh$nodes[domain$nodesmarkers==3,], col="pink", pch=16)
points(mesh$nodes[domain$nodesmarkers==4,], col="green3", pch=16)
mtext(side=1, text="1", col="red")
mtext(side=2, text="4", col="green3", las=2)
mtext(side=3, text="3", col="pink")
mtext(side=4, text="2", col="blue",las=2)



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
mesh <- create.mesh.3D(nodes=domain$nodes, 
                       tetrahedrons = domain$elements, 
                       nodesattributes = domain$boundary)
plot(mesh)

#nodes <- mesh$nodes

# da documentazione FreeFem
# faces == 1 -> y = 0
# faces == 2 -> x = 1
# faces == 3 -> y = 1
# faces == 4 -> x = 0
# faces == 5 -> z = 0
# faces == 6 -> z = 1
# Dobbiamo farci passare anche la numerazione delle facce!!!
faces1 <- as.vector(t(domain$faces[domain$facesmarkers==1,])) 
faces2 <- as.vector(t(domain$faces[domain$facesmarkers==2,]))
faces3 <- as.vector(t(domain$faces[domain$facesmarkers==3,]))
faces4 <- as.vector(t(domain$faces[domain$facesmarkers==4,]))
faces5 <- as.vector(t(domain$faces[domain$facesmarkers==5,]))
faces6 <- as.vector(t(domain$faces[domain$facesmarkers==6,]))

#aux_mesh <- create.mesh.2.5D(nodes=nodes, triangles=x$faces[x$facesmarkers,], order=1)
#edges <- as.vector(t(aux_mesh$edges))

{
  open3d()
  axes3d()
  pop3d("lights")
  light3d(specular="black")
  col = c("red", "blue", "pink", "green3", "white", "gray")
  
  triangles3d(mesh$nodes[faces1,1],mesh$nodes[faces1,2],mesh$nodes[faces1,3],col=col[1])
  triangles3d(mesh$nodes[faces2,1],mesh$nodes[faces2,2],mesh$nodes[faces2,3],col=col[2])
  triangles3d(mesh$nodes[faces3,1],mesh$nodes[faces3,2],mesh$nodes[faces3,3],col=col[3])
  triangles3d(mesh$nodes[faces4,1],mesh$nodes[faces4,2],mesh$nodes[faces4,3],col=col[4])
  triangles3d(mesh$nodes[faces5,1],mesh$nodes[faces5,2],mesh$nodes[faces5,3],col=col[5])
  triangles3d(mesh$nodes[faces6,1],mesh$nodes[faces6,2],mesh$nodes[faces6,3],col=col[6])
  aspect3d("iso")
}
#points3d(nodes[,1], nodes[,2], nodes[,3], col="black", ...)
# ... mmm ...
table(domain$nodesmarkers)
mesh$nodes[domain$nodesmarkers==1,] # yz, x = 0
mesh$nodes[domain$nodesmarkers==2,] # yz, x = 1
mesh$nodes[domain$nodesmarkers==4,] # xz, y = 0
mesh$nodes[domain$nodesmarkers==8,] # xz, y = 1
mesh$nodes[domain$nodesmarkers==16,] # xy, z = 0
mesh$nodes[domain$nodesmarkers==32,] # xy, z = 1

mesh$nodes[domain$nodesmarkers==5,] # x=0, y=0
mesh$nodes[domain$nodesmarkers==6,] # x=1, y=0

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

