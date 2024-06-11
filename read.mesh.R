# Reads .mesh file (FreeFem++)

strtodouble <- function(x){
  as.numeric(unlist(strsplit(x, " ")))
}

strtointeger <- function(x){
  as.integer(unlist(strsplit(x, " ")))
}

read_freefem <- function(filename){
  nodes            <- NULL
  edges            <- NULL
  edgesmarkers     <- NULL
  triangle         <- NULL
  trianglesmarkers <- NULL
  tetrahedra       <- NULL
  boundary         <- NULL
  
  file <- readLines(filename)
  
  format.version <- as.integer(unlist(strsplit(file[grep("MeshVersion", file, fixed=T)]," "))[2])
  
  if(format.version == 1){
    #2D 
    embedding_dimension <- strtoi(file[4])
  }else if(format.version == 0 || format.version == 2){
    
    # 2D, 2.5D or 3D
    embedding_dimension <- strtoi(strsplit(file[3],split=" ")[[1]][2])
    
  }
  
  # setting local dimension
  local_dimension <- as.integer(ifelse( identical(which(file=="Tetrahedra"),integer(0)), 2, 3))
  
  # Nodes and boundary markers
  # Nodes
  idx.nodes <- grep("Vertices", file, fixed = T)
  num_nodes <- strtoi(file[idx.nodes+1])
  nodes <- matrix(strtodouble(file[(idx.nodes+2):(idx.nodes+1+num_nodes)]),
                  nrow=num_nodes, ncol=(embedding_dimension+1), byrow=T)
  
  # nodes boundary markers
  nodesmarkers <- as.matrix(nodes[, embedding_dimension+1])
  nodes <- nodes[,1:embedding_dimension]
  
  # Triangles
  idx.triangles <- grep("Triangles", file, fixed = T)
  num_triangles <- strtoi(file[idx.triangles+1])
  triangles <- matrix(strtodouble(file[(idx.triangles+2):(idx.triangles+1+num_triangles)]),
                      nrow=num_triangles, ncol=4, byrow=T)
  
  if((local_dimension == 2 & embedding_dimension == 2) |
     (local_dimension == 2 & embedding_dimension == 3)){
    idx.edges <- which(file == "Edges")
    num_edges <- strtoi(file[idx.edges+1]) 
    edges <- matrix(strtodouble(file[(idx.edges+2):(idx.edges+1+num_edges)]),
                    nrow=num_edges, ncol=3, byrow=T)
    edgesmarkers <- as.matrix(edges[,3])
    edges <- edges[,1:2]
  }else if(local_dimension == 3 & embedding_dimension ==3){
    # Tetrahedra
    idx.tetrahedra <- which(file == "Tetrahedra")
    num_tetrahedra <- strtoi(file[idx.tetrahedra+1]) 
    tetrahedra <- matrix(strtodouble(file[(idx.tetrahedra+2):(idx.tetrahedra+1+num_tetrahedra)]),
                         nrow=num_tetrahedra, ncol=5, byrow=T)
    tetrahedra <- tetrahedra[,1:4]
    
    trianglesmarkers <- as.matrix(triangles[,4])
  }
  
  triangles <- triangles[,1:3]
  
  storage.mode(nodesmarkers) <- "integer" 
  storage.mode(triangles) <- "integer"
  storage.mode(trianglesmarkers) <- "integer"
  storage.mode(edges) <- "integer"
  storage.mode(edgesmarkers) <-"integer"
  storage.mode(tetrahedra) <- "integer"
  
  if(local_dimension == 2){
    return(list(nodes=nodes, 
                elements=triangles,
                faces = edges,
                facesmarkers = edgesmarkers,
                nodesmarkers=nodesmarkers))
  }else{
    tetrahedrons <- mesh$tetrahedrons
    storage.mode(tetrahedrons) <- "integer"
    return(list(nodes=nodes, 
                elements=tetrahedra,
                faces = triangles,
                facesmarkers = trianglesmarkers,
                nodesmarkers=nodesmarkers))
  }
}

read.mesh <-function(filename){
  
  # Reading file
  # check .mesh :)
  if( grep(".mesh", filename, fixed=T)){
    read_freefem(filename)
  }
}

# read.mesh <-function(filename){
#   
#   # Reading file
#   file <- readLines(filename)
#   
#   format.version <- strtoi(strsplit(file[1],split=" ")[[1]][2])
#   
#   if(format.version == 1){
#     
#     #2D 
#     embedding_dimension <- strtoi(file[4])
#   }else if(format.version == 2){
#     
#     # 2.5D or 3D
#     embedding_dimension <- strtoi(strsplit(file[3],split=" ")[[1]][2])
#     
#   }
#   
#   # setting local dimension
#   local_dimension <- as.integer(ifelse( identical(which(file=="Tetrahedra"),integer(0)), 2, 3))
#   
#   # Nodes and boundary markers
#   idx.nodes <- which(file == "Vertices")
#   num_nodes <- strtoi(file[idx.nodes+1])
#   nodes <- matrix(nrow=num_nodes, ncol=embedding_dimension)
#  
#   boundary <- matrix(nrow=num_nodes, ncol=1)
#   
#   for( i in (idx.nodes+2):(idx.nodes+1+strtoi(file[idx.nodes+1]))){
#     line <- strsplit(file[i], " ")[[1]]
#     nodes[i-(idx.nodes+1),] <- as.numeric(line)[1:embedding_dimension] # coordinates
#     boundary[i-(idx.nodes+1)] <- as.numeric(line)[(embedding_dimension+1)] # boundary markers -> label
#   }
#   
#   # Tringles
#   idx.triangles <- which(file == "Triangles")
#   num_triangles <- strtoi(file[idx.triangles+1]) 
#   triangles <- matrix(nrow=num_triangles, ncol=3)
#   for( i in (idx.triangles+2):(idx.triangles+1+strtoi(file[idx.triangles+1]))){
#     line <- strsplit(file[i], " ")[[1]]
#     triangles[i-(idx.triangles+1),] <- as.integer(line)[1:3]
#   }
#   
#   if(local_dimension == 3 & embedding_dimension ==3){
#     # Tetrahedra
#     idx.tetrahedra <- which(file == "Tetrahedra")
#     num_tetrahedra <- strtoi(file[idx.tetrahedra+1]) 
#     tetrahedra <- matrix(nrow=num_tetrahedra, ncol=(local_dimension+1))
#     for( i in (idx.tetrahedra+2):(idx.tetrahedra+1+strtoi(file[idx.tetrahedra+1]))){
#       line <- strsplit(file[i], " ")[[1]]
#       tetrahedra[i-(idx.tetrahedra+1),] <- as.integer(line)[1:(local_dimension+1)]
#     }
#   }
#   
#   storage.mode(boundary) <- "integer" 
#   
#   if( embedding_dimension == 2 & local_dimension == 2){
#     mesh <- fdaPDE::create.mesh.2D(nodes=nodes, triangles=triangles, 
#                                    nodesattribute = as.integer(boundary != 0))
#   }else if(embedding_dimension==3 & local_dimension == 2){
#     mesh <- fdaPDE::create.mesh.2.5D(nodes=nodes, triangles=triangles,
#                                      nodesattribute = as.integer(boundary != 0))
#   }else if(embedding_dimension==3 & local_dimension == 3){
#     mesh <- fdaPDE::create.mesh.3D(nodes=nodes, tetrahedrons = tetrahedra,
#                                    nodesattribute = as.integer(boundary != 0)) 
#   }
#   nodes <- mesh$nodes
#   elements <- mesh$elements
#   neigh <- mesh$neighbors
#   #boundary <- mesh$nodesmarkers
#   triangles <- mesh$triangles 
#   
#   storage.mode(elements) <- "integer"
#   storage.mode(neigh)    <- "integer"
#   boundary <- matrix((boundary != 0))
#   storage.mode(boundary) <- "integer"
#   storage.mode(triangles)<- "integer"
#   
#   if(local_dimension == 2){
#     return(list(nodes=nodes, 
#                 elements=triangles,
#                 neigh=neigh,
#                 boundary=boundary))
#   }else{
#     tetrahedrons <- mesh$tetrahedrons
#     storage.mode(tetrahedrons) <- "integer"
#     return(list(nodes=nodes, 
#                 elements=tetrahedrons,
#                 neigh=neigh,
#                 boundary=boundary))
#   }
# }

# x = fdaPDE mesh 
write.mesh <- function(x, folder){
  write.csv(x$nodes, file=paste0(folder,"points.csv"))
  write.csv(x$neigh, file=paste0(folder,"neigh.csv"))
  write.csv(x$nodesmarkers, file=paste0(folder,"boundary.csv"))
  
  if(class(x) == "mesh.2D" | class(x) == "mesh.2.5D"){
    elements <- x$triangles
    edges <- x$edges
  }else if(class(x) == "mesh.3D"){
    elements <- x$tetrahedrons
    edges <- x$faces
  }
  write.csv(elements, file=paste0(folder,"elements.csv"))
  write.csv(edges, file=paste0(folder,"edges.csv"))
}


read_mesh <- function(filename) {
  # read data
  dat <- readLines(filename)
  
  # check file format version
  gmshver <- strtodouble(dat[grep("$MeshFormat", dat, fixed = T)+1])[1]
  if (floor(gmshver) != 2)
    stop("Only MSH file format version 2 is supported!", call. = F)
  
  # get node points
  inodes <- grep("$Nodes", dat, fixed = T)
  npoin <- as.numeric(dat[inodes + 1])
  nodes <- sapply(dat[seq(inodes + 2, inodes + npoin + 1)], tonum, USE.NAMES = F)
  
  # analyse elements
  ielem <- grep("$Elements", dat, fixed = T)
  nelem_t <- as.numeric(dat[ielem + 1])
  jelem <- seq(ielem + 2, ielem + nelem_t + 1)
  types <- sapply(dat[jelem], function(x) tonum(x)[2], USE.NAMES = F)
  if (any(types > 2))
    stop("Can only handle element types 'line' and 'triangle' in MSH file!", call. = F)
  
  # get boundary nodes (ipobo)
  nbnd <- length(which(types == 1))
  ipobo <- sapply(dat[jelem[types == 1]], function(x) as.integer(tonum(x)[1]), USE.NAMES = F)
  
  # get triangle definitions (ikle)
  nelem <- length(which(types == 2))
  ikle <- sapply(dat[jelem[types == 2]], function(x) as.integer(tail(tonum(x), 3)), USE.NAMES = F)
  
  # output
  return(list(nelem = as.integer(nelem),
              npoin = as.integer(npoin),
              ikle = t(ikle),
              ipobo = ipobo,
              x = nodes[2,],
              y = nodes[3,]))
}