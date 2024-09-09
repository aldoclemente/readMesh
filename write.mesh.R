write_freefem = function(mesh, filename="mesh.mesh"){
  system(paste0("touch ", filename))
  system(paste0("echo 'MeshVersionFormatted 0 ' >", filename))
  system(paste0("echo ' ' >>", filename))
  system(paste0("echo 'Dimension 2' >>", filename))
  system(paste0("echo ' ' >>", filename))
  system(paste0("echo 'Vertices' >>", filename))
  system(paste0("echo '", nrow(mesh$nodes),"' >>", filename))
  
  boundary = as.integer(mesh$nodesmarkers)
  write.table(cbind(mesh$nodes, boundary), file = filename,
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  system(paste0("echo ' ' >>", filename))
  system(paste0("echo 'Edges' >>", filename))
  system(paste0("echo '", nrow(mesh$edges),"' >>", filename))  
  boundary_edges = as.integer(mesh$edgesmarkers)
  write.table(cbind(mesh$edges, boundary_edges), file = filename,
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  system(paste0("echo ' ' >>", filename))
  system(paste0("echo 'Triangles' >>", filename))
  system(paste0("echo '", nrow(mesh$triangles),"' >>", filename))  
  write.table(cbind(mesh$triangles, rep(0,nrow(mesh$triangles))), file = filename,
              append = TRUE, col.names = FALSE, row.names = FALSE)
  system(paste0("echo ' ' >>", filename))
  system(paste0("echo 'End' >>", filename))
}

write_txt = function(mesh, folder){
    write.table(mesh$nodes,  row.names = FALSE, col.names = FALSE,
                file=paste0(folder,"points.txt"))
    write.table(mesh$neigh,  row.names = FALSE, col.names = FALSE,
                file=paste0(folder,"neigh.txt"))
    write.table(mesh$nodesmarkers,  row.names = FALSE, col.names = FALSE,
                file=paste0(folder,"boundary.txt"))
    
    if(class(mesh) == "mesh.2D" | class(mesh) == "mesh.2.5D"){
      elements <- mesh$triangles
      edges <- mesh$edges
    }else if(class(mesh) == "mesh.3D"){
      elements <- mesh$tetrahedrons
      edges <- mesh$faces
    }
    write.table(elements,  row.names = FALSE, col.names = FALSE,
                file=paste0(folder,"elements.txt"))
    write.table(edges,  row.names = FALSE, col.names = FALSE,
                file=paste0(folder,"edges.txt"))
}
