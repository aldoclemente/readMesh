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
