# Import functions and data
library(fdaPDE)
library(simcausal)
library(latex2exp)
library(raster)

SST<-raster("../ParameterCascading/Data/SeaSurfaceTemp10042020.nc", varname="analysed_sst")
SST<-crop(SST,extent(c(-99,-79,18,32)))
OBSER<-readRDS("../ParameterCascading/Data/GoMSSTData.rds")

R_betax<-function(x,y)
{
  load("../ParameterCascading/Data/Interpu.Rdata")
  
  Value<-extract(Interpu,cbind(x,y))
  
  Value
}

R_betay<-function(x,y)
{
  
  load("../ParameterCascading/Data/Interpv.Rdata")
  
  Value<-extract(Interpv,cbind(x,y))
  
  
  Value
}

################################################################################
### Diffusion estimate from actual buoys data - Florida
################################################################################

load("../ParameterCascading/Data/meshfine_cut.Rdata")
load("../ParameterCascading/Data/Poly_cut.Rdata")

# Plot data
SST_flo = crop(SST, extent(-85.275,-78.975,23.305,30.805))
par(mar=c(0,1,0,1))
plot(SST_flo - 273.15, col = viridis::inferno(100), bty="n", box=FALSE, axes=FALSE, xlim=c(-85.275,-78.975), ylim=c(23.305,30.805), legend.args=list(text=TeX('Temperature [°C]'),line=1,cex=1.28))

# First mesh
{
  # Set boundary nodes
  florida_bound = c(0,0)
  for(i in 1:122)
  {
    if(meshfine$nodes[i,1] > -85.3)
    {
      florida_bound = rbind(florida_bound, meshfine$nodes[i,])
    }
  }
  florida_bound = florida_bound[-1,]
  left_bound_y = seq(23.805,29.225,1)
  left_bound_x = rep(-85.275,length(left_bound_y))
  left_bound = cbind(left_bound_x, left_bound_y)
  Sbound = rbind(florida_bound[1:12,],left_bound,florida_bound[13:dim(florida_bound)[1],])
  
  # Select sampled points
  Points<-meshfine$nodes
  interior_pts =OBSER[,1:2]
  data_vals = OBSER[,3]
  OBS<-data.frame(interior_pts,data_vals)
  OBS<-SpatialPointsDataFrame(cbind(OBS$LON,OBS$LAT), OBS, coords.nrs = numeric(0),
                              proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), bbox = NULL)
  SpatialcoordDS <- spTransform(OBS, CRS(proj4string(Poly_cut)))
  SpatialcoordDS<- SpatialcoordDS[Poly_cut, ]
  interior_pts<-SpatialcoordDS@data
  data_vals <-extract(SST,interior_pts[,1:2]) - 273.15
  interior_pts<-interior_pts[which(is.na(data_vals)==FALSE),]
  data_vals<-data_vals[which(is.na(data_vals)==FALSE)]
  
  florida_pts = c(0,0)
  florida_vals = 0
  for(i in 1:dim(interior_pts)[1])
  {
    if(interior_pts$LON[i] > -85.3)
    {
      florida_pts = rbind(florida_pts, c(interior_pts$LON[i],interior_pts$LAT[i]))
      florida_vals = c(florida_vals, data_vals[i])
    }
  }
  florida_pts = florida_pts[-1,]
  florida_vals = florida_vals[-1]
  
  # Define mesh and FEMbasis
  c<-length(Sbound[,1])
  edge<-cbind(matrix(c(2:c-1,c), nrow=c,ncol=1 ),matrix(c(2:c,1), nrow=c, ncol=1))
  Bound_and_nodes<-rbind(as.matrix(Sbound),florida_pts)
  
  mesh_1 <- create.mesh.2D(nodes = Bound_and_nodes, segments = edge, order = 1, verbosity=0)
  mesh<-refine.mesh.2D(mesh_1,minimum_angle=0.25,maximum_area=0.05)
  plot(mesh)
  
  plot(mesh,asp=1, pch=".")
  box()
  points(mesh$nodes[which(mesh$nodesmarkers==0),], pch=16,cex=1.2)
  points(mesh$nodes[which(mesh$nodesmarkers==1),], pch=16, col='red',cex=1.2)
  
  FEMbasis <-create.FEM.basis(mesh)
}

MESH_tomasetto_lasca = mesh
plot(MESH_tomasetto_lasca)


# Regular mesh with refinement for dominant transport
{
  florida_bound = c(0,0)
  for(i in 1:122)
  {
    if(meshfine$nodes[i,1] > -85.3)
    {
      florida_bound = rbind(florida_bound, meshfine$nodes[i,])
    }
  }
  florida_bound = florida_bound[-(1:13),]
  florida_bound = florida_bound[(1:27),]
  
  bound1_x = seq(-81.485,-78.975,len=8)
  bound1_y = rep(30.805,length(bound1_x))
  
  bound2_y = seq(30.805,23.305,len=15)
  bound2_x = rep(-78.975,length(bound2_y))
  bound2_x = bound2_x[-1]
  bound2_y = bound2_y[-1]
  
  bound3_x = seq(-78.975,-85.275,len=20)
  bound3_y = rep(23.305,length(bound3_x))
  bound3_x = bound3_x[-1]
  bound3_y = bound3_y[-1]
  
  bound4_y = seq(23.305,29.225,len=15)
  bound4_x = rep(-85.275,length(bound4_y))
  bound4_x = bound4_x[-1]
  bound4_y = bound4_y[-1]
  
  bounds_x = c(bound1_x, bound2_x, bound3_x, bound4_x)
  bounds_y = c(bound1_y, bound2_y, bound3_y, bound4_y)
  bounds = cbind(bounds_x, bounds_y)
  Sbound = rbind(florida_bound,bounds)
  
  # Select data points
  Points<-meshfine$nodes
  interior_pts =OBSER[,1:2]
  data_vals = OBSER[,3]
  OBS<-data.frame(interior_pts,data_vals)
  OBS<-SpatialPointsDataFrame(cbind(OBS$LON,OBS$LAT), OBS, coords.nrs = numeric(0),
                              proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), bbox = NULL)
  SpatialcoordDS <- spTransform(OBS, CRS(proj4string(Poly_cut)))
  SpatialcoordDS<- SpatialcoordDS[Poly_cut, ]
  interior_pts<-SpatialcoordDS@data
  data_vals <-extract(SST,interior_pts[,1:2]) - 273.15
  interior_pts<-interior_pts[which(is.na(data_vals)==FALSE),]
  data_vals<-data_vals[which(is.na(data_vals)==FALSE)]
  
  florida_pts = c(0,0)
  florida_vals = 0
  for(i in 1:dim(interior_pts)[1])
  {
    if(interior_pts$LON[i] > -85.3)
    {
      florida_pts = rbind(florida_pts, c(interior_pts$LON[i],interior_pts$LAT[i]))
      florida_vals = c(florida_vals, data_vals[i])
    }
  }
  florida_pts = florida_pts[-1,]
  florida_vals = florida_vals[-1]
  
  # Define mesh and FEMbasis
  l<-length(Sbound[,1])
  edge<-cbind(matrix(c(2:l-1,l), nrow=l,ncol=1 ),matrix(c(2:l,1), nrow=l, ncol=1))
  Bound_and_nodes<-rbind(as.matrix(Sbound),florida_pts)
  mesh_1 <- create.mesh.2D(nodes = Bound_and_nodes, segments = edge, order = 1, verbosity=0)
  mesh<-refine.mesh.2D(mesh_1,minimum_angle=20,maximum_area=0.005)
  
  # Rejection sampling to select nodes with dominant transport
  new_nodes = c(0,0)
  internal_nodes = mesh$nodes[which(mesh$nodesmarkers==0),]
  for(i in 1:dim(internal_nodes)[1]){
    beta_mod = sqrt(R_betax(internal_nodes[i,1],internal_nodes[i,2])^2 + R_betay(internal_nodes[i,1],internal_nodes[i,2])^2)
    if(beta_mod > 0.9){
      new_nodes = rbind(new_nodes, internal_nodes[i,])
    }else if(rbern(1, prob = beta_mod/8) == 1){
      new_nodes = rbind(new_nodes, internal_nodes[i,])
    }
  }
  new_nodes = new_nodes[-1,]
  new_nodes = rbind(as.matrix(Sbound), new_nodes)
  mesh_2 <- create.mesh.2D(nodes = new_nodes, segments = edge, order = 1)
  mesh <- refine.mesh.2D(mesh_2, minimum_angle=30, maximum_area = 0.1)
  plot(mesh,asp=1, pch=".")
}

MESH_tomasetto = mesh
LOCATIONS_tomasetto = florida_pts

# per il campo ?
MESH_fine = refine.mesh.2D(MESH_tomasetto_lasca,minimum_angle=0.25,maximum_area=0.025)

plot(MESH_fine)
points(florida_pts, pch=16, col="red")

plot(MESH_tomasetto)
points(florida_pts, pch=16, col="red")

plot(MESH_tomasetto_lasca)
points(florida_pts, pch=16, col="red")

# write freefem
source("read.mesh.R")
source("write.mesh.R")

if(!dir.exists("data/emanuele/")) dir.create("data/emanuele/")

# mesh lasca -------------------------------------------------------------------
foldername = "data/emanuele/lasca/"
if(!dir.exists(foldername)) dir.create(foldername)

write_freefem(MESH_tomasetto_lasca, 
              filename = paste0(foldername, "mesh.mesh"))

write.table(florida_pts, 
            file=paste0(foldername, "locations.txt"), 
            row.names = FALSE, col.names = FALSE)

write.table(florida_vals, 
            file=paste0(foldername, "observations.txt"), 
            row.names = FALSE, col.names = FALSE)

FEMbasis = create.FEM.basis(MESH_tomasetto_lasca)
transport_x = vector(mode="numeric", length=nrow(MESH_tomasetto_lasca$nodes))
transport_y = vector(mode="numeric", length=nrow(MESH_tomasetto_lasca$nodes))

for(i in 1:nrow(MESH_tomasetto_lasca$nodes)){
  transport_x[i] = R_betax(MESH_tomasetto_lasca$nodes[i,1], 
                           MESH_tomasetto_lasca$nodes[i,2])
  
  transport_y[i] = R_betay(MESH_tomasetto_lasca$nodes[i,1], 
                           MESH_tomasetto_lasca$nodes[i,2])
}

write.table(cbind(transport_x, transport_y), 
            file=paste0(foldername, "Beta.txt"), 
            row.names = FALSE, col.names = FALSE )
write.table(transport_x, 
            file=paste0(foldername, "Beta_X.txt"), 
            row.names = FALSE, col.names = FALSE)
write.table(transport_y, 
            file=paste0(foldername, "Beta_Y.txt"), 
            row.names = FALSE, col.names = FALSE)

# mesh fine --------------------------------------------------------------------

foldername = "data/emanuele/fine/"
if(!dir.exists(foldername)) dir.create(foldername)


write_freefem(MESH_fine, 
              filename = paste0(foldername, "mesh.mesh"))

write.table(florida_pts, 
            file=paste0(foldername, "locations.txt"), 
            row.names = FALSE, col.names = FALSE)

write.table(florida_vals, 
            file=paste0(foldername, "observations.txt"), 
            row.names = FALSE, col.names = FALSE)

FEMbasis = create.FEM.basis(MESH_fine)
transport_x = vector(mode="numeric", length=nrow(MESH_fine$nodes))
transport_y = vector(mode="numeric", length=nrow(MESH_fine$nodes))

for(i in 1:nrow(MESH_fine$nodes)){
  transport_x[i] = R_betax(MESH_fine$nodes[i,1], 
                           MESH_fine$nodes[i,2])
  
  transport_y[i] = R_betay(MESH_fine$nodes[i,1], 
                           MESH_fine$nodes[i,2])
}

write.table(cbind(transport_x, transport_y), 
            file=paste0(foldername, "Beta.txt"), 
            row.names = FALSE, col.names = FALSE )
write.table(transport_x, 
            file=paste0(foldername, "Beta_X.txt"), 
            row.names = FALSE, col.names = FALSE)
write.table(transport_y, 
            file=paste0(foldername, "Beta_Y.txt"), 
            row.names = FALSE, col.names = FALSE)

# mesh tomasetto ---------------------------------------------------------------

foldername = "data/emanuele/tomasetto/"
if(!dir.exists(foldername)) dir.create(foldername)

write_freefem(MESH_tomasetto, 
              filename = paste0(foldername,"mesh.mesh"))

write_freefem(MESH_fine, 
              filename = paste0(foldername, "mesh.mesh"))

write.table(florida_pts, 
            file=paste0(foldername, "locations.txt"), 
            row.names = FALSE, col.names = FALSE)

write.table(florida_vals, 
            file=paste0(foldername, "observations.txt"), 
            row.names = FALSE, col.names = FALSE)

FEMbasis = create.FEM.basis(MESH_tomasetto)
transport_x = vector(mode="numeric", length=nrow(MESH_tomasetto$nodes))
transport_y = vector(mode="numeric", length=nrow(MESH_tomasetto$nodes))

for(i in 1:nrow(MESH_tomasetto$nodes)){
  transport_x[i] = R_betax(MESH_tomasetto$nodes[i,1], 
                           MESH_tomasetto$nodes[i,2])
  
  transport_y[i] = R_betay(MESH_tomasetto$nodes[i,1], 
                           MESH_tomasetto$nodes[i,2])
}

write.table(cbind(transport_x, transport_y), 
            file=paste0(foldername, "Beta.txt"), 
            row.names = FALSE, col.names = FALSE )
write.table(transport_x, 
            file=paste0(foldername, "Beta_X.txt"), 
            row.names = FALSE, col.names = FALSE)
write.table(transport_y, 
            file=paste0(foldername, "Beta_Y.txt"), 
            row.names = FALSE, col.names = FALSE)


# mesh completa ----------------------------------------------------------------
foldername = "data/emanuele/golfo/"
if(!dir.exists(foldername)) dir.create(foldername)

library(femR)
library(sf)
boundary_nodes = meshfine$nodes[which(meshfine$nodesmarkers==TRUE),]
boundary_edges = cbind(1:(nrow(boundary_nodes)-1), 2:nrow(boundary_nodes))
boundary_edges = rbind(boundary_edges, c(nrow(boundary_nodes), 1))

domain = femR::Domain(list(nodes=boundary_nodes, 
                           edges = boundary_edges))
plot(st_as_sfc(domain))

mesh_completa <- build_mesh(domain, maximum_area = 0.05, minimum_angle = 25)
plot(st_as_sfc(mesh_completa))
points(cbind(interior_pts$LON, interior_pts$LAT), pch=16, col="red")

mesh_completa <- fdaPDE::create.mesh.2D(nodes=mesh_completa$nodes(),
                                        triangles = mesh_completa$elements())

write_freefem(mesh_completa, 
              filename = paste0(foldername, "mesh.mesh"))

write.table(cbind(interior_pts$LON, interior_pts$LAT), 
            file=paste0(foldername, "locations.txt"), 
            row.names = FALSE, col.names = FALSE)

write.table(interior_pts$data_vals, 
            file=paste0(foldername, "observations.txt"), 
            row.names = FALSE, col.names = FALSE)

FEMbasis = create.FEM.basis(mesh_completa)
transport_x = vector(mode="numeric", length=nrow(mesh_completa$nodes))
transport_y = vector(mode="numeric", length=nrow(mesh_completa$nodes))

for(i in 1:nrow(mesh_completa$nodes)){
  transport_x[i] = R_betax(mesh_completa$nodes[i,1], 
                           mesh_completa$nodes[i,2])
  
  transport_y[i] = R_betay(mesh_completa$nodes[i,1], 
                           mesh_completa$nodes[i,2])
}

write.table(cbind(transport_x, transport_y), 
            file=paste0(foldername,"Beta.txt"), 
            row.names = FALSE, col.names = FALSE )
write.table(transport_x, file=paste0(foldername, "Beta_X.txt"), 
            row.names = FALSE, col.names = FALSE)
write.table(transport_y, file=paste0(foldername, "Beta_Y.txt"), 
            row.names = FALSE, col.names = FALSE)