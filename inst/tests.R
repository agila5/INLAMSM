# packages
options("rgdal_show_exportToProj4_warnings" = "none")
library("rgdal")
library("spdep")
library("INLA")
# library("INLAMSM")

# data
nc.sids <- readOGR(system.file("shapes/sids.shp", package = "spData"))
proj4string(nc.sids) <- CRS("+proj=longlat +ellps=clrk66")

#Compute adjacency matrix, as nb object 'adj' and sparse matrix 'W'
W <- igraph::as_adjacency_matrix(igraph::graph_from_adj_list(sf::st_touches(st_as_sf(nc.sids))))

# First time period
# Compute expected cases
r74 <- sum(nc.sids$SID74) / sum(nc.sids$BIR74)
nc.sids$EXP74 <- r74 * nc.sids$BIR74
# SMR
nc.sids$SMR74 <- nc.sids$SID74 / nc.sids$EXP74
# Proportion of non-white births
nc.sids$NWPROP74 <- nc.sids$NWBIR74 / nc.sids$BIR74

# Second time period
# Compute expected cases
r79 <- sum(nc.sids$SID79) / sum(nc.sids$BIR79)
nc.sids$EXP79 <- r79 * nc.sids$BIR79
# SMR
nc.sids$SMR79 <- nc.sids$SID79 / nc.sids$EXP79
# Proportion of non-white births
nc.sids$NWPROP79 <- nc.sids$NWBIR79 / nc.sids$BIR79

d <- data.frame(
  OBS = c(nc.sids$SID74, nc.sids$SID79),
  PERIOD = c(rep("74", 100), rep("79", 100)), 
  NWPROP = c(nc.sids$NWPROP74, nc.sids$NWPROP79),
  EXP = c(nc.sids$EXP74, nc.sids$EXP79)
)
# County-period index
d$idx <- 1:length(d$OBS)

# Number of variables (i.e., periods)
k <- 2

# Proper MCAR model
# Define range for the autocorrelation parameter
alpha.min <- 0
alpha.max <- 1
model_mcar <- inla.MCAR.model(k = k, W = W, alpha.min = alpha.min, alpha.max = alpha.max)

#Fit model
IMCAR <- inla(
  OBS ~ 0 + PERIOD + f(idx, model = model_mcar), 
  data = d, 
  E = EXP, 
  family = "poisson", 
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE), 
  verbose = TRUE
)

# Test my model
my_model <- inla.MCAR.model_v2(k = k, W = W, alpha.min = alpha.min, alpha.max = alpha.max)

# Fit my model
my_IMCAR <- inla(
  OBS ~ 0 + PERIOD + f(idx, model = my_model),
  data = d,
  E = EXP,
  family = "poisson",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE),
  verbose = TRUE
)

summary(IMCAR)
summary(my_IMCAR)







