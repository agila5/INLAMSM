# packages
library("spdep")
library("Matrix")
library("INLA")
library("INLAMSM")

#Load SIDS data
nc.sids <- st_read(system.file("shapes/sids.shp", package = "spData"))
st_crs(nc.sids) <- "+proj=longlat +ellps=clrk66"

#Compute adjacency matrix, as nb object 'adj' and sparse matrix 'W'
adj <- poly2nb(nc.sids)
W <- igraph::as_adjacency_matrix(igraph::graph_from_adj_list(st_touches(st_as_sf(nc.sids))))

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

# define effect
my_effect <- MCAR.model_a1_ak_Lambda(W = W, alpha.min = 0, alpha.max = 1, k = k)

# estimate model
my_model <- inla(
  OBS ~ 0 + PERIOD + f(idx, model = my_effect),
  data = d,
  E = EXP,
  family = "poisson",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)
summary(my_model)

# define inv_eps_min for the case 3
fnMatSqrtInverse = function(mA) {
  ei = eigen(mA)
  d = ei$values
  d = (d + abs(d)) / 2
  d2 = 1/sqrt(d)
  d2[d == 0] = 0
  return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}

D_m05 <- fnMatSqrtInverse(Diagonal(nrow(W), rowSums(W)))
inv_eps_min <- 1 / min(eigen(D_m05 %*% W %*% D_m05)$values)

# define effect for model - case 3
my_effect <- MCAR.model_case3(W = W, inv_eps_min = inv_eps_min)

# estimate model
my_model <- inla(
  OBS ~ 0 + PERIOD + f(idx, model = my_effect), 
  data = d,
  E = EXP, 
  family = "poisson", 
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE), 
  verbose = TRUE
)
summary(my_model)
