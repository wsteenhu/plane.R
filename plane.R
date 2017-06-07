# Author R-code: W.A.A. de Steenhuijsen Piters, based on https://github.com/ElDeveloper/reference-plane
# Date: 6-6-2017
# E-mail: w.a.a.desteenhuijsenpiters@umcutrecht.nl

compute_coefficients <- function(xyz) {
  x = xyz[, 1]
  y = xyz[, 2]
  z = xyz[, 3]
  
  A = cbind(x, y)
  
  mod <- lm(z ~. , data=data.frame(A))
  abcd <- c(coef(mod)[[2]], coef(mod)[[3]], -1, coef(mod)[[1]])
  return(abcd) 
}

point_to_plane_distance <- function(abcd, point) {
  abc = abcd[1:3]
  d = abcd[4]
  dist = abs((abc %*% point) + d)
  dist/norm(matrix(abc), "F")
}

point_to_segment_distance <- function(abcd, point, xyz) {
  plane <- function(abcd, xy) {
    a = abcd[1];b = abcd[2];c = abcd[3];d = abcd[4]
    x = xy[1];y = xy[2]
    return (a*x + b*y + d)/(-1*c) }
  
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2)^2))
  
  a = abcd[1];b = abcd[2];c = abcd[3];d = abcd[4]
  p = point[1];q = point[2];r = point[3]
  l = ((d*-1) - p*a - b*q - c*r) / (a^2 + b^2 + c^2)
  extreme = c(p + l*a, q + l*b, r + l*c)
  
  for (i in 1:ncol(xyz)) {
    vector = xyz[,i]
    ranges = range(vector)
    if(extreme[i] < ranges[1]) { extreme[i] <- ranges[1]
    } else if(extreme[i] > ranges[2]) { extreme[i] <- ranges[2] } }
    extreme[length(extreme)] = plane(abcd, extreme[-length(extreme)])
    return(euc.dist(point, extreme))
}

distance_to_reference_plane <- function(ordination, meta, category1=NULL, column1=NULL, category2=NULL, column2=NULL) {
  if (is.null(category2)) { reference_ids <- which(meta[[column1]]==category1) }
  if (!is.null(category2)) { reference_ids <- which(meta[[column1]]==category1 & meta[[column2]]==category2) }
  reference <- ordination$vectors[reference_ids,1:3]
  abcd = compute_coefficients(reference)
  
  apply(ordination$vectors[,1:3], 1, function(x) point_to_segment_distance(abcd, x, reference))
}
#Note: ordination is the result of the pcoa-function in ape (see below), the order of samples in ordination$vectors should be the 
#same as the ordering in meta (samples in rows). Filtering of reference samples can be done based on two columns/categories.

#Example: distance_to_reference_plane(ord, meta, category1="Control", column1="Treatment")

pcoa_custom <- function(otu_table_raw, method="bray", log=T) {
  require(ape)
  if(log==T) { otu_table_raw <- log10(1 + otu_table_raw) }
  d <- vegdist(t(otu_table_raw), method=method)
  return(pcoa(d)) }
#otu_table_raw has taxa in rows, samples in columns.

#Example: ord <- pcoa_custom(raw_rare, log=T)
