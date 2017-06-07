# plane.R
Simple collection of functions implementing the calculation of a reference-plane (as described by Halfvarson, Nat Microbiol, 2017) in R instead of python.

Typical usage would look like this:

source(file="/your/path/here/plane.R")  
ord <- pcoa_custom(otu_table)
dtrp <- distance_to_reference_plane(ord, meta, column1 = "Treatment", category1 = "control", column2 = "Timepoint", category2 = "baseline")
