# plane.R
Collection of functions implementing the calculation of a reference-plane (as described by Halfvarson, Nat Microbiol, 2017) in R instead of Python.

Typical usage would look like this:

source(file="/your/path/here/plane.R")  
ord <- pcoa_custom(otu_table)  
dtrp <- distance_to_reference_plane(ord, meta, column1 = "Treatment", category1 = "control")

Please find the original Python-implementation here: https://github.com/ElDeveloper/reference-plane.
