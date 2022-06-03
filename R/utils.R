clusterFactory <- function(num_cores, outfile = "") {
  if(tolower(.Platform$OS.type) != "windows"){
    cl <- makeCluster(spec=num_cores, type="FORK", outfile=outfile)  
  } else
    cl <- makeCluster(spec=num_cores, outfile=outfile)
}
