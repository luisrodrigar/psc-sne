
# Data
small_rna <- read.table("small_rna.csv", header = FALSE, sep = ",")

# Labels based on 8:9 angles, later ignored
cluster1 <- small_rna[small_rna[, 9] > -60, ]
cluster2 <- small_rna[(small_rna[, 8] > 0) & (small_rna[, 9] < -60), ]
cluster5 <- small_rna[(small_rna[, 8] < 0) & (small_rna[, 9] < -60), ]
clusters_ordered <- rbind(cluster1, cluster2, cluster5)
small_rna_angles <- clusters_ordered[, 1:7]
small_rna_torsion <- clusters_ordered[, 8:9] # Ignored later

# Clustering variables
col1 <- rep(1, 59)
col2 <- rep(2, 88)
col3 <- rep(3, 43)
col <- c(col1, col2, col3)
plot(clusters_ordered[, 8:9], col = col) # Plot pseudotorsion angles

# Dataset with angles 1:7
small_rna_angles <- as.matrix(small_rna_angles)
small_rna_angles <- (small_rna_angles / 180) * pi # Degrees to radians

# Torsion angles
small_rna_torsion <- (small_rna_torsion / 180) * pi # Degrees to radians

# View data on (S^1)^7 = T^7
pairs(small_rna_angles, col = col)

# Construct data frame and save
smallrna <- data.frame("angles" = I(small_rna_angles),
                       "torsion" = I(small_rna_torsion),
                       "clusters" = col)
save(list = "smallrna", file = "smallrna.rda", compress = "bzip2")
