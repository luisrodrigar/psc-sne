## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(pscsne)
library(rotasym)
library(sphunif)
library(abind)
library(rgl)
library(mvtnorm)
num_cores_param <- 7

## -----------------------------------------------------------------------------
res_pscsne_11 <- psc_sne(X = sc_unif_mix_data, d = 1, rho_psc_list = rho_30_1,
                         num_iteration = 500, show_prog = TRUE,
                         colors = sc_unif_mix_colors,
                         parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
res_pscsne_12 <- psc_sne(X = sc_unif_mix_data, d = 2, rho_psc_list = rho_30_1,
                         num_iteration = 500, show_prog = TRUE,
                         colors = sc_unif_mix_colors,
                         parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
Y <- res_pscsne_12$best_Y
rgl::plot3d(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
            radius = 1, type = "s", col = "lightblue", alpha = 0.8, lit = FALSE)
rgl::points3d(Y, col = sc_unif_mix_colors)

## -----------------------------------------------------------------------------
# Visualize some features
n <- 100
g <- 5
p <- 20
x <- r_block(n = n, g = g, p = p)
pairs(x[, c(1:2, p + 1:2, (2 * p) + 1:2)],
      labels = c("Var 1", "Var 2", paste("Var", p + 1), paste("Var", p + 2),
                 paste("Var", (2 * p) + 1), paste("Var", (2 * p) + 2)))
# Standardize variables -- now the vectors of observations for each variable
# (the columns) live on \sqrt{n - 1} * S^{n - 1}!
x_sca <- scale(x)
# Make the features live on S^{n - 1}
x_sca <- x_sca / sqrt(n - 1)
# Transpose matrix (features become observations)
feat_data <- t(x_sca)
# Colors of the groups
cols <- rep(1:g, each = p)
# Set an array of dimension c(g * p, n, 1)
dim(feat_data) <- c(dim(feat_data), 1)
# Change the order of the data
indexes <- sample(1:(g * p))
feat_data <- feat_data[indexes, , , drop = FALSE]
cols <- cols[indexes]

## -----------------------------------------------------------------------------
rho_list <- rho_optim_bst(x = feat_data, perp_fixed = 30,
                          num_cores = num_cores_param)

## -----------------------------------------------------------------------------
res_pscsne_21 <- psc_sne(X = feat_data, d = 1, rho_psc_list = rho_list,
                      colors = cols, show_prog = TRUE, num_iteration = 775,
                      eta = 35, parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
psc_sne_res_22 <- psc_sne(feat_data, d = 2, rho_psc_list = rho_list,
                          colors = cols, num_iteration = 775, show_prog = TRUE,
                          eta = 35, parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
rgl::plot3d(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
            radius = 1, type = "s", col = "lightblue", alpha = 0.8, lit = FALSE)
rgl::points3d(psc_sne_res_22$best_Y, col = cols)

## -----------------------------------------------------------------------------
# Visualize some features
n <- 100
g <- 5
p <- 20
rho_values_neg <- rep(c(-0.9, 0.9), times = g)[1:g]
x <- r_block(n = n, g = g, p = p, rho = rho_values_neg)
pairs(x[, c(1:2, p + 1:2, (2 * p) + 1:2)],
      labels = c("Var 1", "Var 2", paste("Var", p + 1), paste("Var", p + 2),
                 paste("Var", (2 * p) + 1), paste("Var", (2 * p) + 2)))
# Standardize variables -- now the vectors of observations for each variable
# (the columns) live on \sqrt{n - 1} * S^{n - 1}!
x_sca <- scale(x)
# Make the features live on S^{n - 1}
x_sca <- x_sca / sqrt(n - 1)
# Transpose matrix (features become observations)
feat_data_2 <- t(x_sca)
# Colors of the groups
cols <- rep(1:g, each = p)
# Set an array of dimension c(g * p, n, 1)
dim(feat_data_2) <- c(dim(feat_data_2), 1)
# Change the order of the data
indexes <- sample(1:(g * p))
feat_data_2 <- feat_data_2[indexes, , , drop = FALSE]
cols <- cols[indexes]

## -----------------------------------------------------------------------------
rho_list_2 <- rho_optim_bst(x = feat_data_2, perp_fixed = 30,
                            num_cores = num_cores_param)

## -----------------------------------------------------------------------------
res_pscsne_21_neg <- psc_sne(X = feat_data_2, d = 1, rho_psc_list = rho_list_2,
                      colors = cols, show_prog = TRUE, num_iteration = 775,
                      eta = 35, parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
psc_sne_res_22_neg <- psc_sne(feat_data_2, d = 2, rho_psc_list = rho_list_2,
                              colors = cols, num_iteration = 775,
                              show_prog = TRUE, eta = 35,
                              parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
rgl::plot3d(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
            radius = 1, type = "s", col = "lightblue", alpha = 0.8, lit = FALSE)
rgl::points3d(psc_sne_res_22_neg$best_Y, col = cols)

## -----------------------------------------------------------------------------
cols <- ifelse(gen_data <= n / 2, 2, 3)
psc_sne_res_31 <- psc_sne(X = x_s1_100, d = 1, rho_psc_list = rho_30_s1_100,
                          colors = cols, show_prog = TRUE, num_iteration = 775,
                          eta = 25, parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
psc_sne_res_32 <- psc_sne(X = x_s1_100, d = 2, rho_psc_list = rho_30_s1_100,
                          colors = cols, show_prog = TRUE, num_iteration = 775,
                          eta = 25, parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
rgl::plot3d(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
             radius = 1, type = "s", col = "lightblue", lit = FALSE)
rgl::points3d(psc_sne_res_32$best_Y, col = cols)

## -----------------------------------------------------------------------------
psc_sne_res_41 <- psc_sne(X = x_s1_100_path, d = 1,
                          rho_psc_list = rho_30_s1_100_path, colors = cols,
                          show_prog = TRUE, num_iteration = 775, eta = 25,
                          parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
psc_sne_res_42 <- psc_sne(X = x_s1_100_path, d = 2,
                          rho_psc_list = rho_30_s1_100_path, colors = cols,
                          show_prog = TRUE, num_iteration = 775, eta = 25,
                          parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
rgl::plot3d(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
             radius = 1, type = "s", col = "lightblue", lit = FALSE)
rgl::points3d(psc_sne_res_42$best_Y, col = cols)

## ----cache = TRUE-------------------------------------------------------------
rho_30_s2_100 <- rho_optim_bst(x = x_s2_100, perp_fixed = 30,
                               num_cores = num_cores_param)

## -----------------------------------------------------------------------------
psc_sne_res_51 <- psc_sne(X = x_s2_100, d = 1, rho_psc_list = rho_30_s2_100,
                          colors = cols, show_prog = TRUE, num_iteration = 775,
                          eta = 25, parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
psc_sne_res_52 <- psc_sne(X = x_s2_100, d = 2, rho_psc_list = rho_30_s2_100,
                          colors = cols, show_prog = TRUE, num_iteration = 775,
                          eta = 25, parallel_cores = num_cores_param)

## -----------------------------------------------------------------------------
rgl::plot3d(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
             radius = 1, type = "s", col = "lightblue", lit = FALSE)
rgl::points3d(psc_sne_res_52$best_Y, col = cols)

