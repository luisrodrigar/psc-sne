## ----setup, include = FALSE---------------------------------------------------
def.chunk.hook <- knitr::knit_hooks$get("chunk")
knitr::knit_hooks$set(chunk = function(x, options) {
  x <- def.chunk.hook(x, options)
  ifelse(options$size != "normalsize",
         paste0("\n \\", options$size, "\n\n", x, "\n\n \\normalsize"), x)
})
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
stopifnot(packageVersion("pscsne") <= "0.0.1.900004")

## ----cache = TRUE-------------------------------------------------------------
sc_unif_mix <- function(n, p, w_sc, w_unif, kappa = 50) {
  if (w_sc + w_unif != 1) {
    stop("w_sc and w_unif must sum 1")
  }
  n1 <- rbinom(1, n, w_sc)
  n2 <- n - n1
  r_1 <- sphunif::r_alt(n = n1, p = p, alt = "SC", kappa = kappa)
  r_2 <- sphunif::r_unif_sph(n = n2, p = p)
  data <- abind(r_1, r_2, along = 1)
  # Change the order of the data
  indexes <- sample(1:n)
  data <- data[indexes, , , drop = FALSE]
  cols <- c(rep(1, times = n1), rep(2, times = n2))
  cols <- cols[indexes]
  return(list("data" = data, "colors" = cols))
}
n <- 800
p <- 101
w_sc <- 0.5
w_unif <- 0.5
kappa <- 1000
set.seed(42)
sc_unif_mix_res <- sc_unif_mix(n = n, p = p, w_sc = w_sc, w_unif = w_unif,
                               kappa = kappa)
sc_unif_mix_data <- sc_unif_mix_res$data
sc_unif_mix_colors <- sc_unif_mix_res$colors

## ----cache = TRUE-------------------------------------------------------------
rho_30_1 <- rho_optim_bst(x = sc_unif_mix_data, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----cache = TRUE, size = "scriptsize", fig.align = 'center'------------------
res_pscsne_11 <- psc_sne(X = sc_unif_mix_data, d = 1, rho_psc_list = rho_30_1,
                         show_prog = TRUE, colors = sc_unif_mix_colors,
                         parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_11$best_Y
plot(Y[, 1], Y[, 2], col = sc_unif_mix_colors, xlim = c(-1, 1),
     ylim = c(-1, 1))
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
res_pscsne_12 <- psc_sne(X = sc_unif_mix_data, d = 2, rho_psc_list = rho_30_1,
                         show_prog = TRUE, colors = sc_unif_mix_colors,
                         parallel_cores = num_cores_param, eta = 100)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
seq_rad <- seq(-pi, pi, by = pi / 30)
meridian <- do.call(rbind, lapply(seq_rad, function(i) c(0, i)))
equator <- do.call(rbind, lapply(seq_rad, function(i) c(i, pi/2)))
Y <- res_pscsne_12$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = sc_unif_mix_colors, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c('+', '-')[ifelse(sign(Y[, 2]) == 1, 1, 2)]
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----cache = TRUE-------------------------------------------------------------
# Visualize some features
n <- 200
g <- 5
p <- 20
set.seed(42)
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
set.seed(42)
indexes <- sample(1:(g * p))
feat_data <- feat_data[indexes, , , drop = FALSE]
cols <- cols[indexes]

## ----cache = TRUE-------------------------------------------------------------
rho_list <- rho_optim_bst(x = feat_data, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----cache = TRUE, size = "scriptsize", fig.align = 'center'------------------
res_pscsne_21 <- psc_sne(X = feat_data, d = 1, rho_psc_list = rho_list,
                         colors = cols, show_prog = TRUE,
                         eta = 35, parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_21$best_Y
plot(Y[, 1], Y[, 2], col = cols, xlim = c(-1, 1), ylim = c(-1, 1))
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_22 <- psc_sne(feat_data, d = 2, rho_psc_list = rho_list,
                          colors = cols, show_prog = TRUE,
                          eta = 35, parallel_cores = num_cores_param)

## ----fig.asp = 1--------------------------------------------------------------
Y <- psc_sne_res_22$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c('+', '-')[ifelse(sign(Y[, 2]) == 1, 1, 2)]
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----cache = TRUE-------------------------------------------------------------
# Visualize some features
n <- 100
g <- 5
p <- 20
rho_values_neg <- rep(c(-0.9, 0.9), times = g)[1:g]
set.seed(42)
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
set.seed(42)
indexes <- sample(1:(g * p))
feat_data_2 <- feat_data_2[indexes, , , drop = FALSE]
cols <- cols[indexes]

## ----cache = TRUE-------------------------------------------------------------
rho_list_2 <- rho_optim_bst(x = feat_data_2, perp_fixed = 30,
                            num_cores = num_cores_param)

## ----cache = TRUE, size = "scriptsize", fig.align = 'center'------------------
res_pscsne_21_neg <- psc_sne(X = feat_data_2, d = 1, rho_psc_list = rho_list_2,
                             colors = cols, show_prog = TRUE, eta = 35,
                             parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_21_neg$best_Y
plot(Y[, 1], Y[, 2], col = cols, xlim = c(-1, 1), ylim = c(-1, 1))
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_22_neg <- psc_sne(feat_data_2, d = 2, rho_psc_list = rho_list_2,
                              colors = cols, show_prog = TRUE, eta = 20,
                              parallel_cores = num_cores_param)

## ----fig.asp = 1--------------------------------------------------------------
Y <- psc_sne_res_22_neg$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c('+', '-')[ifelse(sign(Y[, 2]) == 1, 1, 2)]
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----cache = TRUE-------------------------------------------------------------
r1 <- 5
r2 <- 95
p <- 1
n <- 100
kappa <- 3
set.seed(42)
gen_data <- rbinom(n = n, size = n, prob = 0.5)
x_vMF_5 <- sapply(seq_len(r1), function(k1) {
  data_vMF <- lapply(seq_len(n), function(i) {
    if (gen_data[i] <= n / 2) {
      set.seed(42)
      rotasym::r_vMF(n = 1, mu = c(0, 1), kappa = kappa)
    } else {
      set.seed(42)
      rotasym::r_vMF(n = 1, mu = c(0, -1), kappa = kappa)
    }
  })
  do.call(rbind, data_vMF)
}, simplify = "array")
set.seed(42)
x_uniform_95 <- sphunif::r_unif_sph(n = n, p = p + 1, M = r2)
x_s1_100 <- abind(x_vMF_5, x_uniform_95, along = 3)

## ----cache = TRUE-------------------------------------------------------------
rho_30_s1_100 <- rho_optim_bst(x = x_s1_100, perp_fixed = 30,
                               num_cores = num_cores_param)

## ----cache = TRUE, size = "scriptsize", fig.align = 'center'------------------
cols <- ifelse(gen_data <= n / 2, 2, 3)
psc_sne_res_31 <- psc_sne(X = x_s1_100, d = 1, rho_psc_list = rho_30_s1_100,
                          colors = cols, show_prog = TRUE, eta = 25,
                          parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- psc_sne_res_31$best_Y
plot(Y[, 1], Y[, 2], col = cols, xlim = c(-1, 1), ylim = c(-1, 1))
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_32 <- psc_sne(X = x_s1_100, d = 2, rho_psc_list = rho_30_s1_100,
                          colors = cols, show_prog = TRUE, eta = 25,
                          parallel_cores = num_cores_param)

## ----fig.asp = 1--------------------------------------------------------------
Y <- psc_sne_res_32$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c('+', '-')[ifelse(sign(Y[, 2]) == 1, 1, 2)]
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----cache = TRUE-------------------------------------------------------------
n <- 100
r <- 100
set.seed(42)
x_s1_100_path <- r_path_s1r(n = n, r = r)
cols <- rainbow(n, alpha = 1)
# Change the order of the data
set.seed(42)
indexes <- sample(1:n)
x_s1_100_path <- x_s1_100_path[indexes, , ]
cols <- cols[indexes]

## ----cache = TRUE-------------------------------------------------------------
rho_30_s1_100_path <- rho_optim_bst(x = x_s1_100_path, perp_fixed = 30,
                                    num_cores = num_cores_param)

## ----cache = TRUE, size = "scriptsize", fig.align = 'center'------------------
psc_sne_res_41 <- psc_sne(X = x_s1_100_path, d = 1,
                          rho_psc_list = rho_30_s1_100_path, colors = cols,
                          show_prog = TRUE, eta = 25,
                          parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- psc_sne_res_41$best_Y
plot(Y[, 1], Y[, 2], col = cols, xlim = c(-1, 1), ylim = c(-1, 1))
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_42 <- psc_sne(X = x_s1_100_path, d = 2,
                          rho_psc_list = rho_30_s1_100_path, colors = cols,
                          show_prog = TRUE, eta = 25,
                          parallel_cores = num_cores_param)

## ----fig.asp = 1--------------------------------------------------------------
Y <- psc_sne_res_42$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c('+', '-')[ifelse(sign(Y[, 2]) == 1, 1, 2)]
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----cache = TRUE-------------------------------------------------------------
n <- 100
r1 <- 1
r2 <- 2
r3 <- 2
r4 <- 95
# Calculate the north pole of the sphere
north_pole <- cbind(c(0, 0, 1))
# small circle in the equator (1 sample)
set.seed(42)
samp_1 <- r_path_s2r(n = n, r = r1, sigma = 0.08, Theta = north_pole)
# small circle rotated 1 and 2 (2 samples)
set.seed(42)
samp_2_3 <- r_path_s2r(n = n, r = r2, sigma = 0.1)
# spherical spiral 1 and 2 (2 samples)
set.seed(42)
samp_4_5 <- r_path_s2r(n = n, r = r3, c = 3, spiral = TRUE, sigma = 0.01)
# Data following an uniform distribution on the sphere (95 samples)
set.seed(42)
samp_95 <- r_unif_sph(n = n, p = 3, M = r4)
# Join the data by the third dimension
x_s2_100 <- abind(samp_1, samp_2_3, samp_4_5, samp_95, along = 3)
# Create rainbow colors, alpha is the level of opacity
cols <- rainbow(n, alpha = 1)
# Change the order of the data
set.seed(42)
indexes <- sample(1:n)
x_s2_100 <- x_s2_100[indexes, , ]
cols <- cols[indexes]

## ----cache = TRUE-------------------------------------------------------------
rho_30_s2_100 <- rho_optim_bst(x = x_s2_100, perp_fixed = 30,
                               num_cores = num_cores_param)

## ----cache = TRUE, size = "scriptsize", fig.align = 'center'------------------
psc_sne_res_51 <- psc_sne(X = x_s2_100, d = 1, rho_psc_list = rho_30_s2_100,
                          colors = cols, show_prog = TRUE, eta = 25,
                          parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- psc_sne_res_51$best_Y
plot(Y[, 1], Y[, 2], col = cols, xlim = c(-1, 1), ylim = c(-1, 1))
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_52 <- psc_sne(X = x_s2_100, d = 2, rho_psc_list = rho_30_s2_100,
                          colors = cols, show_prog = TRUE, eta = 25,
                          parallel_cores = num_cores_param)

## ----fig.asp = 1--------------------------------------------------------------
Y <- psc_sne_res_52$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c('+', '-')[ifelse(sign(Y[, 2]) == 1, 1, 2)]
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
set.seed(42)
psc_sne_res_52 <- psc_sne(X = x_s2_100, d = 2, rho_psc_list = rho_30_s2_100,
                          colors = cols, show_prog = TRUE, eta = 30,
                          parallel_cores = num_cores_param, init = "random")

## ----fig.asp = 1--------------------------------------------------------------
Y <- psc_sne_res_52$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c('+', '-')[ifelse(sign(Y[, 2]) == 1, 1, 2)]
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

