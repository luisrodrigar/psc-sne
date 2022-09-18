## ----setup, include = FALSE---------------------------------------------------
def_chunk_hook <- knitr::knit_hooks$get("chunk")
knitr::knit_hooks$set(chunk = function(x, options) {
  x <- def_chunk_hook(x, options)
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
num_cores_param <- 2
stopifnot(packageVersion("pscsne") >= "0.0.1.900004")

## ----cache = TRUE-------------------------------------------------------------
MvMF_mix <- function(n, n1, n2, n3, mu_mat, kappa = 50) {
  if ((n1 + n2 + n3) != n) {
    stop("(n1 + n2 + n3) must sum 1")
  }
  r_1 <- rotasym::r_vMF(n1, mu = mu_mat[, 1], kappa = kappa)
  r_2 <- rotasym::r_vMF(n2, mu = mu_mat[, 2], kappa = kappa)
  r_3 <- rotasym::r_vMF(n3, mu = mu_mat[, 3], kappa = kappa)
  r_mvmf <- rbind(r_1, r_2, r_3)
  colors <- c(rep(1, n1), rep(2, n2), rep(3, n3))
  return(list("data" = array(r_mvmf, dim = c(dim(r_mvmf), 1)), "colors" = colors))
}
n <- 500
n1 <- ceiling(n / 3)
n2 <- ceiling(n / 3)
n3 <- floor(n / 3)
p <- 2
kappa <- 100
mu_mat <- cbind(c(0, 1, 0), c(0, 0, 1), c(0, 0, -1))
set.seed(42)
mvmf_mix_res <- MvMF_mix(n = n, n1 = n1, n2 = n2, n3 = n3,
                         mu_mat = mu_mat, kappa = kappa)
mvmf_mix_data <- mvmf_mix_res$data
mvmf_mix_colors <- mvmf_mix_res$colors
# shuffle the sample
set.seed(42)
indexes <- sample(1:n)
mvmf_mix_data <- mvmf_mix_data[indexes, , , drop = FALSE]
mvmf_mix_colors <- mvmf_mix_colors[indexes]

## ----fig.asp = 1, fig.align='center'------------------------------------------
seq_rad <- seq(-pi, pi, by = pi / 30)
meridian <- do.call(rbind, lapply(seq_rad, function(i) c(0, i)))
equator <- do.call(rbind, lapply(seq_rad, function(i) c(i, pi / 2)))
sd3 <- scatterplot3d::scatterplot3d(
  mvmf_mix_data[, , 1], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = mvmf_mix_colors, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(mvmf_mix_data[, 2, 1]) == 1, 1, 2)],
  mar = c(0,0,0,0), grid = FALSE
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----cache = TRUE-------------------------------------------------------------
rho_30_1 <- rho_optim_bst(x = mvmf_mix_data, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----cache = TRUE, size = "scriptsize", fig.align = 'center'------------------
res_pscsne_11 <- psc_sne(X = mvmf_mix_data, d = 1, rho_psc_list = rho_30_1,
                         show_prog = TRUE, colors = mvmf_mix_colors, eta = 50,
                         parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_11$best_Y
plot(Y[, 1], Y[, 2], col = mvmf_mix_colors, xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "", ylab = "", axes = FALSE)
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
res_pscsne_12 <- psc_sne(X = mvmf_mix_data, d = 2, rho_psc_list = rho_30_1,
                         show_prog = TRUE, colors = mvmf_mix_colors,
                         parallel_cores = num_cores_param, eta = 100)

## ----fig.asp = 1, fig.align='center'------------------------------------------
Y <- res_pscsne_12$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = mvmf_mix_colors, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(Y[, 2]) == 1, 1, 2)], grid = FALSE,
  mar = c(0, 0, 0, 0)
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## -----------------------------------------------------------------------------
n <- 500
n1 <- ceiling(n / 4)
n2 <- ceiling(n / 4)
n3 <- floor(n / 4)
n4 <- floor(n / 4)
kappa <- 50
north_pole <- cbind(c(0, 0, 1))
set.seed(42)
spiral_sph <- r_path_s2r(n = n4, r = 1, c = 3, spiral = TRUE, sigma = 0.01,
                         Theta = north_pole)
mu_mat <- cbind(c(0, 1, 0), c(0, 0, 1), c(0, 0, -1))
set.seed(42)

mvmf_mix_res <- MvMF_mix(n = n1 + n2 + n3, n1 = n1, n2 = n2, n3 = n3,
                         mu_mat = mu_mat, kappa = kappa)
mvmf_mix_data <- mvmf_mix_res$data
mvmf_mix_colors <- mvmf_mix_res$colors
mvmf_spiral_mix_data <- abind(mvmf_mix_data, spiral_sph, along = 1)
mvmf_spiral_mix_colors <- c(mvmf_mix_colors, rep(4, times = n4))
# shuffle the sample
set.seed(42)
indexes <- sample(1:n)
mvmf_spiral_mix_data <- mvmf_spiral_mix_data[indexes, , , drop = FALSE]
mvmf_spiral_mix_colors <- mvmf_spiral_mix_colors[indexes]

## ----fig.asp = 1, fig.align='center'------------------------------------------
sd3 <- scatterplot3d::scatterplot3d(
  mvmf_spiral_mix_data[, , 1], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = mvmf_spiral_mix_colors, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(mvmf_spiral_mix_data[, 2, 1]) == 1, 1, 2)], 
  mar = c(0,0,0,0), grid = FALSE
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----cache = TRUE-------------------------------------------------------------
rho_30_2 <- rho_optim_bst(x = mvmf_spiral_mix_data, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----cache = TRUE, size = "scriptsize", fig.align = 'center'------------------
res_pscsne_21 <- psc_sne(X = mvmf_spiral_mix_data, d = 1,
                         rho_psc_list = rho_30_2, eta = 50,
                         colors = mvmf_spiral_mix_colors, show_prog = TRUE,
                         parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_21$best_Y
plot(Y[, 1], Y[, 2], col = mvmf_spiral_mix_colors, xlim = c(-1, 1),
     ylim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_22 <- psc_sne(mvmf_spiral_mix_data, d = 2, rho_psc_list = rho_30_2,
                          colors = mvmf_spiral_mix_colors, show_prog = TRUE,
                          parallel_cores = num_cores_param, eta = 100)

## ----fig.asp = 1, fig.align='center'------------------------------------------
Y <- psc_sne_res_22$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = mvmf_spiral_mix_colors, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(Y[, 2]) == 1, 1, 2)], grid = FALSE,
  mar = c(0, 0, 0, 0)
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## -----------------------------------------------------------------------------
n <- 500
n1 <- ceiling(n / 3)
n2 <- ceiling(n / 3)
n3 <- floor(n / 3)
set.seed(42)
r_1 <- r_path_s2r(n = n1, r = 1, Theta = rbind(0, 0, 1), sigma = 0.05,
                  t = -0.75)
set.seed(42)
r_2 <- r_path_s2r(n = n2, r = 1, Theta = rbind(0, 0, 1), sigma = 0.05, t = 0)
set.seed(42)
r_3 <- r_path_s2r(n = n3, r = 1, Theta = rbind(0, 0, 1), sigma = 0.05,
                  t = 0.75)
mix_sc <- abind(r_1, r_2, r_3, along = 1)
mix_sc_colors <- c(rep(1, n1), rep(2, n2), rep(3, n3))
set.seed(42)
indexes <- sample(1:n)
mix_sc <- mix_sc[indexes, , , drop = FALSE]
mix_sc_colors <- mix_sc_colors[indexes]

## ----fig.asp = 1, fig.align='center'------------------------------------------
sd3 <- scatterplot3d::scatterplot3d(
  mix_sc[, , 1], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = mix_sc_colors, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(mix_sc[, 2, 1]) == 1, 1, 2)],
  mar = c(0,0,0,0), grid = FALSE
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----cache = TRUE-------------------------------------------------------------
rho_30_3 <- rho_optim_bst(x = mix_sc, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----size = "scriptsize", fig.align = 'center'--------------------------------
res_pscsne_31 <- psc_sne(X = mix_sc, d = 1, rho_psc_list = rho_30_3, eta = 50,
                         colors = mix_sc_colors, show_prog = TRUE, 
                         perplexity = 30, parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_31$best_Y
plot(Y[, 1], Y[, 2], col = mix_sc_colors, xlim = c(-1, 1), ylim = c(-1, 1),
     axes = FALSE, xlab = "", ylab = "")
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_32 <- psc_sne(mix_sc, d = 2, rho_psc_list = rho_30_3, eta = 50,
                          colors = mix_sc_colors, show_prog = TRUE,
                          parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align='center'------------------------------------------
Y <- psc_sne_res_32$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = mix_sc_colors, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(Y[, 2]) == 1, 1, 2)], grid = FALSE,
  mar = c(0, 0, 0, 0)
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## -----------------------------------------------------------------------------
n <- 500
n1 <- ceiling(n / 4)
n2 <- ceiling(n / 4)
n3 <- floor(n / 4)
n4 <- floor(n / 4)
mix_3_mvMF <- function(n, n1, n2, n3, n4, mu_mat, kappa_vec) {
  if ((n1 + n2 + n3 + n4) != n) {
    stop("(n1 + n2 + n3) must sum 1")
  }
  r_1 <- rotasym::r_vMF(n = n1, mu = mu_mat[, 1], kappa = kappa_vec[1])
  r_2 <- rotasym::r_vMF(n = n2, mu = mu_mat[, 2], kappa = kappa_vec[2])
  r_3 <- rotasym::r_vMF(n = n3, mu = mu_mat[, 3], kappa = kappa_vec[3])
  r_4 <- rotasym::r_vMF(n = n4, mu = mu_mat[, 4], kappa = kappa_vec[4])
  data <- rbind(r_1, r_2, r_3, r_4)
  return(array(data, dim = c(dim(data), 1)))
}
mu_mat_1 <- cbind(drop(DirStats::to_cir(0)),
                drop(DirStats::to_cir(-pi / 2)),
                drop(DirStats::to_cir(pi / 2)),
                drop(DirStats::to_cir(3 * pi / 2.5)))
set.seed(42)
sph_1 <- mix_3_mvMF(n, n1, n2, n3, n4, mu_mat = mu_mat_1,
                    kappa_vec = c(30, 25, 10, 50))
mu_mat_2 <- cbind(drop(DirStats::to_cir(0)),
                drop(DirStats::to_cir(-pi)),
                drop(DirStats::to_cir(0)),
                drop(DirStats::to_cir(-3 * pi / 2.5)))
set.seed(44)
sph_2 <- mix_3_mvMF(n, n1, n2, n3, n4, mu_mat = mu_mat_2,
                    kappa_vec = c(25, 32, 12, 100))
data <- abind(sph_1, sph_2, along = 3)
set.seed(42)
indexes <- sample(1:n)
mix_3_mvMF_data <- data[indexes, , , drop = FALSE]
mix_3_mvMF_cols <- c(rep(1, n1), rep(2, n2), rep(3, n3), rep(4, n4))
mix_3_mvMF_cols <- mix_3_mvMF_cols[indexes]

## ----fig.asp = 1--------------------------------------------------------------
mix_3_mvMF_torus <- sdetorus::toPiInt(cbind(
  DirStats::to_rad(mix_3_mvMF_data[, , 1]),
  DirStats::to_rad(mix_3_mvMF_data[, , 2])))

plot(mix_3_mvMF_torus, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
     col = rep(mix_3_mvMF_cols, 2), pch = 16,
     xlab = "", ylab = "")
sdetorus::torusAxis()

## ----cache = TRUE-------------------------------------------------------------
rho_30_4 <- rho_optim_bst(x = mix_3_mvMF_data, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----size = "scriptsize", fig.align = 'center'--------------------------------
res_pscsne_41 <- psc_sne(X = mix_3_mvMF_data, d = 1, rho_psc_list = rho_30_4,
                         eta = 50, colors = mix_3_mvMF_cols, show_prog = TRUE,
                         parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_41$best_Y
plot(Y[, 1], Y[, 2], col = mix_3_mvMF_cols, xlim = c(-1, 1), ylim = c(-1, 1),
     axes = FALSE, xlab = "", ylab = "")
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_42 <- psc_sne(mix_3_mvMF_data, d = 2, rho_psc_list = rho_30_4,
                          eta = 25, colors = mix_3_mvMF_cols,
                          show_prog = TRUE, parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align='center'------------------------------------------
Y <- psc_sne_res_42$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = mix_3_mvMF_cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(Y[, 2]) == 1, 1, 2)], grid = FALSE,
  mar = c(0, 0, 0, 0)
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
n <- 500
n1 <- ceiling(n / 3)
n2 <- ceiling(n / 3)
n3 <- floor(n / 3)
mean1 <- c(1, 1.25)
mean2 <- c(-1.8, -0.5)
mean3 <- c(2.5, 2.75)
sigma1 <- matrix(c(0.2, -0.12, -0.12, 0.2), byrow = TRUE, nrow = 2)
sigma2 <- toeplitz(c(0.12, 0.1))
sigma3 <- toeplitz(c(0.2, -0.15))
mix_3_mvnorm <- function(n, n1, n2, n3, mean_mat, sigma1, sigma2, sigma3,
                       kappa = 10) {
  if ((n1 + n2 + n3) != n) {
    stop("(n1 + n2 + n3) must sum 1")
  }
  r_1 <- mvtnorm::rmvnorm(n = n1, mean = mean_mat[, 1], sigma = sigma1)
  r_2 <- mvtnorm::rmvnorm(n = n2, mean = mean_mat[, 2], sigma = sigma2)
  r_3 <- mvtnorm::rmvnorm(n = n3, mean = mean_mat[, 3], sigma = sigma3)
  data <- rbind(r_1, r_2, r_3)
  return(data)
}
set.seed(5)
x_5_bvnorm <- mix_3_mvnorm(n = n, n1 = n1, n2 = n2, n3 = n3,
                           mean_mat = cbind(mean1, mean2, mean3),
                           sigma1 = sigma1, sigma2 = sigma2, sigma3 = sigma3)
x_5_torus <- sdetorus::toPiInt(x_5_bvnorm)

samp_5 <- abind(DirStats::to_cir(x_5_torus[, 1]),
                DirStats::to_cir(x_5_torus[, 2]), along = 3)
samp_5_cols <- c(rep(1, n1), rep(2, n2), rep(3, n3))

plot(x_5_torus, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
     col = samp_5_cols, pch = 16,
     xlab = "x", ylab = "y")
sdetorus::torusAxis()

## -----------------------------------------------------------------------------
set.seed(42)
indexes <- sample(1:n)
samp_5 <- samp_5[indexes, , , drop = FALSE]
samp_5_cols <- samp_5_cols[indexes]

## ----cache = TRUE-------------------------------------------------------------
rho_30_5 <- rho_optim_bst(x = samp_5, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----size = "scriptsize", fig.align = 'center'--------------------------------
res_pscsne_51 <- psc_sne(X = samp_5, d = 1, rho_psc_list = rho_30_5,
                         eta = 100, colors = samp_5_cols, show_prog = TRUE,
                         parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_51$best_Y
plot(Y[, 1], Y[, 2], col = samp_5_cols, xlim = c(-1, 1), ylim = c(-1, 1),
     axes = FALSE, xlab = "", ylab = "")
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_52 <- psc_sne(samp_5, d = 2, rho_psc_list = rho_30_5,
                          eta = 50, colors = samp_5_cols,
                          show_prog = TRUE, parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align='center'------------------------------------------
Y <- psc_sne_res_52$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = samp_5_cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(Y[, 2]) == 1, 1, 2)], grid = FALSE,
  mar = c(0, 0, 0, 0)
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----fig.asp = 1--------------------------------------------------------------
n <- 500
set.seed(3)
samp_1 <- r_path_s1r(n = n, r = 2, k = c(1, 2), angles = TRUE)
plot(samp_1, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(n, alpha = 1),
     axes = FALSE, xlab = "", ylab = "", pch = 16)
sdetorus::torusAxis()

## -----------------------------------------------------------------------------
x_1 <- DirStats::to_cir(samp_1[, 1])
x_2 <- DirStats::to_cir(samp_1[, 2])
x_path_1 <- abind(x_1, x_2, along = 3)
x_path_1_cols <- rainbow(n, alpha = 1)
set.seed(42)
indexes <- sample(1:n)
x_path_1 <- x_path_1[indexes, , , drop = FALSE]
x_path_1_cols <- x_path_1_cols[indexes]

## -----------------------------------------------------------------------------
rho_30_6 <- rho_optim_bst(x = x_path_1, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----size = "scriptsize", fig.align = 'center'--------------------------------
set.seed(42)
res_pscsne_61 <- psc_sne(X = x_path_1, d = 1, rho_psc_list = rho_30_6,
                         eta = 25, colors = x_path_1_cols, show_prog = TRUE,
                         parallel_cores = num_cores_param, init = "random")

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_61$best_Y
plot(Y[, 1], Y[, 2], col = x_path_1_cols, xlim = c(-1, 1), ylim = c(-1, 1),
     axes = FALSE, xlab = "", ylab = "")
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_62 <- psc_sne(x_path_1, d = 2, rho_psc_list = rho_30_6,
                          eta = 15, colors = x_path_1_cols,
                          show_prog = TRUE, parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align='center'------------------------------------------
Y <- psc_sne_res_62$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = x_path_1_cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(Y[, 2]) == 1, 1, 2)], grid = FALSE,
  mar = c(0, 0, 0, 0)
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## -----------------------------------------------------------------------------
n <- 500
n1 <- ceiling(n / 4)
n2 <- ceiling(n / 4)
n3 <- floor(n / 4)
n4 <- floor(n / 4)
mu_mat_1 <- cbind(drop(DirStats::to_cir(3 * pi / 2.5)),
                drop(DirStats::to_cir(-pi / 2)),
                drop(DirStats::to_cir(pi / 2)),
                drop(DirStats::to_cir(0)))
kappa_vec1 <- c(50, 100, 30, 15)
set.seed(3)
sph_1 <- mix_3_mvMF(n, n1 = n1,  n2 = n2, n3 = n3, n4 = n4,
                    mu_mat = mu_mat_1, kappa_vec = kappa_vec1)
mu_mat_2 <- cbind(drop(DirStats::to_cir(0)),
                drop(DirStats::to_cir(-pi)),
                drop(DirStats::to_cir(-pi / 2)),
                drop(DirStats::to_cir(3 * pi / 2.5)))
kappa_vec2 <- c(40, 85, 35, 20)
set.seed(4)
sph_2 <- mix_3_mvMF(n, n1 = n1,  n2 = n2, n3 = n3, n4 = n4,
                    mu_mat = mu_mat_2, kappa_vec = kappa_vec2)
mu_mat_3 <- cbind(drop(DirStats::to_cir(0)),
                drop(DirStats::to_cir(-pi / 4)),
                drop(DirStats::to_cir(pi / 4)),
                drop(DirStats::to_cir(3 * pi / 2.5)))
kappa_vec3 <- c(65, 105, 28, 25)
set.seed(5)
sph_3 <- mix_3_mvMF(n, n1 = n1,  n2 = n2, n3 = n3, n4 = n4,
                    mu_mat = mu_mat_3, kappa_vec = kappa_vec3)
data <- abind(sph_1, sph_2, sph_3, along = 3)
set.seed(42)
indexes <- sample(1:n)
mix_3_mvMF_data_7 <- data[indexes, , , drop = FALSE]
mix_3_mvMF_cols_7 <- c(rep(1, n1), rep(2, n2), rep(3, n3), rep(4, n4))
mix_3_mvMF_cols_7 <- mix_3_mvMF_cols_7[indexes]

## ----fig.asp = 1--------------------------------------------------------------
mix_3_mvMF_3dtorus <- sdetorus::toPiInt(cbind(
  DirStats::to_rad(mix_3_mvMF_data_7[, , 1]),
  DirStats::to_rad(mix_3_mvMF_data_7[, , 2]),
  DirStats::to_rad(mix_3_mvMF_data_7[, , 3])))

pairs(mix_3_mvMF_3dtorus, col = mix_3_mvMF_cols_7)

## ----fig.asp = 1--------------------------------------------------------------
scatterplot3d::scatterplot3d(mix_3_mvMF_3dtorus, xlim = c(-pi, pi),
                             ylim = c(-pi, pi), zlim = c(-pi, pi),
                             color = mix_3_mvMF_cols_7, xlab = "", ylab = "",
                             zlab = "", tick.marks = FALSE)

## ----cache = TRUE-------------------------------------------------------------
rho_30_7 <- rho_optim_bst(x = mix_3_mvMF_data_7, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----size = "scriptsize", fig.align = 'center'--------------------------------
res_pscsne_71 <- psc_sne(X = mix_3_mvMF_data_7, d = 1, rho_psc_list = rho_30_7,
                         eta = 50, colors = mix_3_mvMF_cols_7, show_prog = TRUE,
                         parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_71$best_Y
plot(Y[, 1], Y[, 2], col = mix_3_mvMF_cols_7, xlim = c(-1, 1), ylim = c(-1, 1),
     axes = FALSE, xlab = "", ylab = "")
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
set.seed(42)
psc_sne_res_72 <- psc_sne(mix_3_mvMF_data_7, d = 2, rho_psc_list = rho_30_7,
                          eta = 10, colors = mix_3_mvMF_cols_7,
                          show_prog = TRUE, parallel_cores = num_cores_param,
                          maxit = 2000)

## ----fig.asp = 1, fig.align='center'------------------------------------------
Y <- psc_sne_res_72$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = mix_3_mvMF_cols_7, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(Y[, 2]) == 1, 1, 2)], grid = FALSE,
  mar = c(0, 0, 0, 0)
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----fig.asp = 1--------------------------------------------------------------
n <- 500
mean1 <- c(1, 1, 2)
mean2 <- c(-2, -1, -1)
mean3 <- c(3, 3, 3.25)
sigma1 <- matrix(c(0.22, -0.10, -0.15, -0.10, 0.28, 0.10, -0.15, 0.10, 0.30),
                 byrow = TRUE, nrow = 3)
sigma2 <- toeplitz(c(0.32, 0.1, 0.25))
sigma3 <- toeplitz(c(0.25, -0.15, 0.21))
n1 <- ceiling(n / 3)
n2 <- ceiling(n / 3)
n3 <- floor(n / 3)
set.seed(42)
x_8_bvnorm <- mix_3_mvnorm(n = n, n1 = n1, n2 = n2, n3 = n3,
                           mean_mat = cbind(mean1, mean2, mean3),
                           sigma1 = sigma1, sigma2 = sigma2, sigma3 = sigma3)
x_8_torus <- sdetorus::toPiInt(x_8_bvnorm)

samp_8 <- abind(DirStats::to_cir(x_8_torus[, 1]),
                DirStats::to_cir(x_8_torus[, 2]),
                DirStats::to_cir(x_8_torus[, 3]), along = 3)
samp_8_cols <- c(rep(1, n1), rep(2, n2), rep(3, n3))


## ----fig.asp = 1, fig.align = 'center'----------------------------------------
pairs(x_8_torus, col = samp_8_cols)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
scatterplot3d::scatterplot3d(x_8_torus, xlim = c(-pi, pi), ylim = c(-pi, pi),
                             zlim = c(-pi, pi), xlab = "", ylab = "", zlab = "",
                             color = samp_8_cols, pch = 16, tick.marks = FALSE)

## -----------------------------------------------------------------------------
set.seed(42)
indexes <- sample(1:n)
samp_8 <- samp_8[indexes, , , drop = FALSE]
samp_8_cols <- samp_8_cols[indexes]

## -----------------------------------------------------------------------------
rho_30_8 <- rho_optim_bst(x = samp_8, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----size = "scriptsize", fig.align = 'center'--------------------------------
set.seed(42)
res_pscsne_81 <- psc_sne(X = samp_8, d = 1, rho_psc_list = rho_30_8,
                         eta = 50, colors = samp_8_cols, show_prog = TRUE,
                         parallel_cores = num_cores_param, init = "random")

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_81$best_Y
plot(Y[, 1], Y[, 2], col = samp_8_cols, xlim = c(-1, 1), ylim = c(-1, 1),
     axes = FALSE, xlab = "", ylab = "")
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_82 <- psc_sne(samp_8, d = 2, rho_psc_list = rho_30_8,
                          eta = 50, colors = samp_8_cols,
                          show_prog = TRUE, parallel_cores = num_cores_param)

## ----fig.asp = 1, fig.align='center'------------------------------------------
Y <- psc_sne_res_82$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = samp_8_cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(Y[, 2]) == 1, 1, 2)], grid = FALSE,
  mar = c(0, 0, 0, 0)
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
n <- 500
set.seed(4)
samp_9 <- r_path_s1r(n = n, r = 3, angles = TRUE, k = c(1, 2, 3))
samp_9_cols <- rainbow(n, alpha = 1)

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
pairs(samp_9, xlim = c(-pi, pi), ylim = c(-pi, pi), col = samp_9_cols,
      pch = 16, labels = c("Sphere 1", "Sphere 2", "Sphere 3"))

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
scatterplot3d::scatterplot3d(samp_9, xlim = c(-pi, pi), ylim = c(-pi, pi),
                             zlim = c(-pi, pi), xlab = "", ylab = "", zlab = "",
                             color = samp_9_cols, pch = 16, tick.marks = FALSE)

## -----------------------------------------------------------------------------
x_1 <- DirStats::to_cir(samp_9[, 1])
x_2 <- DirStats::to_cir(samp_9[, 2])
x_3 <- DirStats::to_cir(samp_9[, 3])
x_path_2 <- abind(x_1, x_2, x_3, along = 3)
set.seed(42)
indexes <- sample(1:n)
x_path_2 <- x_path_2[indexes, , , drop = FALSE]
samp_9_cols <- samp_9_cols[indexes]

## ----cache = TRUE-------------------------------------------------------------
rho_30_9 <- rho_optim_bst(x = x_path_2, perp_fixed = 30,
                          num_cores = num_cores_param)

## ----size = "scriptsize", fig.align = 'center'--------------------------------
set.seed(42)
res_pscsne_91 <- psc_sne(X = x_path_2, d = 1, rho_psc_list = rho_30_9,
                         eta = 25, colors = samp_9_cols, show_prog = TRUE,
                         parallel_cores = num_cores_param, init = "random")

## ----fig.asp = 1, fig.align = 'center'----------------------------------------
Y <- res_pscsne_91$best_Y
plot(Y[, 1], Y[, 2], col = samp_9_cols, xlim = c(-1, 1), ylim = c(-1, 1),
     axes = FALSE, xlab = "", ylab = "")
th <- seq(0, 2 * pi, length.out = 100)
polygon(x = cos(th), y = sin(th))

## ----cache = TRUE, fig.height = 8, fig.width = 8, size = "scriptsize", fig.align = 'center'----
psc_sne_res_92 <- psc_sne(x_path_2, d = 2, rho_psc_list = rho_30_9,
                          eta = 10, colors = samp_9_cols,
                          show_prog = TRUE, parallel_cores = num_cores_param, 
                          maxit = 3000)

## ----fig.asp = 1, fig.align='center'------------------------------------------
Y <- psc_sne_res_92$best_Y
sd3 <- scatterplot3d::scatterplot3d(
  Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = samp_9_cols, xlab = "", ylab = "", zlab = "", axis = FALSE,
  pch = c("+", "-")[ifelse(sign(Y[, 2]) == 1, 1, 2)], grid = FALSE,
  mar = c(0, 0, 0, 0)
)
sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
             type = "l", lty = 3)
sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
             type = "l", lty = 3)

