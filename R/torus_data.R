library(sphunif)
library(DirStats)
library(sdetorus)
library(rotasym)
library(scatterplot3d)

# Sample on the (S^1)^2
n <- 200
x1 <- sdetorus::toPiInt(cbind(
  DirStats::to_rad(rotasym::r_vMF(n = n, mu = drop(DirStats::to_cir(th = 0)), kappa = 10)),
  DirStats::to_rad(rotasym::r_vMF(n = n, mu = drop(DirStats::to_cir(th = pi)), kappa = 10))
))
x2 <- sdetorus::toPiInt(cbind(
  DirStats::to_rad(rotasym::r_vMF(n = n, mu = drop(DirStats::to_cir(th = pi / 2)), kappa = 5)),
  DirStats::to_rad(rotasym::r_vMF(n = n, mu = drop(DirStats::to_cir(th = 0)), kappa = 5))
))
x <- rbind(x1, x2)
plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
     col = rep(c(1, 2), each = n), pch = 16)
sdetorus::torusAxis()

# Cartesian coordinates
x <- DirStats::to_cir(x)
x_array <- array(dim = c(2 * n, 2, 2))
x_array[, , 1] <- x[, 1:2]
x_array[, , 2] <- x[, 3:4]

# Sample spherical Cauchy with mean e_p = (0, 0, …, 1) (well, one that is
# very similar) on S^d
r_sc <- function(n, d = 1, rho = 0.5) {
  sphunif::r_alt(n = n, p = d + 1, M = 1, alt = "C",
                 kappa = rho / (1 - rho^2))[, , 1]
}

# Sample polyspherical Cauchy on S^{d,r}
r_psc <- function(n, d = 1, r = 2, rho = rep(0.5, r)) {
  samp <- matrix(nrow = n, ncol = r * (d + 1))
  for (k in 1:r) {
    ind_k <- (d + 1) * (k - 1) + 1:(d + 1)
    samp[, ind_k] <- r_sc(n = n, d = d, rho = rho[k])
  }
  return(samp)
}

# Sample data on S^1 x S^1
x <- r_psc(n = 1000, d = 1, r = 2, rho = rep(0.8, 2))
x <- sdetorus::toPiInt(cbind(DirStats::to_rad(x[, 1:2]),
                             DirStats::to_rad(x[, 3:4])))
plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE, pch = 16)
sdetorus::torusAxis()

# Von Mises with any mean
n <- 1e3
samp1 <- rotasym::r_vMF(n = n, mu = drop(DirStats::to_sph(th = 0, ph = 0.5)),
                        kappa = 50)
samp2 <- rotasym::r_vMF(n = n, mu = drop(DirStats::to_sph(th = 2, ph = -1.5)),
                        kappa = 50)
samp3 <- rotasym::r_vMF(n = n, mu = drop(DirStats::to_sph(th = -1, ph = 0)),
                        kappa = 50)
scatterplot3d::scatterplot3d(rbind(samp1, samp2, samp3),
                             xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
                             color = rep(1:3, each = n))

# Mixture of von Mises
samp <- sphunif::r_alt(n = n, p = 3, alt = "MvMF", kappa = 50)[, , 1]
scatterplot3d::scatterplot3d(samp, xlim = c(-1, 1), ylim = c(-1, 1),
                             zlim = c(-1, 1))

# Mixing small circles
samp <- sphunif::r_alt(n = n, p = 3, alt = "SC", kappa = 50, nu = 0)[, , 1]
scatterplot3d::scatterplot3d(samp, xlim = c(-1, 1), ylim = c(-1, 1),
                             zlim = c(-1, 1))
scatterplot3d::scatterplot3d(rbind(samp, samp[, c(3, 1, 2)]), xlim = c(-1, 1),
                             ylim = c(-1, 1), zlim = c(-1, 1),
                             color = rep(1:2, each = n))

# Other generators:
# ?rotasym::r_ACG(), ?sphunif::r_alt()