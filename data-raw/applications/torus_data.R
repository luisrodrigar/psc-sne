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
plot(x,
  xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
  col = rep(c(1, 2), each = n), pch = 16
)
sdetorus::torusAxis()

# Cartesian coordinates
# x_array <- array(dim = c(2 * n, 2, 2))
# x_array[, , 1] <- DirStats::to_cir(x[, 1])
# x_array[, , 2] <- DirStats::to_cir(x[, 2])

# Sample spherical Cauchy with mean e_p = (0, 0, …, 1) (well, one that is
# very similar) on S^d
r_sc <- function(n, d = 1, rho = 0.5) {
  sphunif::r_alt(
    n = n, p = d + 1, M = 1, alt = "C",
    kappa = rho / (1 - rho^2)
  )[, , 1]
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
x <- sdetorus::toPiInt(cbind(
  DirStats::to_rad(x[, 1:2]),
  DirStats::to_rad(x[, 3:4])
))
plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE, pch = 16)
sdetorus::torusAxis()

# Von Mises with any mean
n <- 1e3
samp1 <- rotasym::r_vMF(
  n = n, mu = drop(DirStats::to_sph(th = 0, ph = 0.5)),
  kappa = 50
)
samp2 <- rotasym::r_vMF(
  n = n, mu = drop(DirStats::to_sph(th = 2, ph = -1.5)),
  kappa = 50
)
samp3 <- rotasym::r_vMF(
  n = n, mu = drop(DirStats::to_sph(th = -1, ph = 0)),
  kappa = 50
)
scatterplot3d::scatterplot3d(rbind(samp1, samp2, samp3),
  xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  color = rep(1:3, each = n)
)

# Mixture of von Mises
samp <- sphunif::r_alt(n = n, p = 3, alt = "MvMF", kappa = 50)[, , 1]
scatterplot3d::scatterplot3d(samp,
  xlim = c(-1, 1), ylim = c(-1, 1),
  zlim = c(-1, 1)
)

# Mixing small circles
samp <- sphunif::r_alt(n = n, p = 3, alt = "SC", kappa = 50, nu = 0)[, , 1]
scatterplot3d::scatterplot3d(samp,
  xlim = c(-1, 1), ylim = c(-1, 1),
  zlim = c(-1, 1)
)
scatterplot3d::scatterplot3d(rbind(samp, samp[, c(3, 1, 2)]),
  xlim = c(-1, 1),
  ylim = c(-1, 1), zlim = c(-1, 1),
  color = rep(1:2, each = n)
)

# Other generators:
# ?rotasym::r_ACG(), ?sphunif::r_alt()

rotate_matrix_z_axis <- function(alpha) {
  deg2rad <- function(deg) {(deg * pi) / (180)}
  rads <- deg2rad(alpha)
  matrix(c(cos(rads), -sin(rads), 0,
           sin(rads),  cos(rads), 0,
           0,                  0, 1),
         byrow = T, nrow = 3)
}

p = 2
A <- matrix(runif((p + 1)^2)*2-0.5, ncol = p + 1)
sigma <- t(A) %*% A
n <- 200
x_1_belt1 <- rotasym::r_ACG(n, sigma)
rotate_mat_90 <- rotate_matrix(90)
x_2_belt1 <- x_1_belt1 %*% rotate_mat_90
rotate_mat_45 <- rotate_matrix(45)
x_3_belt1 <- x_1_belt1 %*% rotate_mat_45
x_belt1 <- rbind(x_1_belt1, x_2_belt1, x_3_belt1)
scatterplot3d::scatterplot3d(x_belt1, xlim = c(-1, 1),
                             ylim = c(-1, 1), zlim = c(-1, 1),
                             color = rep(c(1,2,3), each=n))

rgl::plot3d(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
            radius = 1, type = "s", col = "lightblue", alpha = 0.25,
            lit = FALSE)
rgl::points3d(x_belt1, col = rep(c(1,2,3), each=n))

x_1_belt2 <- rotasym::r_ACG(n, sigma)
x_2_belt2 <- x_1_belt2 %*% rotate_mat_90
x_3_belt2 <- x_1_belt1 %*% rotate_mat_45
x_belt2 <- rbind(x_1_belt2, x_2_belt2, x_3_belt2)
scatterplot3d::scatterplot3d(x_belt2, xlim = c(-1, 1),
                             ylim = c(-1, 1), zlim = c(-1, 1),
                             color = rep(c(1,2,3), each=n))

rgl::plot3d(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
            radius = 1, type = "s", col = "lightblue", alpha = 0.25,
            lit = FALSE)
rgl::points3d(x_belt2, col = rep(c(1,2,3), each=n))

x_1_belt3 <- rotasym::r_ACG(n, sigma)
x_2_belt3 <- x_1_belt3 %*% rotate_mat_90
x_3_belt3 <- x_1_belt1 %*% rotate_mat_45
x_belt3 <- rbind(x_1_belt3, x_2_belt3, x_3_belt3)
scatterplot3d::scatterplot3d(x_belt3, xlim = c(-1, 1),
                             ylim = c(-1, 1), zlim = c(-1, 1),
                             color = rep(c(1,2,3), each=n))

rgl::plot3d(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
            radius = 1, type = "s", col = "lightblue", alpha = 0.25,
            lit = FALSE)
rgl::points3d(x_belt3, col = rep(c(1,2,3), each=n))

