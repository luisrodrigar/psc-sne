library(sdetorus)
n <- 10
x1 <- sdetorus::toPiInt(cbind(
  DirStats::to_rad(rotasym::r_vMF(n = n, mu = c(0, 1), kappa = 10)),
  DirStats::to_rad(rotasym::r_vMF(n = n, mu = c(1, 0), kappa = 10))
  ))
x2 <- sdetorus::toPiInt(cbind(
  DirStats::to_rad(rotasym::r_vMF(n = n, mu = c(-1, 0), kappa = 3)),
  DirStats::to_rad(rotasym::r_vMF(n = n, mu = c(-1, 0), kappa = 3))
))
x <- rbind(x1, x2)
plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
     col = rep(c(1, 2), each = n), pch = 16)
sdetorus::torusAxis()

# Sample spherical Cauchy with mean e_p = (0, 0, …, 1) (well, one that is
# very similar) on S^d
r_sc <- function(n, d = 1, rho = 0.5) {
  r_alt(n=1e3,  p=3, alt='MvMF', kappa=50)
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
