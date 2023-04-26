
library(devtools)
library(covr)
library(lintr)
library(rcmdcheck)
library(here)

# Package name
pkg <- "pscsne"

# # Install from GitHub
# remove.packages(pkg)
# devtools::install_github("https://github.com/luisrodrigar/psc-sne/",
#                           auth_token = paste0("ghp_GmbcUYxjKPJaDez5e",
#                                               "5hTf8TC0BdTzO0qZJJX"),
#                          ref = "create-package")
# packageVersion(pkg)

# Set wd
setwd(here())

# Load package
devtools::load_all(path = ".")

# Documentation
document(pkg = ".")

# Tests
test(pkg = ".")

# Run examples
run_examples(".", run_dontrun = TRUE, run_donttest = TRUE, fresh = TRUE)

# Check
# check(pkg = ".", cran = TRUE)
# check(pkg = ".", cran = TRUE, args = c("--no-tests")) # No tests
check(pkg = ".", args = "--no-examples", vignettes = FALSE) # No examples nor vignettes
# check(pkg = pkg, document = FALSE) # Examples and vignettes
# check(pkg = pkg, document = FALSE, args = "--run-donttest") # Extended examples

# R CMD check
chk <- rcmdcheck(".")

# lintr for flagging issues
lint_package(path = ".",
             exclusions = list(paste0(pkg, "/R/RcppExports.R")))

# Tests coverage
(covr_pkg <- package_coverage(path = "."))
report(covr_pkg)

# Check spelling
spell_check(pkg = ".")

# Build pdf
path <- find.package(pkg)
file.remove(paste0(pkg, ".pdf"))
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD", "Rd2pdf", shQuote(path)))
# build_manual(path = ".")

# Build
build(".", vignettes = TRUE)

# Install locally
install(".", args = c("--no-multiarch", "--no-test-load"),
        build_vignettes = FALSE)
packageVersion(pkg)

# Load package
library(devtools)
setwd(paste("/Users/Eduardo/GitHub/", pkg, "/", sep = ""))
do.call(library, list(pkg))

# Checks before release()
# as.data.frame(rhub::platforms())[, 1]
plts <- c("windows-x86_64-release", # Windows Server 2008 R2 SP1, R-release, 32/64 bit
          "windows-x86_64-devel", # Windows Server 2008 R2 SP1, R-devel, 32/64 bit
          "windows-x86_64-oldrel", # Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
          "windows-x86_64-patched", # Windows Server 2008 R2 SP1, R-patched, 32/64 bit
          "macos-highsierra-release", # macOS 10.13.6 High Sierra, R-release, brew
          "macos-highsierra-release-cran", # macOS 10.13.6 High Sierra, R-release, CRAN's setup
          "ubuntu-gcc-release", # Ubuntu Linux 20.04.1 LTS, R-release, GCC
          "ubuntu-gcc-devel", # Ubuntu Linux 20.04.1 LTS, R-devel, GCC
          "debian-gcc-release", # Debian Linux, R-release, GCC
          "debian-gcc-devel" # Debian Linux, R-devel, GCC
)
check_rhub(pkg = pkg, interactive = FALSE, platforms = plts,
           check_args = "--as-cran",
           env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false",
                        `_R_CHECK_DONTTEST_EXAMPLES_` = "false"))
check_win_release(pkg = pkg)
check_win_devel(pkg = pkg)
compare_to_cran(chk)

# Fingers crossed
release(pkg)



png(filename = "file.png", width = 7, height = 7, units = "in", res = 300)
plot(1:5)
dev.off()

