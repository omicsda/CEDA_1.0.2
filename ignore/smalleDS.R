load("./mda231.RData")
class(mda231)
length(mda231)
names(mda231)
lapply(mda231, dim)

set.seed(39755)
N <- nrow(mda231$sgRNA)
target <- round(N/3)
picker <- sample(N, target)
packer <- sample(350, 200)
smallmda231 <- list(sgRNA = mda231$sgRNA[picker,],
                    neGene = mda231$neGene[packer, , drop = FALSE])
lapply(smallmda231, dim)
mda231 <- smallmda231
rm(smallmda231)
save(mda231, file = "../data/mda231.RData", compress = "xz")
