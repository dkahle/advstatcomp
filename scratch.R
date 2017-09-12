f <- deriv(~ x^2 + y^2 + a * x * y, c("x", "y"), function.arg = TRUE)
a <- 1

n <- 40
xpts <- seq(-3, 2, len = n)
ypts <- seq(-2, 3, len = n)
gr <- expand.grid(x = xpts, y = ypts)
feval <- with(gr, f(x, y))
z <- matrix(feval, nrow = n, ncol = n)

par(mar = c(5, 4, 1, 1))
contour(xpts, ypts, z, nlevels = 20)
x0 <- -2.5
y0 <- 1.2
points(x0, y0, pch = 19, cex = 2)

f0 <- f(x0, y0)
p0 <- drop(-attr(f0, "gradient"))
f.sub <- function(alpha) {
        ff <- f(x0 + alpha * p0[1], y0 + alpha * p0[2])
        as.numeric(ff)
}
op <- optimize(f.sub, c(0, 4))
alpha <- op$minimum

arrows(x0, y0, 0, 0, lwd = 3, col = "grey")

x1 <- x0 + alpha * p0[1]
y1 <- y0 + alpha * p0[2]
arrows(x0, y0, x1, y1, lwd = 2)


f1 <- f(x1, y1)
## Fletcher-Reeves
f1g <- drop(attr(f1, "gradient"))
beta <- drop(crossprod(f1g) / crossprod(p0))
p1 <- -f1g + beta * p0
f.sub <- function(alpha) {
        ff <- f(x1 + alpha * p1[1], y1 + alpha * p1[2])
        as.numeric(ff)
}
op <- optimize(f.sub, c(0, 4))
alpha <- op$minimum
x2 <- x1 + alpha * p1[1]
y2 <- y1 + alpha * p1[2]
arrows(x1, y1, x2, y2, lwd = 2)

## Steepest descent only
p1 <- -f1g
op <- optimize(f.sub, c(0, 4))
alpha <- op$minimum

x2 <- x1 + alpha * p1[1]
y2 <- y1 + alpha * p1[2]
arrows(x1, y1, x2, y2, col = "red", lwd = 2, lty = 2)



#########################################################################
## EM Algorithm

mu1 <- 1
s1 <- 2
mu2 <- 4
s2 <- 1

lambda0 <- 0.4
n <- 100
set.seed(2017-09-12)
z <- rbinom(n, 1, lambda0)
x <- rnorm(n, mu1 * z + mu2 * (1-z), s1 * z + (1-z) * s2)
hist(x)
rug(x)

f <- function(x, lambda) {
        lambda * dnorm(x, mu1, s1) + (1-lambda) * dnorm(x, mu2, s2)
}
nll <- function(lambda) {
        sum(log(f(x, lambda)))
}
nll <- Vectorize(nll, "lambda")

curve(nll, 0.01, 0.95, n = 200)

lam0 <- 0.2
minor <- function(lambda) {
        p1 <- sum(log(f(x, lam0)))
        pi <- lam0 * dnorm(x, mu1, s1) / (lam0 * dnorm(x, mu1, s1) 
                                          + (1 - lam0) * dnorm(x, mu2, s2))
        p2 <- sum(pi * dnorm(x, mu1, s1, log = TRUE) 
                  + (1-pi) * dnorm(x, mu2, s2, log = TRUE)
                  + pi * log(lambda)
                  + (1-pi) * log(1-lambda))
        p3 <- sum(pi * dnorm(x, mu1, s1, log = TRUE) 
                  + (1-pi) * dnorm(x, mu2, s2, log = TRUE)
                  + pi * log(lam0)
                  + (1-pi) * log(1-lam0))
        p1 + p2 - p3
}
minor <- Vectorize(minor, "lambda")

curve(minor, 0.01, 0.99, add = TRUE, col = "red")
lam0 <- 0.1
curve(minor, 0.01, 0.99, add = TRUE, col = "red")
lam0 <- 0.9
curve(minor, 0.01, 0.99, add = TRUE, col = "blue")
lam0 <- 0.4
curve(minor, 0.01, 0.99, add = TRUE, col = "blue")
abline(v = op$minimum)








