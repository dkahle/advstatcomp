f <- function(x, y, a = 1) {
        x^2 + y^2 + a * x * y
}

f <- deriv(~ x^2 + y^2 + a * x * y, c("x", "y"), function.arg = TRUE)
a <- 1

n <- 40
xpts <- seq(-3, 3, len = n)
ypts <- seq(-3, 3, len = n)
gr <- expand.grid(x = xpts, y = ypts)
feval <- with(gr, vf(x, y, 1))
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

arrows(x0, y0, 0, 0, lwd = 3)

x1 <- x0 + alpha * p0[1]
y1 <- y0 + alpha * p0[2]
arrows(x0, y0, x1, y1, col = "grey", lwd = 2)


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
arrows(x1, y1, x2, y2, col = "gray", lwd = 2)

## Steepest descent only
p1 <- -f1g
op <- optimize(f.sub, c(0, 4))
alpha <- op$minimum

x2 <- x1 + alpha * p1[1]
y2 <- y1 + alpha * p1[2]
arrows(x1, y1, x2, y2, col = "red", lwd = 2, lty = 2)






