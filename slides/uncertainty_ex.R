


my = mean(d$y)

nx = c(-2, -1,1)
lm_res = lm(z ~ x, data = d)

set.seed(12)

f = rstanarm::stan_glm(z ~ x, data = d) 

ddf = f |> rstanarm::as_draws()

xv = seq(min(d$x), max(d$x), length.out = 100)

m = predict(f, xv)

bds = mapply(\(i,b,s) i + b*xv, ddf$`(Intercept)`, ddf$x, ddf$sigma) |> 
    dapply(fquantile, probs = c(.025, .975),
           MARGIN = 1)

pts_df = ddf |> head() |> mtt(pts = mapply(\(i,b,s) i + b*c(-2,-1,1) + rnorm(3)*s, `(Intercept)`, x, s, SIMPLIFY = FALSE)) |> 
    get_elem("pts") |> 
    unlist2d() |> 
    pivot(".id") |> 
    mtt(x = rep(c(-2,-1,1), each = 6))

png("prop_ex.png",
    width=960*2,height=720*2, res = 300)

  plt(d$x, d$z, 
    xlab = "", 
    ylab = "", 
    col = mblue,
    cex = 1.8,
    ylim = yr)
abline(coef(f)[1], coef(f)[2],
       col = mred,
       lwd=1.3)

points(pts_df$x, pts_df$value,
       col = mred,
       cex=2)
polygon(c(xv,rev(xv)), c(bds[,1], rev(bds[,2])),
        border = NA, col = mred |> scales::alpha(.2))

dev.off()
