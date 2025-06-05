## peak_simu.R. -- function for doing coale scenario and extensions
## this one has 2 parameters k1 and k2

## to do:
## (1) add check box for NRR slope in Births

leslie1 <- function(Lx, Fx)
{
    nLx <- Lx
    nFx <- Fx
    ## subdiagonal
    nLxpn <- c(nLx[-1], 0)              # nL(x+n), next age-group
    subdi <- (nLxpn/nLx)[-length(nLx)]  # survivorship sub-diagonal of
                                        # Leslie matrix
    ## first row
    nFxpn <- c(nFx[-1], 0)              # nF(x+n), next age group of fert.
    nL0 <- nLx[1]
    firstrow <- (nL0/2)*(nFx + nFxpn * nLxpn/nLx)*0.4886
    ## put subdiagonal and first row into matrix format
    A <- rbind(firstrow,
               cbind(diag(subdi), 0))
    dimnames(A) <- NULL                 # remove messy labels
    return(A)
}

project.leslie1 <- function(nsteps,
                            nrr.vec,
                            e0.vec,
                            Lx.zero, Fx.zero)
{
    ## make K.mat (we do this with a stable pop from starting A)
    new.Lx <- ezero2Lx(e0 = e0.vec[1],
                       Lx.zero)
    new.Fx <- nrr2Fx(nrr = nrr.vec[1],
                     Fx.zero,
                     Lx = new.Lx)

    ## make A, (leslie1)
    new.A <- leslie1(Lx = new.Lx, Fx = new.Fx)
    Kx.zero <- eigen(new.A)$vec[,1]
    Kx.zero <- Kx.zero/sum(Kx.zero)
    K.mat <- matrix(0,
                    nrow = length(Fx.zero),
                    ncol = length(nrr.vec))
    K.mat[,1] <- Re(Kx.zero)



    ## do time varying projection
    for (i in 1:(nsteps-1))
    {
        ## get Lx and Fx to match e0 and nrr
        new.Lx <- ezero2Lx(e0 = e0.vec[i],
                           Lx.zero)
        new.Fx <- nrr2Fx(nrr = nrr.vec[i],
                         Fx.zero,
                         Lx = new.Lx)

        ## make A, (leslie1)
        new.A <- leslie1(Lx = new.Lx, Fx = new.Fx)

        ## project
        K.mat[,i+1] <- new.A %*% K.mat[,i]
    }
    return(K.mat)
}

nrr2Fx <- function(nrr, Fx.zero, Lx)
{
    nrr.zero <- sum(Fx.zero * Lx * .4886)
    rat <- nrr/nrr.zero
    Fx <- Fx.zero * rat
    Fx
}

if(0) {
library(data.table)
dt.mort <- fread("~/Documents/hmd/hmd_statistics/lt_both/bltper_1x1/SWE.bltper_1x1.txt")
dt.fert <- fread("~/Documents/hfd/zip_w/asfrRR.txt")
## need to hardcode!!
Lx <- dt.mort[Year == 2014, Lx]/10^5
x <- 0:110
Fx.orig <- dt.fert[Code == "SWE" &  Year == 1900, ASFR]
Fx.ages <- dt.fert[Code == "SWE" &  Year == 1900, Age]
## put in 0:110 ages
Fx <- rep(0, length(x))
Fx[x %in% 12:55] <- Fx.orig
dump(c("Fx","Lx"), file = "data.R")
}

source("data.R")



ezero2Lx <- function(e0, Lx.zero)
{
    f <- function(theta)
    {
        (sum(Lx.zero^theta) - e0)
    }
    theta.hat <- uniroot(f, lower = .00001, upper = 100)$root
    ## do only mort change over age 50, so neutral
    Lx <- Lx.zero^c(rep(0, 50), rep(theta.hat, length(Lx.zero)-50))
    return(Lx)
}



center.diff <- function(x, end.fill = F)
{
    ## approximate derivatives with discrete data by taking central
    ## differences d(x) = ([f(x+1) - f(x)] + [f(x) - f(x-1)])/2 =
    ## [f(x+1) - f(x-1)]/2 if end.fill = T, then first and last
    ## differences are not centered.

    ## useful for Bongaarts-Feeney, Gompertz fertility, and other
    ## fitting of models that are developed in continuous time and
    ## involve derivatives

    forward.diff = c(diff(x), NA)
    backward.diff = c(NA, diff(x))
    out = (forward.diff + backward.diff)/2
    if(end.fill)
    {
        out[1] = diff(x)[1]
        out[length(out)] = diff(x)[length(diff(x))]
    }
    ## preserve names if exist
    if(!is.null(names(x)))
        names(out) = names(x)
    return(out)
}

get_approx_peak <- function(t, ft)
{
    dft <- center.diff(ft)
    peak <- approx(x = dft, y = t, xout = 0)$y
    return(peak)
}

## debug(coale_simu_movie_frame_4)
## coale_simu_movie_frame_4(k1 = -.02, k2 = 0,
##                          rho = .002,
##                          Lx = Lx,
##                          Fx = Fx,
##                          t = -100:100,
##                          myt = 50,
##                          nrr_slope = FALSE,
##                          show_population = TRUE,
##                          show_prediction = TRUE,
##                          curvature = FALSE)

coale_simu_movie_frame_4 <-
    function(k1, k2, rho = 0,
             Lx, Fx, t, myt,
             nrr_slope, show_population,
             show_prediction,
             curvature)
{
    k1 <- as.numeric(k1)
    k2 <- as.numeric(k2)
    radix = 1000
    ## do projection
    nrr.t <- exp(k1 * t  + k2 * t^2)

    ## rho
    rho <- as.numeric(rho)
    e0.vec <- exp(rho*t) * 80
    ## get prediction
    if (show_prediction) {
    x <- seq(Fx) - 1
    mu0 = sum((x + .5)*Lx*Fx)/sum(Lx*Fx)
    A0 = sum((x + .5)*Lx)/sum(Lx)
    A2 = sum((x + .5)^2*Lx)/sum(Lx)
    S2 <- A2 - A0^2
    print(mu0)
    print(S2)
    print(k1)
    print(k2)
    t0 = 0
    tB.hat = t0 - mu0/2
    tN.hat = tB.hat + A0 + (-rho/k1) * mu0 -
        (2*k2/k1)*S2/2
    print(tN.hat)
    }
    K.mat.out <- project.leslie1(nsteps = length(t),
                                 nrr.vec = nrr.t,
                                 e0.vec = e0.vec,
                                 Lx.zero = Lx,
                                 Fx.zero = Fx)
    ## normalize population size to 1.0 at time t = 0
    Nt.unnorm <- colSums(K.mat.out)
    K.mat.out <- K.mat.out/   Nt.unnorm[t == 0]
    K.mat.out <- radix*K.mat.out
    ## Calculate pop size Nt, births Bt
    Nt <- colSums(K.mat.out)
    Bt <- K.mat.out[1,]
    ## Get dates of max pop and max births
    ## first decline
    tm <- t[min(which(center.diff(Nt) < 0))]
##    print(tm)
    tN.exact <- get_approx_peak(t=t[t < tm], ft=Nt[t < tm])
    tB.exact <- get_approx_peak(t=t[t < 0], ft=Bt[t < 0])
    tN.exact <- round(tN.exact, 2)
    tB.exact <- round(tB.exact, 2)

    ##############
    ## plotting ##
    ##############

    my.xlim = c(-50, 50)
    s <- t <= myt
    par(mfrow = c(2,2),
        cex = 1.3,
        las = 1,
        mar = c(5.1-1, 4.1, 4.1+1, 2.1))

    #########
    ## NRR ##
    #########
    plot(t, nrr.t, type = 'l',
         lwd = 3, col = "grey",
         xlim = my.xlim,
         ylim = c(0, 2.7), axes = F,
         ylab = "",
         xlab = "time")
    lines(t, exp(k1*t), lty = 2, col = "grey", lwd = 3)
    axis(1, at = seq(-50, 50, 25))
    axis(2, at = c(0, 1, 2, 3))
    lines(t[s], nrr.t[s],
          lwd = 3, col = "turquoise")
    abline(v = 0)
    abline(h = 1)
    grid()
    title("Net reproduction rate")
    if(!curvature)
    title(line = +.8,
          bquote(NRR(t)==1*e^{k*t}),
          cex.main = .8)
    else
    title(line = +.8,
          bquote(NRR(t)==1*e^{k[1]*t~+~k[2]*t^2}),
          cex.main = .8)


    ## title(paste0("Net Reproduction Rate: ",
    ##              "NRR(t) = 1*exp(k*t)"))
    points(x = myt, nrr.t[t == myt])
    text(-30, 1, "replacement", pos = 1, cex = .8)

    ########
    ## e0 ##
    ########
    plot(t, e0.vec, type = 'l', lwd = 3,
         xlim = my.xlim,axes = F,
         ylab = "", xlab = "time",
         ylim = c(60,100))
    axis(1, at = seq(-50, 50, 25))
    axis(2)
##     abline(v = 0)
     grid()
    title("Life expectancy")
    title(line = +.8,
          bquote(e[0](t)==80*e^{rho*t}),
          cex.main = .8)

    ########
    ## Bt ##
    ########
    plot(t, Bt, axes = F, type = 'n',
         lwd = 3, col = "red",
         xlim = my.xlim,
         ylab = "",
         xlab = "time",
         ##          ylim = radix*c(.005, .02))
         ylim = radix*c(.005, .02))
    lines(t[s], Bt[s], lwd = 3, col = "red")
        points(x = myt, Bt[t == myt])
    if (nrr_slope == TRUE) {
    segments(x0 = myt-30, x1 = myt -0,
             y0 = Bt[ t == myt-30],
             y1 = Bt[t == myt-0],
             lty = 2, col = "turquoise", lwd = 3)
    text(x = myt - c(30, 0), y = Bt[t %in% c(myt - c(30,0))],
         c("B(t-30)", "B(t)"), pos = 3, col = "black")

    points(x = myt-30, Bt[ t == myt-30], pch = 19, col = "turquoise")

    }
    axis(1, at = seq(-50, 50, 25))
    axis(2)
    grid()
    abline(v = 0)

    if(myt >= round(tB.exact))
    {
        segments(x0 = tB.exact,
                 y0 = 0, y1 = max(Bt), lty = 2)
        ##         abline(v = tN.exact, lty = 2)
        text(x = tB.exact, y = 7, pos = 2,
             round(tB.exact,1))
        ## abline(v = tB.exact, lty = 2)
        ## axis(3, at = round(tB.exact,1),
        ##      line = -1, tick = FALSE,
        ##      labels = format(round(tB.exact,1), nsmall = 1))
    }
    if(show_prediction & myt >= round(tB.exact))
    {
        segments(x0 = tB.hat,
                 y0 = 0, y1 = max(Bt),
                 lty = 2, col = "orange")
    }
##     title("Births: B(t)")
    ## text(+30, 18.5, expression(
    ##     paste(B(t)==integral(B(t-x)*phi(x,t)*dx))))
    title("Births")
    title(line = +.8,
          bquote(B(t)==integral(B(t-x)*phi(x,t)*dx)),
          cex.main = .8)

    ########
    ## Nt ##
    ########
    if(show_population) {
    plot(t, Nt, axes = F, type = 'n',
         lwd = 3, col = "blue",
         ylab = "",
         ##       ylim = c(0,max(Nt)*1.1))
         xlim = my.xlim,
         xlab = "time",
         ylim = radix*c(.2,1.3))
    lines(t[s], Nt[s],
          lwd = 3, col = "blue")
    points(x = myt, Nt[t == myt])
    axis(1, at = seq(-50, 50, 25))
    axis(2)
    abline(v = 0)
    tN.exact <- get_approx_peak(t=t, ft=Nt)
    ##     print(tN.exact)
    tN <- t[which.max(Nt)]
    ##     print(tN)
    if(myt >= round(tN.exact))
    {
        segments(x0 = tN.exact,
                 y0 = 0, y1 = max(Nt), lty = 2)
        ##         abline(v = tN.exact, lty = 2)
        text(x = tN.exact, y = 400, pos = 4,
             round(tN.exact,1))
        ## axis(3, at = round(tN.exact,1),
        ##      labels = format(round(tN.exact,1), nsmall = 1),
        ##      line = -2, tick = FALSE)
    }
    if(show_prediction & myt >= round(tN.exact))
    {
        segments(x0 = tN.hat,
                 y0 = 0, y1 = max(Nt),
                 lty = 2, col = "orange")

    }

    grid()
##     abline(h = max(Nt), col = "grey")
##     text(-90, max(Nt), round(max(Nt),2),
    ##          pos = 1, col = "grey",)
    ## title("Population: N(t)")
    ## text(-35, 1000, expression(
    ##     paste(N(t)==integral(B(t-x)*l(x,t)*dx))))
##     my_expression <- paste(expression(N(t)==integral(B(t-x)*l(x,t)*dx)))
##     text(-35, 1000, my_expression)
    title("Population")
    title(line = +.8,
          bquote(N(t)==integral(B(t-x)*l(x,t)*dx)),
          cex.main = .8)

    }
    return(NULL)
}


library(shiny)
## runApp() to debug


## Define server logic required to draw plot
library(shiny)
shinyServer(function(input, output) {
    ## Expression that generates a histogram. The expression is
    ## wrapped in a call to renderPlot to indicate that:
    ##
    ##  1) It is "reactive" and therefore should re-execute automatically
    ##    when inputs change
    ##  2) Its output type is a plot
    output$peakPlot <- renderPlot({
        coale_simu_movie_frame_4(k1 = input$k1,
                                 k2 = input$k2,
                                 rho = input$rho,
                               Lx = Lx,
                               Fx = Fx,
                               t = -100:100,
                               myt = input$myt,
                               nrr_slope = input$nrr_slope,
 ##          show_population = input$show_population,
                               show_population = TRUE,
                               show_prediction = input$show_prediction,
                               curvature = input$curvature)

    })
})
