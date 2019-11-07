library(FluMoDL)

cat("Loading temperature data... ")
load("data/mdTemp.RData")
mdTemp$temp <- linterp(mdTemp$temp, 5)
mdTemp$yearweek <- isoweek(mdTemp$date)
cat(sprintf("(temperatures available until week %s)\n", max(mdTemp$yearweek)))


# Maybe trim what is excessive?



cat("Reading in surveillance data... ")
load("../influenza_data/fluData.RData")
cat(sprintf("(ILI rates available until week %s)\n", max(resAll$yearweek)))

load("output/latest_report.RData")
tgtweek <- 201931  # Manually set tgtweek if necessary
rm(isoweek)

totDeaths <- subset(totDeaths, yearweek<=tgtweek)

labConfDeaths <- c(table(subset(datICU, outcome=="death" & season>=2013 & season<tgtyear)$season), tgtyear=nrow(totDeaths))
names(labConfDeaths)[length(labConfDeaths)] <- tgtyear

datLab <- subset(datLabAll, yearweek <= tgtweek)
datLab$pct <- datLab$Positive/datLab$Total
datLab$pctH3 <- datLab[["A(H3N2)"]]/datLab$Total
datLab$pctH1 <- datLab[["A(H1N1)pdm09"]]/datLab$Total
datLab$pctB <- datLab[["B"]]/datLab$Total

datLab <- datLab[,c("yearweek","Positive","pct","pctH1","pctH3","pctB")]



cat("Loading deaths data... ")

load("../MOMO/input/allDeaths.RData")
big <- subset(big, is.na(violentDeathType))

deaths <- as.data.frame(table(big$DoD), responseName="diedAll")
deaths <- merge(deaths, as.data.frame(table(subset(big, age>=75)$DoD), 
    responseName="died75"), by="Var1", all=TRUE)
deaths <- merge(deaths, as.data.frame(table(subset(big, age>=65)$DoD), 
    responseName="died65"), by="Var1", all=TRUE)
deaths <- merge(deaths, as.data.frame(table(subset(big, age>=15 & age<65)$DoD), 
    responseName="died1564"), by="Var1", all=TRUE)
deaths <- merge(deaths, as.data.frame(table(subset(big, age>=5 & age<15)$DoD), 
    responseName="died0514"), by="Var1", all=TRUE)
deaths <- merge(deaths, as.data.frame(table(subset(big, age<5)$DoD), 
    responseName="died04"), by="Var1", all=TRUE)

names(deaths)[1] <- "date"
deaths$date <- as.Date(as.character(deaths$date))
cat(sprintf("(deaths data available until week %s)\n", isoweek(max(deaths$date))))

deaths$yearweek <- isoweek(deaths$date)
weeklyDeaths <- aggregate(deaths[,grep("died", names(deaths))], deaths[,"yearweek", drop=FALSE], sum, na.rm=TRUE)



cat("Merging and filling extra columns...\n")


weekly <- Reduce(function(x,y) merge(x,y,all.x=TRUE), list(resAll[,c("yearweek","gri")], datLab))
weekly$GoldH3 <- weekly$pctH3 * weekly$gri
weekly$GoldH3[is.na(weekly$GoldH3)] <- 0
weekly$GoldH1 <- weekly$pctH1 * weekly$gri
weekly$GoldH1[is.na(weekly$GoldH1)] <- 0
weekly$GoldB <- weekly$pctB * weekly$gri
weekly$GoldB[is.na(weekly$GoldB)] <- 0

weekly$gri[weekly$yearweek<201440] <- weekly$gri[weekly$yearweek<201440]*4.5


daily <- Reduce(function(x,y) merge(x,y,all.x=TRUE), list(mdTemp, deaths))
daily <- daily[order(daily$date),]

weekly <- subset(weekly, yearweek>=min(daily$yearweek))
weekly <- merge(weekly, weeklyDeaths, all.x=TRUE)

daily <- subset(daily, yearweek<=tgtweek)
weekly <- subset(weekly, yearweek<=tgtweek)




cat("\n\nNow analyzing attributable mortality...\n\n")


ageGroups <- c("All","75","65","1564")
popAgeGroups <- c(10768193, 1192225, 2319741, 6893783)

fluMoDLs <- lapply(ageGroups, function(a)
  fitFluMoDL(deaths = daily[[paste0("died",a)]],
             temp = daily$temp, dates = daily$date,
             proxyH1 = weekly$GoldH1, proxyH3 = weekly$GoldH3,
             proxyB = weekly$GoldB, yearweek = weekly$yearweek))
names(fluMoDLs) <- ageGroups

cat("Calculating weekly estimates for the latest season....\n")
inflAttr_weekly <- lapply(ageGroups, function(a) {
  cat(sprintf("Calculating for group %s...\n", a))
  attrMort(fluMoDLs[[a]], sel="week", from=tgtyear*100+40, to=tgtyear*100+120)
})
names(inflAttr_weekly) <- ageGroups


cat("\nCalculating entire-season estimates for the current and past seasons....\n")
inflAttr_seasonal <- lapply(ageGroups, function(a) {
  cat(sprintf("Calculating for group %s...\n", a))
  attrMort(fluMoDLs[[a]], sel="season", to=tgtyear)
})
names(inflAttr_seasonal) <- ageGroups


save(weekly, daily, ageGroups, popAgeGroups, 
    fluMoDLs, inflAttr_weekly, inflAttr_seasonal, 
    totDeaths, labConfDeaths, tgtyear,
    file="output/latest_attrMort.RData")






cat("\nNow plotting output...\n")

library(plotrix)

plotYrFluDeaths <- function(family="Fira Sans") {
  par(family=family)
  ylim <- range(pretty(range(inflAttr_seasonal$All[,c("AllFlu","AllFlu.lo","AllFlu.hi")])))  
  ps <- barplot(rbind(inflAttr_seasonal$All$AllFlu, labConfDeaths),
    ylim=ylim, xlab="Περίοδος επιτήρησης γρίπης",
    border="white", col=c("purple","orange"),
    ylab="Αριθμός θανάτων", beside=TRUE,
    names.arg=rep(NA,length(labConfDeaths)))

  legend("topright", c("Εκτιμώμενοι θάνατοι αποδιδόμενοι στη γρίπη", "Καταγεγραμμένοι θάνατοι με εργαστηριακά επιβεβαίωμένη γρίπη"),
    bty="n", border=NA, fill=c("purple","orange"), xpd=TRUE, inset=c(0,-0.2))
  
  with(inflAttr_seasonal$All, plotCI(x=ps[1,], y=AllFlu, 
      li=AllFlu.lo, ui=AllFlu.hi, add=TRUE, cex=0.0001, col="magenta", lwd=2))
  # Annotate with number of deaths
  text(x=ps[2,], y=labConfDeaths+diff(ylim)/40, labConfDeaths, cex=0.9)
  text(x=ps[1,]+diff(ps[,1])/2.5, y=inflAttr_seasonal$All$AllFlu+diff(ylim)/40, 
    inflAttr_seasonal$All$AllFlu, cex=0.9)
  
  axis(1, at=(ps[2,]+ps[1,])/2, sprintf("%s-%02d", names(labConfDeaths),
    as.integer(names(labConfDeaths))%%100+1), lwd=0, font=2, cex.axis=1)
}



pie_ <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
    init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
    col = NULL, border = NULL, lty = NULL, main = NULL, ticks = TRUE, showValues = FALSE, ...) 
{
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop("'x' values must be positive.")
    if (is.null(labels)) 
        labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)
    vals <- x; vals[vals==0] <- NA
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    dev.hold()
    on.exit(dev.flush())
    plot.window(xlim, ylim, "", asp = 1)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("white", "lightblue", "mistyrose", "lightcyan", 
                "lavender", "cornsilk")
        else par("fg")
    if (!is.null(col)) 
        col <- rep_len(col, nx)
    if (!is.null(border)) 
        border <- rep_len(border, nx)
    if (!is.null(lty)) 
        lty <- rep_len(lty, nx)
    angle <- rep(angle, nx)
    if (!is.null(density)) 
        density <- rep_len(density, nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
        P <- t2xy(mean(x[i + 0:1]))
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            if (ticks) lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
            text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
                adj = ifelse(P$x < 0, 1, 0), ...)
            if (showValues) text(0.7 * P$x, 0.7 * P$y, vals[i])
        }
    }
    title(main = main, ...)
    invisible(NULL)
}


plotPropTypeDeaths <- function(family="Fira Sans") {
  par(family=family)
  tgtyear <- (totDeaths$yearweek[1] %/% 100) + as.integer((totDeaths$yearweek[1] %% 100)<40)
  
  cols <- c("B" = "orchid3", "A/H3N2" = "springgreen3", "A/H1N1" = "deepskyblue3", "A" = "turquoise", "gray")

  lbd <- table(totDeaths$flutypef)
  names(lbd) <- c("B", "A/H3N2", "A/H1N1", "A", "")
  lbd <- lbd[lbd>0]
  
  est <- unlist(inflAttr_seasonal$All[as.character(tgtyear), c("FluB","FluH3","FluH1")])
  names(est) <- c("B", "A/H3N2", "A/H1N1")

  par(mfrow=c(1,2), family=family, oma=c(0,0,2,3))
  pie_(est, showValues=TRUE, border="white", ticks=FALSE, col=cols[names(est)], cex=1.2,
    main="Εκτιμώμενοι θάνατοι\nαποδιδόμενοι στη γρίπη", radius=1, xpd=NA, mar=c(0,0,0,0))
  pie_(lbd, showValues=TRUE, border="white", ticks=FALSE, col=cols[names(lbd)], cex=1.2,
    main="Καταγεγραμμένοι θάνατοι με\nεργαστηριακά επιβεβαιωμένη γρίπη", radius=1, xpd=NA, mar=c(0,0,0,0))
}



plotMomoTypeDeaths <- function(family="Fira Sans") {
  momo$yearweek <- as.integer(gsub("-","",as.character(momo$wk2)))
  momopl <- subset(momo, yearweek>=tgtyear*100+40 & yearweek<=tgtyear*100+120)[,c("yearweek","nbc","Pnb","UPIb2", "UPIb4")]
  momopl$inflAttr <- inflAttr_weekly$All$AllFlu[match(momopl$yearweek, as.integer(rownames(inflAttr_weekly$All)))]
  par(mar=c(7,4,2,2), family=family)
  plot(0, type="n", bty="l", 
    ylim=range(pretty(range(momopl[,2:5]))),
    xlim=c(1,nrow(momopl)), 
    xaxt="n", ylab="Αριθμός θανάτων (από όλες τις αιτίες)", xlab=NA)
  mtext("Έτος - Αριθμός εβδομάδας", line=5, side=1)
  polygon(x=c(1:nrow(momopl),nrow(momopl):1),
    y=c(momopl$nbc-momopl$inflAttr, rev(momopl$nbc)),
    col="magenta", border=NA)
  points(y=momopl$UPIb4, x=1:nrow(momopl), col="yellow3", type="l", lwd=2)
  points(y=momopl$UPIb2, x=1:nrow(momopl), col="orange2", type="l", lwd=2)
  points(y=momopl$Pnb, x=1:nrow(momopl), col="firebrick2", type="l", lwd=2)
  points(y=momopl$nbc, x=1:nrow(momopl), col="steelblue4", type="l", lwd=2)
  axis(1, at=1:nrow(momopl), las=2,
    labels=paste0(substr(momopl$yearweek,1,4), "-", substr(momopl$yearweek,5,6)))
  legend("topright", "Αποδιδόμενη στη γρίπη θνησιμότητα", pch=15, bty="n", col="magenta", pt.cex=3)
}




cairo_pdf("output/fluDeaths_yearly.pdf", width=8, height=5)
#par(bg="#fafafa")
plotYrFluDeaths()
dev.off()


cairo_pdf("output/fluDeaths_prop.pdf", width=10, height=6)
# par(bg="#fafafa")
plotPropTypeDeaths(family="Fira sans")
dev.off()


cairo_pdf("output/fluDeaths_momo.pdf", width=10, height=6)
# par(bg="#fafafa")
plotMomoTypeDeaths(family="Fira Sans")
dev.off()

