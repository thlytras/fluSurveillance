library(gdata)
library(foreign)

path_input = "./data/"
path_output = "./output/"

# Συνάρτηση υπολογισμού της εβδομάδας κατά ISO
isoweek <- function(x, type="week") {
  alts=c("week","year","both_text","both_num")
  if(!(type %in% alts)) stop("Unknown isoweek type requested!")
  x.date<-as.Date(x)
  x.weekday<-as.integer(format(x.date,"%w"))
  x.weekday[x.weekday==0]=7
  x.nearest.thu<-x.date-x.weekday+4
  x.isoyear<-as.integer(substring(x.nearest.thu,1,4)) # Μπορεί οι πρώτες μέρες του χρόνου να ανήκουν (κατά ISO) στην προηγούμενη χρονιά!
  x.isoweek<-(as.integer(x.nearest.thu-as.Date(paste(x.isoyear,"-1-1",sep="")))%/%7)+1
  switch(type,
    week = x.isoweek,
    year = x.isoyear,
    both_text = ifelse((is.na(x.isoyear) | is.na(x.isoweek)),NA,paste(x.isoweek,x.isoyear,sep="/")),
    both_num = ifelse((is.na(x.isoyear) | is.na(x.isoweek)),NA,x.isoyear*100+x.isoweek)
    )
  }

isoweekStart <- function(x) {
  year <- x %/% 100
  week <- x %% 100
  x.date <- as.Date(paste(year,"-6-1", sep=""))
  x.weekday <- as.integer(format(x.date,"%w"))
  x.weekday[x.weekday==0]=7
  x.nearest.thu <- x.date-x.weekday+4
  x.isoweek <- isoweek(x.nearest.thu)
  res <- x.nearest.thu + 7*(week-x.isoweek) - 3
  if (sum(isoweek(res, type="both_num") != x)>0) stop("Error specifying ISO week number")
  return(res)
}

formatDate <- function(x) {
  x.year <- format(x, format="%Y")
  x.month <- as.integer(format(x, format="%m"))
  x.day <- format(x, format="%d")
  months <- c("Ιανουαρίου", "Φεβρουαρίου", "Μαρτίου",
      "Απριλίου", "Μαϊου", "Ιουνίου", "Ιουλίου", "Αυγούστου",
      "Σεπτεμβρίου", "Οκτωβρίου", "Νοεμβρίου", "Δεκεμβρίου")
  x.month <- months[x.month]
  if (length(x) != 2) return(paste(x.day, x.month, x.year))
  if ((x.year[1] == x.year[2]) & (x.month[1] == x.month[2]))
      return(paste(x.day[1], "–", x.day[2], " ", x.month[1], " ", x.year[1], sep=""))
  if (x.year[1] == x.year[2])
      return(paste(x.day[1], x.month[1], "–", x.day[2], x.month[2], x.year[1]))
  return(paste(x.day[1], x.month[1], x.year[1], "–", x.day[2], x.month[2], x.year[2]))
}


curweek<-isoweek(Sys.Date()-7,"both_num")

repeat {
  input<-readline(paste("\nΓια ποιά εβδομάδα να βγάλω τα διαγράμματα? (YYYYWW) [",curweek,"] ",sep=""))
  if(input=="") { tgtweek<-curweek; break }
  else {
    suppressWarnings(input<-as.integer(input))
    if (!is.na(input) && (input%%100<54) && input>200427 && input<=curweek) { tgtweek<-input; break }
    cat("\nΕσφαλμένη εισαγωγή - ξαναπροσπαθήστε!\n")
    }
  }

tgtyear <- tgtweek%/%100 - (tgtweek%%100 < 40)


formats<-c()
if (capabilities("png")) formats <- append(formats,"png")
if (capabilities("tiff")) formats <- append(formats,"tiff")
if (capabilities("jpeg")) formats <- append(formats,"jpg")
if (length(formats)>0) {
  repeat {
    graphtype<-readline(paste("\nΣε τι μορφή τα θέλετε τα διαγράμματα? (",paste(formats,collapse=","),") [",formats[1],"] ",sep=""))
    if (graphtype=="") graphtype<-formats[1]
    if (graphtype %in% formats) break
    cat("\nΕσφαλμένη εισαγωγή - ξαναπροσπαθήστε!\n")
    }
} else {
  cat("\nΣφάλμα! Δεν βρέθηκαν ρουτίνες εξαγωγής γραφικών!")
  cat("\nΔεν δύναμαι να εξάγω γραφήματα!")
  cat("\nΠαρακαλώ εγκαταστήστε τη βιβλιοθήκη cairo, ή μια νεότερη έκδοση του R.\n\n")
  stop("Αδύνατη η δημιουργία γραφημάτων.\n\n")
}

maxwk <- ifelse(sum(as.integer(isoweek(as.Date(paste(tgtyear,"-12-31",sep="")))==53))>0, 53, 52)
weekSel <- c(tgtyear*100+40:maxwk, (tgtyear+1)*100+1:20)

all_swabs <- Reduce(`+`, lapply(c("BOREIA ELLADA", "EKPA", "NOTIA.ELLADA"), function(x){
    a <- read.xls(paste(path_input, "GRIPI_", x, ".xls", sep=""), sheet=4)
    a <- sapply(a[3:(nrow(a)-1), -1], function(y)as.integer(as.character(y)))
    rownames(a) <- a[,1]
    colnames(a) <- c("Week", "Total", "Positive", "Negative", "Pending", 
			"Α", "Β", "C", "A.H3N2.", "A.H1N1.pdm09", "Other", "Pending.1")
    a[,-1]
}))
all_swabs[,"Α"] <- all_swabs[,"A.H3N2."] + all_swabs[,"A.H1N1.pdm09"] + all_swabs[,"Other"] + all_swabs[,"Pending.1"]
all_swabs[,"Positive"] <- all_swabs[,"Α"] + all_swabs[,"Β"] + all_swabs[,"C"]
all_swabs[,"Total"] <- all_swabs[,"Positive"] + all_swabs[,"Negative"] + all_swabs[,"Pending"]



meth <- read.xls(paste(path_input, "ICU_Full.xls", sep=""), sheet=2)
outOfMeth <- read.xls(paste(path_input, "ICU_Full.xls", sep=""), sheet=3)
names(meth)[c(8, 22, 13, 28, 29, 4, 5)] <- 
	c("flutype", "ICUadmDate", "hospAdmDate", "outcome", "deathdate", "sex", "age")
meth$yearweek <- isoweek(as.Date(as.character(meth$ICUadmDate), format="%d/%m/%y"), type="both_num")
meth$yearweekf <- factor(meth$yearweek, levels=weekSel)
meth$flutypef <- factor(meth$flutype, levels=c("B", "A(H3N2)", "A(H1N1)pdm09", "A", "χωρίς τύπο"))
meth$flutypef[is.na(meth$flutypef)] <- "χωρίς τύπο"

names(outOfMeth)[c(8, 13, 26, 27, 4, 5)] <- c("flutype", "hospAdmDate", "outcome", "deathdate", "sex", "age")
outOfMeth$flutypef <- factor(outOfMeth$flutype, levels=c("B", "A(H3N2)", "A(H1N1)pdm09", "A", "χωρίς τύπο"))
outOfMeth$flutypef[is.na(outOfMeth$flutypef)] <- "χωρίς τύπο"
outOfMeth$yearweek <- isoweek(as.Date(as.character(outOfMeth$hospAdmDate), format="%d/%m/%y"), type="both_num")
meth$meth <- TRUE; outOfMeth$meth <- FALSE
totDeaths <- rbind(subset(meth, outcome=="Θάνατος")[,c("flutypef", "deathdate", "meth", "sex", "age", "HRCG")], 
	subset(outOfMeth, outcome=="Θάνατος")[,c("flutypef", "deathdate", "meth", "sex", "age", "HRCG")])
totDeaths$yearweek <- isoweek(as.Date(as.character(totDeaths$deathdate), format="%d/%m/%y"), type="both_num")
totDeaths$yearweekf <- factor(totDeaths$yearweek, levels=weekSel)



momo <- read.dta(paste(path_input, "momoOutput.dta", sep=""))
momo$wy <- with(momo, paste(YoDi, formatC(WoDi, width=2, flag="0"), sep="-"))
momo$wy[!(momo$WoDi %in% c(1,26))] <- NA
momo <- subset(momo, nbc>0)


swabPlot <- function(limweek=tgtweek, ymax=NA){
    allSwabs <- all_swabs
    allSwabs[-(1:match(limweek, weekSel)),] <- 0
    if (is.na(ymax)) ymax <- max(allSwabs[,1])+20
    ymax <- ceiling(ymax/20)*20
    swCol <- c("dodgerblue3", "sandybrown", "red3", "orangered", "lightpink3", "darkgrey")
    barplot(t(allSwabs[,c("Β", "A.H3N2.", "A.H1N1.pdm09", "Other", "Pending.1", "Negative")]), 
	beside=FALSE, border=NA, col=swCol, las=2, axisnames=F, cex.axis=0.9, ylim=c(0, ymax),
	ylab="Φαρυγγικά δείγματα και απομονωθέντα στελέχη", font.lab=2, cex.lab=0.9)
    abline(h=seq(0,500,20), col="lightgrey", lwd=0.5)
    abline(h=0)
    bpos <- barplot(t(allSwabs[,c("Β", "A.H3N2.", "A.H1N1.pdm09", "Other", "Pending.1", "Negative")]), 
	beside=FALSE, border=NA, col=swCol, add=TRUE, axes=F, axisnames=F, ylim=c(0, ymax),
	legend.text=c("B", "A(H3N2)", "A(H1N1)pdm09", "A", "Αναμένεται τυποποίηση", "Αρνητικά δείγματα"),
	args.legend=list(bty="o", box.col="white", bg="white", border=NA, cex=0.8, x="topright"))
    axis(1, at=bpos[seq(1,length(bpos),2)], labels=rownames(allSwabs)[seq(1,length(bpos),2)], 
	lwd=0, cex.axis=0.8, line=-1)
    axis(1, at=bpos[seq(2,length(bpos),2)], labels=rownames(allSwabs)[seq(2,length(bpos),2)], 
	lwd=0, cex.axis=0.8, line=-1)
    axis(2, at=c(-10, 10000))
    mtext("Εβδομάδα", side=1, cex=0.9, font=2, line=1.5)
    return()
}

methPlot <- function(limweek=tgtweek){
    swCol <- c("dodgerblue3", "sandybrown", "red3", "orangered", "lightpink3", "darkgrey")
    bylim <- max(colSums(with(subset(meth, yearweek<=limweek), table(flutypef, yearweek))))
    bylim <- ifelse(bylim<=10, 20, 10+2*(bylim-10))
    bpos <- barplot(with(subset(meth, yearweek<=limweek), table(flutypef, yearweekf)), border=NA, col=c(swCol[1:4], "darkslategrey"), 
      axisnames=F, axes=F, ylab="Αριθμός κρουσμάτων", font.lab=2, cex.lab=0.9,
      ylim=c(0, bylim))
    abline(h=seq(0,bylim,2), col="lightgrey", lwd=0.5)
    abline(h=0)
    bpos <- barplot(with(subset(meth, yearweek<=limweek), table(flutypef, yearweekf)), border=NA, col=c(swCol[1:4], "darkslategrey"), 
      legend.text=TRUE, ylim=c(0,bylim), axisnames=F, axes=F, font.lab=2, add=TRUE,
      args.legend=list(bty="o", box.col="white", bg="white", border=NA, cex=0.8, x="topright", inset=c(0,-0.03)))
    axis(1, at=bpos[seq(1,length(bpos),2)], labels=(weekSel %% 100)[seq(1,length(bpos),2)], 
	lwd=0, cex.axis=0.8, line=-1)
    axis(1, at=bpos[seq(2,length(bpos),2)], labels=(weekSel %% 100)[seq(2,length(bpos),2)], 
	lwd=0, cex.axis=0.8, line=-1)
    axis(2, at=seq(0,bylim,2), las=2, cex.axis=0.9)
    mtext("Εβδομάδα εισαγωγής στη ΜΕΘ", side=1, cex=0.9, font=2, line=1.5)
    return()
}

deathPlot <- function(limweek=tgtweek){
    swCol <- c("dodgerblue3", "sandybrown", "red3", "orangered", "lightpink3", "darkgrey")
    bylim <- max(colSums(with(subset(totDeaths, yearweek<=limweek), table(flutypef, yearweek))))
    bylim <- ifelse(bylim<=10, 20, 10+2*(bylim-10))
    bpos <- barplot(with(subset(totDeaths, yearweek<=limweek), table(flutypef, yearweekf)), border=NA, col=c(swCol[1:4], "darkslategrey"), 
      axisnames=F, axes=F, ylab="Αριθμός κρουσμάτων", font.lab=2, cex.lab=0.9,
      ylim=c(0, bylim))
    abline(h=seq(0,bylim,2), col="lightgrey", lwd=0.5)
    abline(h=0)
    bpos <- barplot(with(subset(totDeaths, yearweek<=limweek), table(flutypef, yearweekf)), border=NA, col=c(swCol[1:4], "darkslategrey"), 
      legend.text=TRUE, ylim=c(0,bylim), axisnames=F, axes=F, font.lab=2, add=TRUE,
      args.legend=list(bty="o", box.col="white", bg="white", border=NA, cex=0.8, x="topright", inset=c(0,-0.03)))
    axis(1, at=bpos[seq(1,length(bpos),2)], labels=(weekSel %% 100)[seq(1,length(bpos),2)], 
	lwd=0, cex.axis=0.8, line=-1)
    axis(1, at=bpos[seq(2,length(bpos),2)], labels=(weekSel %% 100)[seq(2,length(bpos),2)], 
	lwd=0, cex.axis=0.8, line=-1)
    axis(2, at=seq(0,bylim,2), las=2, cex.axis=0.9)
    mtext("Εβδομάδα θανάτου", side=1, cex=0.9, font=2, line=1.5)
    return()
}

methDeathAgePlot <- function(limweek=tgtweek){
    a <- cbind("Εισαγωγές σε ΜΕΘ"=table(cut(subset(meth, yearweek<=limweek)$age, breaks=seq(0,100,10), right=FALSE, labels=paste(seq(0,90,10), seq(9,99,10), sep="-"))),
	"Θάνατοι"=table(cut(subset(totDeaths, yearweek<=limweek)$age, breaks=seq(0,100,10), right=FALSE, labels=paste(seq(0,90,10), seq(9,99,10), sep="-"))))
    barplot(t(a), border=NA, beside=TRUE, col=c("darkred", "darkgreen"), axisnames=F, axes=F, ylab=NA, font.lab=2, cex.lab=0.9, ylim=c(0,max(a)+6))
    abline(h=seq(0, max(a)+6, 2), col="lightgrey", lwd=0.5)
    abline(h=0)
    bpos <- barplot(t(a), border=NA, beside=TRUE, col=c("darkred", "darkgreen"), axisnames=F, axes=F, ylab="Αριθμός κρουσμάτων", font.lab=2, cex.lab=0.9, ylim=c(0,max(a)+6),
	add=TRUE, legend.text=TRUE, args.legend=list(bty="o", box.col="white", bg="white", border=NA, cex=0.8, x="topright", inset=c(0,-0.03)))
    axis(2, at=seq(0, max(a)+6, 2), las=2, cex.axis=0.9)

    axis(1, at=apply(bpos, 2, mean)[seq(1,ncol(bpos),2)], labels=rownames(a)[seq(1,ncol(bpos),2)], 
	lwd=0, cex.axis=0.8, line=-1)
    axis(1, at=apply(bpos, 2, mean)[seq(2,ncol(bpos),2)], labels=rownames(a)[seq(2,ncol(bpos),2)], 
	lwd=0, cex.axis=0.8, line=-1)

    mtext("Ηλικιακή ομάδα", side=1, cex=0.9, font=2, line=1.5)
    return(bpos)
}

momoPlot <- function() {
    par(mar=c(7,4,2,10))
    plot(0, type="n", bty="l", ylim=c(100,700), xlim=c(1,nrow(momo)), 
	xaxt="n", ylab="Αριθμός θανάτων", xlab=NA)
    mtext("Έτος - Αριθμός εβδομάδας", line=5, side=1)
    points(y=momo$UPIb4, x=1:nrow(momo), col="yellow3", type="l", lwd=2)
    points(y=momo$UPIb2, x=1:nrow(momo), col="orange2", type="l", lwd=2)
    points(y=momo$Pnb, x=1:nrow(momo), col="firebrick2", type="l", lwd=2)
    points(y=momo$nbc, x=1:nrow(momo), col="steelblue4", type="l", lwd=2)
    axis(1, at=which(!is.na(momo$wy)), labels=momo$wy[!is.na(momo$wy)], las=2)
    legend("topright", lwd=2, xpd=NA, bty="n", inset=c(-0.28,0.15), y.intersp=3, cex=0.8,
	col=c("yellow3", "orange2", "firebrick2", "steelblue4"),
	legend=c("+4 σταθερές αποκλίσεις\nαπό το αναμενόμενο",
	    "+4 σταθερές αποκλίσεις\nαπό το αναμενόμενο",
	    "Αναμενόμενος αριθμός\nθανάτων",
	    "Παρατηρούμενος αριθμός\nθανάτων"))
    return()
}


if (graphtype=="svg") {
  graph2 <- call("svg", filename = paste(path_output,"swabs.svg",sep=""), width=10, height=6, res=288)
  graph3 <- call("svg", filename = paste(path_output,"meth.svg",sep=""), width=10, height=6, res=288)
  graph3s <- call("svg", filename = paste(path_output,"meth-short.svg",sep=""), width=10, height=13*10/28, res=288)
  graph4 <- call("svg", filename = paste(path_output,"deaths.svg",sep=""), width=10, height=6, res=288)
  graph4s <- call("svg", filename = paste(path_output,"deaths-short.svg",sep=""), width=10, height=13*10/28, res=288)
  graph5 <- call("svg", filename = paste(path_output,"methDeathAges.svg",sep=""), width=10, height=6, res=288)
  graph6 <- call("svg", filename = paste(path_output,"momoPlot.svg",sep=""), width=10, height=6, res=288)
} else if (graphtype=="jpg") {
  graph2 <- call("jpeg", filename = paste(path_output,"swabs.jpg",sep=""), width=2800, height=1680, res=288)
  graph3 <- call("jpeg", filename = paste(path_output,"meth.jpg",sep=""), width=2800, height=1680, res=288)
  graph3s <- call("jpeg", filename = paste(path_output,"meth-short.jpg",sep=""), width=2800, height=1300, res=288)
  graph4 <- call("jpeg", filename = paste(path_output,"deaths.jpg",sep=""), width=2800, height=1680, res=288)
  graph4s <- call("jpeg", filename = paste(path_output,"deaths-short.jpg",sep=""), width=2800, height=1300, res=288)
  graph5 <- call("jpeg", filename = paste(path_output,"methDeathAges.jpg",sep=""), width=2800, height=1680, res=288)
  graph6 <- call("jpeg", filename = paste(path_output,"momoPlot.jpg",sep=""), width=2800, height=1680, res=288)
} else if (graphtype=="png") {
  graph2 <- call("png", filename = paste(path_output,"swabs.png",sep=""), width=2800, height=1680, res=288)
  graph3 <- call("png", filename = paste(path_output,"meth.png",sep=""), width=2800, height=1680, res=288)
  graph3s <- call("png", filename = paste(path_output,"meth-short.png",sep=""), width=2800, height=1300, res=288)
  graph4 <- call("png", filename = paste(path_output,"deaths.png",sep=""), width=2800, height=1680, res=288)
  graph4s <- call("png", filename = paste(path_output,"deaths-short.png",sep=""), width=2800, height=1300, res=288)
  graph5 <- call("png", filename = paste(path_output,"methDeathAges.png",sep=""), width=2800, height=1680, res=288)
  graph6 <- call("png", filename = paste(path_output,"momoPlot.png",sep=""), width=2800, height=1680, res=288)
} else if (graphtype=="tiff") {
  graph2 <- call("tiff", filename = paste(path_output,"swabs.tif",sep=""), width=2800, height=1680, res=288, compression="lzw")
  graph3 <- call("tiff", filename = paste(path_output,"meth.tif",sep=""), width=2800, height=1680, res=288, compression="lzw")
  graph3s <- call("tiff", filename = paste(path_output,"meth-short.tif",sep=""), width=2800, height=1300, res=288, compression="lzw")
  graph4 <- call("tiff", filename = paste(path_output,"deaths.tif",sep=""), width=2800, height=1680, res=288, compression="lzw")
  graph4s <- call("tiff", filename = paste(path_output,"deaths-short.tif",sep=""), width=2800, height=1300, res=288, compression="lzw")
  graph5 <- call("tiff", filename = paste(path_output,"methDeathAges.tif",sep=""), width=2800, height=1680, res=288, compression="lzw")
  graph6 <- call("tiff", filename = paste(path_output,"momoPlot.tif",sep=""), width=2800, height=1680, res=288, compression="lzw")
}

eval(graph2); swabPlot(); dev.off()
eval(graph3); methPlot(); dev.off()
eval(graph3s); par(mar=c(3.5,4,1,2)); methPlot(); dev.off()
eval(graph4); deathPlot(); dev.off()
eval(graph4s); par(mar=c(3.5,4,1,2)); deathPlot(); dev.off()
eval(graph5); methDeathAgePlot(); dev.off()
eval(graph6); momoPlot(); dev.off()



rb <- list()
rb$wklab <- paste(tgtweek %% 100, "/", tgtweek %/%100, " (", 
	formatDate(c(isoweekStart(tgtweek), isoweekStart(tgtweek)+6)), ")", sep="")
rb$seasEndLab <- paste(20, "/", tgtyear+1, " (", 
	formatDate(c(isoweekStart(tgtyear*100 + 120), isoweekStart(tgtyear*100 + 120) + 6)), ")", sep="")

rb$summSwab <- cbind(nosok = rowSums(mapply(function(x,i){
	a <- read.xls(paste(path_input, "GRIPI_", x, ".xls", sep=""), sheet=i)
	a <- as.integer(as.character(unlist(a[a[,(4-i/2)]==as.character(tgtweek%%100),])))[-(1:(4-i/2))]
	names(a) <- colnames(all_swabs)
	a
    }, c("BOREIA ELLADA", "EKPA", "NOTIA.ELLADA"), c(2,4,2)), na.rm=TRUE),
    sentinel = Reduce(`+`, lapply(c("BOREIA ELLADA", "NOTIA.ELLADA"), function(x){
	rowSums(sapply(c(1,3), function(i) {
	    a <- read.xls(paste(path_input, "GRIPI_", x, ".xls", sep=""), sheet=i)
	    a <- as.integer(as.character(unlist(a[a[,3]==as.character(tgtweek%%100),])))[-(1:3)]
	    names(a) <- colnames(all_swabs)
	    a
	}), na.rm=TRUE)
    })))
rb$summSwab["Α",] <- colSums(rb$summSwab[8:11,])
rb$summSwab["Positive",] <- colSums(rb$summSwab[5:7,])
rb$summSwab["Total",] <- colSums(rb$summSwab[2:4,])


rb$summSwabTot <- cbind(nosok = rowSums(mapply(function(x,i){
	a <- read.xls(paste(path_input, "GRIPI_", x, ".xls", sep=""), sheet=i)
	a <- as.integer(as.character(unlist(a[nrow(a),])))[-(1:(4-i/2))]
	names(a) <- colnames(all_swabs)
	a
    }, c("BOREIA ELLADA", "EKPA", "NOTIA.ELLADA"), c(2,4,2)), na.rm=TRUE),
    sentinel = Reduce(`+`, lapply(c("BOREIA ELLADA", "NOTIA.ELLADA"), function(x){
	rowSums(sapply(c(1,3), function(i) {
	    a <- read.xls(paste(path_input, "GRIPI_", x, ".xls", sep=""), sheet=i)
	    a <- as.integer(as.character(unlist(a[nrow(a),])))[-(1:3)]
	    names(a) <- colnames(all_swabs)
	    a
	}), na.rm=TRUE)
    })))
rb$summSwabTot["Α",] <- colSums(rb$summSwabTot[8:11,])
rb$summSwabTot["Positive",] <- colSums(rb$summSwabTot[5:7,])
rb$summSwabTot["Total",] <- colSums(rb$summSwabTot[2:4,])


rb$showPct <- function(x, y) { paste(x, " (", round(x*100/y, 1), "%)", sep="") }

rb$summF <- function(x,y="Total") rb$showPct(rowSums(rb$summSwab)[x], rowSums(rb$summSwab)[y])
rb$summFtot <- function(x,y="Total") rb$showPct(rowSums(rb$summSwabTot)[x], rowSums(rb$summSwabTot)[y])

rb$meth <- subset(meth, yearweek<=tgtweek)
rb$outOfMeth <- subset(outOfMeth, yearweek<=tgtweek)
rb$totDeaths <- subset(totDeaths, yearweek<=tgtweek)




cat("\nΈτοιμα τα διαγράμματα!")
cat("\nΜπορώ να ετοιμάσω και την εβδομαδιαία έκθεση.\n")
cat("Απαιτείται το πακέτο odfWeave, και να έχει τρέξει\n")
cat("(ή να είναι ήδη αποθηκευμένη) η ανάλυση του sentinel.\n")

repeat {
  input<-readline(paste("\nΝα ετοιμάσω την εβδομαδιαία έκθεση?\n (1=Ναι, 2=Όχι) [1] ",sep=""))
  if(input=="") { makeReport <- TRUE; break }
  else {
    suppressWarnings(input<-as.integer(input))
    if (is.na(input)) input <- 0
    if (input==1) { makeReport <- TRUE; break }
    if (input==2) { makeReport <- FALSE; break }
    cat("\nΕσφαλμένη εισαγωγή - ξαναπροσπαθήστε!\n")
    }
  }

if (makeReport) {
  cat("\nΕτοιμάζω την έκθεση...\n\n")
  library(odfWeave)
  if (!exists("sentinel_graph")) {
    tgtweek.bak <- tgtweek
    load(paste(path_output, "latest_analysis.RData", sep=""))
    tgtweek <- tgtweek.bak; rm(tgtweek.bak)
  }
  options("OutDec" = ",")
  #unlink(paste(path_output, "flureport-out.odt", sep=""))
  odfWeave(paste(path_input, "flureport.odt", sep=""), paste(path_output, "flureport-out.odt", sep=""))
  options("OutDec" = ".")
  dev.off()
  cat("\nΈτοιμη η αναφορά!\n")
}

cat("\nΟλοκληρώθηκε\n\n")

