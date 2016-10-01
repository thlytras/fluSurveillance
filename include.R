
# Function to print warnings within a function (even before the top-level function returns)
# This is used in run_labDeaths.R
showMeWarnings <- function(expr) 
{
    .list_of_warnings <- c()
    frame_number <- sys.nframe()
    withCallingHandlers(expr, warning = function(w) 
    {
      .list_of_warnings <<- c(.list_of_warnings, w$message)
      invokeRestart("muffleWarning")
    })
    if (length(.list_of_warnings)>0) cat("\nΠΡΟΣΟΧΗ:\n")
    cat(paste(.list_of_warnings, collapse=""))
    .list_of_warnings
}


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

# Συνάρτηση εύρεσης της ημ/νίας έναρξης μιας ISO εβδομάδας
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

# Συνάρτηση μετατροπής των ημερομηνιών του SPSS σε κατανοητή από το R μορφή
spssdate<-function(x){as.Date(ISOdate(1582,10,14)+x)}

# Συνάρτηση για recoding τιμών σε μια μεταβλητή
recode <- function(target,from,to) {
  if (!is.vector(target)) stop("Argument \"target\" is not a vector")
  if (!is.vector(from)) stop("Argument \"from\" is not a vector")
  if (!is.vector(to)) stop("Argument \"to\" is not a vector")
  if (length(from)!=length(to)) stop("Arguments \"from\" and \"to\" must have same length")
  for (i in 1:length(from)) target<-gsub(from[i],to[i],target)
  target
  }


# ΣΥΝΑΡΤΗΣΕΙΣ ΕΠΕΞΕΡΓΑΣΙΑΣ ΔΕΔΟΜΕΝΩΝ & ΑΛΓΟΡΙΘΜΩΝ ΑΝΑΛΥΣΗΣ

# Function to fit the main model (whole country weekly ILI rate)
# Data.frame big should have the following columns: gritit, totvis, codeiat, nuts, astikot
fitMainModel <- function(big, NUTSpop, verbose=FALSE, returnModels=FALSE) {
  # Match stratum information
    big$stratum <- with(big, paste(nuts, astikot, sep=""))
    big$prop <- NUTSpop$prop[match(big$stratum, NUTSpop$stratum)]
  # Adjust records (if missing ILI cases replace with 0)
    big$gritot[is.na(big$gritot)] <- 0   # If ILI cases missing, set to zero
    big <- subset(big, !is.na(totvis) & totvis>0)   # Discard if total visits zero or missing
  # Fit the models
    suppressWarnings({   # Don't echo glmer warnings
      aggrm <- lapply(sort(unique(big$yearweek)), function(w) {
        x <- subset(big, yearweek==w)
        x$prop2 <- x$prop / table(x$stratum)[x$stratum]
        if (verbose) cat(".")
        x$wgt <- x$prop2/sum(x$prop2)*nrow(x)
        m <- glmer(gritot ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(totvis), data=x, family="poisson", weights=wgt)
        # If convergence warnings, try "bobyqa" optimizer
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gritot ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(totvis), data=x, family="poisson", weights=wgt,
                control=glmerControl(optimizer="bobyqa"))
          }
        # If still convergence warnings, set calc.derivs=FALSE
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gritot ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(totvis), data=x, family="poisson", weights=wgt,
                control=glmerControl(calc.derivs=FALSE))
          }
        return(m)
        })
    })
    names(aggrm) <- sort(unique(big$yearweek))
  # Extract point estimate (**exp**, per 1000) and **log**SD
    res <- as.data.frame.matrix(t(sapply(aggrm, function(m) c(
                gri=unname(1000*exp(fixef(m))), 
                log.gri.sd=unname(sqrt(vcov(m)[1]))
    ))))
    res$yearweek <- sort(unique(big$yearweek))
    if (returnModels) { attr(res, "aggrm") <- aggrm }
    res
}


# Descriptives by week, whole country
aggrByWeek <- function(big) {
    aggr1a <- aggregate(big[,c("totvis", "gritot")], by=list(yearweek=big$yearweek), sum, na.rm=TRUE)
    aggr1b <- as.data.frame.matrix(with(big, table(yearweek, factor(eid, levels=1:2))))
    aggr1c <- as.data.frame.matrix(with(big, table(yearweek, factor(nuts, levels=1:4))))

    names(aggr1b) <- c("pa","pd")
    names(aggr1c) <- paste("nuts", 1:4, sep="")

    rownames(aggr1a) <- sort(unique(big$yearweek))
    aggr1a$yearweek <- sort(unique(big$yearweek))
    aggr1b$yearweek <- sort(unique(big$yearweek))
    aggr1c$yearweek <- sort(unique(big$yearweek))
    
    merge(merge(aggr1a, aggr1b), aggr1c)
}



# Γράφημα για διαχρονική τάση του rate
diax_graph <- function(years,col="darkred") {
  ratechart <- resAll$gri; names(ratechart) <- resAll$yearweek
  set <- ratechart[names(ratechart)>=years[1] & names(ratechart)<=years[2]]
  limrate <- (max(set,na.rm=TRUE)%/%20+1)*20
  labelsel <- c()
  labelsel[1] = ((years[1]%/%100)+1)*100+1
  labelsel[2] = (years[2]%/%100)*100+1
  if (years[1]%%100>26) labelsel[3]=((years[1]%/%100)+1)*100+26 else labelsel[3]=(years[1]%/%100)*100+26
  if (years[2]%%100>26) labelsel[4]=(years[2]%/%100)*100+26 else labelsel[4]=((years[2]%/%100)-1)*100+26
  labelsel2 <- seq(labelsel[1], labelsel[2], by=100)
  labelsel <- sort(c(seq(labelsel[1], labelsel[2], by=100), seq(labelsel[3], labelsel[4], by=100)))
  labelpos <- 1:length(set)
  names(labelpos) <- names(set)
  par(mar=c(7.1,4.1,2.1,2.1), mgp=c(5,1,0))
  plot(0, type="n", bty="n", axes=FALSE, xlab="Εβδομάδα")
  par(new=TRUE, mgp=c(3,1,0))
  plot(0, type="n", bty="l", xaxt="n", yaxt="n", ylim=c(0,limrate), xlim=c(0,length(set)+30), ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις", xlab=NA)
  abline(v=labelpos[as.character(labelsel2)], lty=3, col="lightgrey")
  abline(h=seq(0,limrate,by=20), lty=3, col="lightgrey")
  points(set, type="l", lwd=2, col=col)
  axis(1, at=labelpos[as.character(labelsel)], labels=paste(labelsel%/%100,sprintf("%02d",labelsel%%100),sep="-"), las=3)
  axis(2, at=seq(0,limrate,by=10), labels=seq(0,limrate,by=10))
}

# Γράφημα που δείχνει μία περίοδο γρίπης (εβδ. 40 έως εβδ. 20)
# NEA συνάρτηση, με δυνατότητα δεύτερου scaled άξονα y. (Βλέπε README.html για οδηγίες χρήσης).
sentinel_graph <- function(years, col=rainbow(length(years)), 
	yaxis2=NA, mult=1, ygrid=0, lty=rep(1,length(years)), lwd=rep(1,length(years)),
	ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις",
	ylab2=NA, ylab2rot=TRUE, ci=FALSE)
{
  ratechart <- resAll$gri; names(ratechart) <- resAll$yearweek
  ratechart_var <- resAll$log.gri.sd; names(ratechart_var) <- resAll$yearweek
  if(length(ci)==1) ci <- rep(ci, length(years))
  maxwk <- ifelse(sum(as.integer(isoweek(as.Date(paste(years,"-12-31",sep="")))==53))>0, 53, 52)
  set <- sapply(years, function(x){ratechart[as.character(c((x*100+40):(x*100+maxwk),((x+1)*100+1):((x+1)*100+20)))]})
  set_var <- sapply(years, function(x){ratechart_var[as.character(c((x*100+40):(x*100+maxwk),((x+1)*100+1):((x+1)*100+20)))]})
  limrate <- (max(set,na.rm=TRUE)%/%10+2)*10
  if(!is.na(yaxis2[1])) {
    i2 <- match(yaxis2, years)
    i1 <- (1:length(years))[-match(yaxis2,years)]
    maxes <- apply(set, 2, max, na.rm=TRUE)
    maxes[i2] <- maxes[i2]/mult
    limrate <- (max(maxes, na.rm=TRUE)%/%10+2)*10
    limrate2 <- limrate*mult
  }
  par(mar=c(5.1,4.1+grepl("\n",ylab),2.1,2.1+(2.5*!is.na(yaxis2[1]))+grepl("\n",ylab2)))
  plot(0, type="n", bty="l", xaxt="n", ylim=c(0,limrate), xlim=c(1,ifelse(maxwk==53,34,33)), ylab=ylab, xlab="Εβδομάδα")
  axis(1, at=1:(maxwk-19), labels=NA)
  mtext(c(40:maxwk,1:20), side=1, at=1:(maxwk-19), cex=0.7, line=0.5)
  if (!is.na(ygrid)) {
    if (ygrid==0)
      ygrid <- ifelse(limrate<120, 10, 20)
    abline(h=seq(0,limrate,by=ygrid), lty=3, col="lightgrey")
  }
  legend("topright", legend=paste(years,years+1,sep="-"), col=col[1:length(years)], lwd=2*lwd, pch=16, 
      box.col="white", box.lwd=10, bg="white", inset=0.01, pt.cex=lwd, lty=lty)
  if (!is.na(yaxis2[1])) {
    sapply(i1, function(i){
      points(set[,i], type="o", col=col[i], lwd=2, pch=16)
      if (ci[i]) {
        plotCI(set[,i], 
            ui=exp(log(set[,i]) + 1.96*sqrt(set_var[,i])),
            li=exp(log(set[,i]) - 1.96*sqrt(set_var[,i])), 
            col=col[i], sfrac=0.005, cex=0.01, xpd=TRUE, add=TRUE)
      }
    })
    par(new=TRUE)
    plot(0, type="n", bty="u", xaxt="n", yaxt="n", ylim=c(0,limrate2), xlim=c(1,ifelse(maxwk==53,34,33)), ylab=NA, xlab=NA)
    if (ylab2rot) {
      text(par("usr")[2] + 0.10*diff(par("usr")[1:2]*(1+0.15*grepl("\n",ylab2))), 
	  par("usr")[4]-diff(par("usr")[3:4])/2, 
	  srt = -90, labels=ylab2, xpd=TRUE)
    } else {
      mtext(ylab2, side=4, line=2.5+grepl("\n",ylab))
    }
    sapply(i2, function(i){
      points(set[,i], type="o", col=col[i], lwd=2*lwd[i], pch=16, lty=lty[i], cex=lwd[i])
      if (ci[i]) {
        plotCI(set[,i], 
            ui=exp(log(set[,i]) + 1.96*sqrt(set_var[,i])),
            li=exp(log(set[,i]) - 1.96*sqrt(set_var[,i])), 
            col=col[i], sfrac=0.005, cex=0.01, xpd=TRUE, add=TRUE)
      }
    })
    axis(4)
    # Οι γραμμές του αριστερού y άξονα θέλουμε να βρίσκονται σε πρώτο πλάνο.
    par(new=TRUE)
    plot(0, type="n", axes=FALSE, ylim=c(0,limrate), xlim=c(1,ifelse(maxwk==53,34,33)), ylab=NA, xlab=NA)
    sapply(i1, function(i){
      points(set[,i], type="o", col=col[i], lwd=2*lwd[i], pch=16, cex=lwd[i])
      if (ci[i]) {
        plotCI(set[,i], 
            ui=exp(log(set[,i]) + 1.96*sqrt(set_var[,i])),
            li=exp(log(set[,i]) - 1.96*sqrt(set_var[,i])), 
            col=col[i], sfrac=0.005, cex=0.01, xpd=TRUE, add=TRUE)
      }
    })
  } else {
    sapply(1:length(years), function(i){
      points(set[,i], type="o", col=col[i], lwd=2*lwd[i], pch=16, lty=lty[i], cex=lwd[i])
      if (ci[i]) {
        plotCI(set[,i], 
            ui=exp(log(set[,i]) + 1.96*sqrt(set_var[,i])),
            li=exp(log(set[,i]) - 1.96*sqrt(set_var[,i])), 
            col=col[i], sfrac=0.005, cex=0.01, xpd=TRUE, add=TRUE)
      }
    })
  }
  return()
}
