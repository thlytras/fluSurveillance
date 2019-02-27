
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

# Function to fit the main flu model (whole country weekly ILI rate)
# Data.frame big should have the following columns: gritit, totvis, codeiat, nuts, astikot
fitFluModel <- function(big, NUTSpop, gri="gritot", vis="totvis", verbose=FALSE, returnModels=FALSE) {
  # Set variables for analysis
    big$gri <- big[,gri]
    big$vis <- big[,vis]
  # Match stratum information
    big$stratum <- with(big, paste(nuts, astikot, sep=""))
    big$prop <- NUTSpop$prop[match(big$stratum, NUTSpop$stratum)]
  # Adjust records (if missing ILI cases replace with 0)
    big$gri[is.na(big$gri)] <- 0   # If ILI cases missing, set to zero
    big <- subset(big, !is.na(vis) & vis>0)   # Discard if total visits zero or missing
  # Fit the models
    suppressWarnings({   # Don't echo glmer warnings
      aggrm <- lapply(sort(unique(big$yearweek)), function(w) {
        x <- subset(big, yearweek==w)
        x$prop2 <- x$prop / table(x$stratum)[x$stratum]
        if (verbose) cat(".")
        x$wgt <- x$prop2/sum(x$prop2)*nrow(x)
        err <- try({
          m <- glmer(gri ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(vis), data=x, family="poisson", weights=wgt,
                control=glmerControl(optimizer="bobyqa", 
                    check.conv.singular=.makeCC(action="ignore", tol=1e-4)))
        }, silent=TRUE)
        if (class(err)=="try-error") return(NA)
        # If convergence warnings, set calc.derivs=FALSE
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gri ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(vis), data=x, family="poisson", weights=wgt,
                control=glmerControl(optimizer="bobyqa", calc.derivs=FALSE,
                    check.conv.singular=.makeCC(action="ignore", tol=1e-4)))
          }
        return(m)
        })
    })
    names(aggrm) <- sort(unique(big$yearweek))
  # Extract point estimate (**exp**, per 1000) and **log**SD
    res <- as.data.frame.matrix(t(sapply(aggrm, function(m) {
        if (class(m)=="glmerMod") {
            c(gri=unname(1000*exp(fixef(m))), log.gri.sd=unname(sqrt(vcov(m)[1])) )
        } else {
            c(gri=NA, log.gri.sd=NA)
        }
    })))
    res$yearweek <- sort(unique(big$yearweek))
    if (returnModels) { attr(res, "aggrm") <- aggrm }
    res
}


# Descriptives by week, whole country
aggrByWeek <- function(big, gri="gritot", vis="totvis") {
    aggr1a <- aggregate(big[,c(vis, gri)], by=list(yearweek=big$yearweek), sum, na.rm=TRUE)
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
diax_graph <- function(years, dataset="resAll", col="darkred", ci=FALSE, alpha=0.25) {
  ratechart <- get(dataset)$gri; names(ratechart) <- get(dataset)$yearweek
  set <- ratechart[names(ratechart)>=years[1] & names(ratechart)<=years[2]]
  limrate <- (max(set,na.rm=TRUE)%/%20+1)*20
  if (ci) {
    ratechart_logsd <- get(dataset)$log.gri.sd; names(ratechart_logsd) <- get(dataset)$yearweek
    set_logsd <- ratechart_logsd[names(ratechart_logsd)>=years[1] & names(ratechart_logsd)<=years[2]]
    limrate <- (max(exp(log(set) + 1.96*set_logsd),na.rm=TRUE)%/%20+1)*20
  }
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
  if (ci) {
    drawBand(x=1:length(set), col=addalpha(col, alpha),
        y.lo=exp(log(set) - 1.96*set_logsd),
        y.hi=exp(log(set) + 1.96*set_logsd))
  }
}



# Συνάρτηση για μετατροπή χρώματος σε διαφανές
# (χρησιμοποιείται στη συνάρτηση sentinel_graph() )
addalpha <- function(colors, alpha=1.0) {
    r <- col2rgb(colors, alpha=T)
    # Apply alpha
    r[4,] <- alpha*255
    r <- r/255.0
    return(rgb(r[1,], r[2,], r[3,], r[4,]))
  }

drawBand <- function(x, y.lo, y.hi, col) {
  if (length(x)!=length(y.lo) || length(x)!=length(y.hi)) {
    stop("Arguments x, y.lo, y.hi must have the same length")
  }
  for (i in 1:length(x)) {
    if (sum(is.na(c(y.lo[i:(i+1)], y.hi[i:(i+1)])))==0) {
      polygon(x=x[i+c(0,1,1,0)], y=c(y.lo[i:(i+1)], y.hi[(i+1):i]), col=col, border=NA)
    } else if (sum(is.na(c(y.lo[i], y.hi[i])))==0 && (i==1 || sum(is.na(c(y.lo[i-1], y.hi[i-1])))>0)) {
      polygon(x=x[c(i,i)], y=c(y.lo[i], y.hi[i]), col=col, border=col)
    }
  }
}


# Γράφημα που δείχνει μία περίοδο γρίπης (εβδ. 40 έως εβδ. 20)
sentinel_graph <- function(years, dataset="resAll", col=rainbow(length(years)), 
	yaxis2=NA, mult=1, ygrid=0, lty=rep(1,length(years)), lwd=rep(1,length(years)),
	ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις",
	ylab2=NA, ylab2rot=TRUE, ci=FALSE, alpha=0.1)
{
  drawCI <- function(i) {
    if (ci[i]) {
      if (alpha[i]==0) {
        plotCI(set[,i], 
            ui=exp(log(set[,i]) + 1.96*set_logsd[,i]),
            li=exp(log(set[,i]) - 1.96*set_logsd[,i]), 
            col=col[i], sfrac=0.005, cex=0.01, xpd=TRUE, add=TRUE)
      } else {
        drawBand(x=1:length(set[,i]), col=addalpha(col[i], alpha[i]),
            y.lo=exp(log(set[,i]) - 1.96*set_logsd[,i]),
            y.hi=exp(log(set[,i]) + 1.96*set_logsd[,i]))
      }
    }
  }
  ratechart <- get(dataset)$gri; names(ratechart) <- get(dataset)$yearweek
  ratechart_logsd <- get(dataset)$log.gri.sd; names(ratechart_logsd) <- get(dataset)$yearweek
  if(length(ci)==1) ci <- rep(ci, length(years))
  if(length(alpha)==1) alpha <- rep(alpha, length(years))
  maxwk <- ifelse(sum(as.integer(isoweek(as.Date(paste(years,"-12-31",sep="")))==53))>0, 53, 52)
  set <- sapply(years, function(x){ratechart[as.character(c((x*100+40):(x*100+maxwk),((x+1)*100+1):((x+1)*100+20)))]})
  set_logsd <- sapply(years, function(x){ratechart_logsd[as.character(c((x*100+40):(x*100+maxwk),((x+1)*100+1):((x+1)*100+20)))]})
  maxes <- sapply(1:ncol(set), function(i){
    suppressWarnings(max(exp(log(set[,i]) + ci[i]*1.96*set_logsd[,i]), na.rm=TRUE))
  })
  limrate <- (max(maxes, na.rm=TRUE)%/%10+2)*10
  if(!is.na(yaxis2[1])) {
    i2 <- match(yaxis2, years)
    i1 <- (1:length(years))[-match(yaxis2,years)]
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
      drawCI(i)
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
      drawCI(i)
    })
    axis(4)
    # Οι γραμμές του αριστερού y άξονα θέλουμε να βρίσκονται σε πρώτο πλάνο.
    par(new=TRUE)
    plot(0, type="n", axes=FALSE, ylim=c(0,limrate), xlim=c(1,ifelse(maxwk==53,34,33)), ylab=NA, xlab=NA)
    sapply(i1, function(i){
      points(set[,i], type="o", col=col[i], lwd=2*lwd[i], pch=16, cex=lwd[i])
      drawCI(i)
    })
  } else {
    sapply(1:length(years), function(i){
      points(set[,i], type="o", col=col[i], lwd=2*lwd[i], pch=16, lty=lty[i], cex=lwd[i])
      drawCI(i)
    })
  }
  return()
}




fitFluGroupModel <- function(grp, big, NUTSpop, verbose=FALSE, returnModels=FALSE) {
  # grp must be one of "asty", "nuts"
    if (!(grp %in% c("asty", "nuts"))) stop("Argument 'grp' must be one of \"asty\", \"nuts\".")
  # Match stratum information
    big$stratum <- with(big, paste(nuts, astikot, sep=""))
    big$prop <- NUTSpop$prop[match(big$stratum, NUTSpop$stratum)]
  # Adjust records (if missing ILI cases replace with 0)
    big$gritot[is.na(big$gritot)] <- 0   # If ILI cases missing, set to zero
    big <- subset(big, !is.na(totvis) & totvis>0)   # Discard if total visits zero or missing
    if (grp=="asty") {
      big$grp <- factor(big$astikot, levels=1:2)
      prop_key <- with(NUTSpop, tapply(prop, astikot, sum))
    } else {
      big$grp <- factor(big$nuts, levels=1:4)
      prop_key <- with(NUTSpop, tapply(prop, nuts, sum))
    }
  # Fit the models
    suppressWarnings({   # Don't echo glmer warnings
      aggrm <- lapply(sort(unique(big$yearweek)), function(w) {
        x <- subset(big, yearweek==w)
        x$prop2 <- x$prop / prop_key[x$grp] / table(x$stratum)[x$stratum]
        if (verbose) cat(".")
        x$wgt <- x$prop2 * table(x$grp)[x$grp]
        if (length(levels(factor(x$grp)))==1) {
          # If all but one grp levels are missing, fit an intercept-only model
          # (otherwise we would get an error....)
          err <- try({
            m <- glmer(gritot ~ 1 + (1|stratum) + (1|codeiat), 
              offset=log(totvis), data=x, family="poisson", weights=wgt,
              control=glmerControl(optimizer="bobyqa", 
                check.conv.singular=.makeCC(action="ignore", tol=1e-4)))
          }, silent=TRUE)
          if (class(err)=="try-error") return(NA)
          attr(m, "fi") <- as.integer(levels(factor(x$grp)))
        } else {
          err <- try({
            m <- glmer(gritot ~ -1 + grp + (1|stratum) + (1|codeiat), 
              offset=log(totvis), data=x, family="poisson", weights=wgt,
              control=glmerControl(optimizer="bobyqa", 
                check.conv.singular=.makeCC(action="ignore", tol=1e-4)))
          }, silent=TRUE)
          if (class(err)=="try-error") return(NA)
          # If convergence warnings, set calc.derivs=FALSE
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gritot ~ -1 + grp + (1|stratum) + (1|codeiat), 
              offset=log(totvis), data=x, family="poisson", weights=wgt,
              control=glmerControl(optimizer="bobyqa", calc.derivs=FALSE, 
                check.conv.singular=.makeCC(action="ignore", tol=1e-4)))
          }
        }
        return(m)
        })
    })
    names(aggrm) <- sort(unique(big$yearweek))
  # Extract point estimate (**exp**, per 1000) and **log**SD
    ncat <- c("asty"=2, "nuts"=4)[grp]
    res <- as.data.frame.matrix(t(sapply(aggrm, function(m) {
      if (class(m)!="glmerMod") {
        return(rep(NA, ncat*2))
      }
      if (length(fixef(m))==1) {
        return(c(rep(NA, (attr(m,"fi")-1)), 
            1000*exp(fixef(m)), rep(NA,ncat-1), sqrt(diag(vcov(m))), 
            rep(NA, (ncat-attr(m,"fi")))))
      } else {
        return(c(
            unname(1000*exp(fixef(m)[paste("grp",1:ncat,sep="")])), 
            unname(coef(summary(m))[,2][paste("grp",1:ncat,sep="")])
        ))
      }
    })))
    names(res) <- c(paste("gri", 1:ncat, sep=""), paste("log.gri.sd", 1:ncat, sep=""))
    res$yearweek <- sort(unique(big$yearweek))
    if (returnModels) { attr(res, "aggrm") <- aggrm }
    attr(res, "grp") <- grp
    res
}




# Γράφημα με το ILI rate ομαδοποιημένο κατά NUTS ή αστικότητα (για μία μόνο χρονιά)
sentinelGraphByGroup <- function(resGrp, year, ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις", ygrid=0, pal=c("magenta", "dodgerblue3", "orange", "green"), ci=FALSE, plot=TRUE, alpha=0.1)
{
  ncat <- (ncol(resGrp)-1)/2
  if (length(ci)==1) ci <- rep(ci, ncat)
  if (length(alpha)==1) alpha <- rep(alpha, ncat)
  if (length(plot)==1) plot <- rep(plot, ncat)
  maxwk <- ifelse(sum(as.integer(isoweek(as.Date(paste(year,"-12-31",sep="")))==53))>0, 53, 52)
  ywk <- as.character(c((year*100+40):(year*100+maxwk),((year+1)*100+1):((year+1)*100+20)))
  set <- resGrp[ywk,match(paste("gri",1:ncat,sep=""), names(resGrp))]
  set_logsd <- resGrp[ywk,match(paste("log.gri.sd",1:ncat,sep=""), names(resGrp))]
  limrate <- (max(sapply(1:ncol(set), function(i) suppressWarnings(max(exp(log(set[,i]) + 1.96*set_logsd[,i]*ci[i])*plot[i], na.rm=TRUE)))) %/% 10 + 2) * 10
  if (is.nan(limrate)) limrate <- 100
  par(mar=c(5.1,4.1+grepl("\n",ylab),2.1,2.1))
  plot(0, type="n", bty="l", xaxt="n", ylim=c(0,limrate), xlim=c(1,ifelse(maxwk==53,34,33)), ylab=ylab, xlab="Εβδομάδα")
  axis(1, at=1:(maxwk-19), labels=NA)
  mtext(c(40:maxwk,1:20), side=1, at=1:(maxwk-19), cex=0.7, line=0.5)
  if (!is.na(ygrid)) {
    if (ygrid==0)
      ygrid <- ifelse(limrate<120, 10, 20)
    abline(h=seq(0,limrate,by=ygrid), lty=3, col="lightgrey")
    abline(v=1:(maxwk-19), lty=3, col="lightgrey")
  }
  pal <- pal[1:ncat]
  jitf <- 0.10 * c(0, 0, 1, -1, 0)
  pallwd <- c(rep(2, 4), 3)
  sapply(1:ncol(set), function(i){
    if (plot[i]) {
      points(1:(maxwk-19) + jitf[i]*ci[i]*(alpha[i]==0),  set[,i], type="l", col=pal[i], lwd=pallwd)
      if (ci[i]==TRUE) {
        if (alpha[i]==0) {
          plotCI(1:(maxwk-19) + jitf[i]*ci[i], set[,i], 
            ui=exp(log(set[,i]) + 1.96*set_logsd[,i]),
            li=exp(log(set[,i]) - 1.96*set_logsd[,i]), 
            col=pal[i], sfrac=0.005, pch=19, cex=0.6, add=TRUE)
        } else {
          drawBand(x=1:(maxwk-19), col=addalpha(pal[i], alpha[i]),
            y.lo = exp(log(set[,i]) - 1.96*set_logsd[,i]),
            y.hi = exp(log(set[,i]) + 1.96*set_logsd[,i]))
        }
      }
    }
  })
  if (ncat==2) {
    legend("topleft", c("Αστικές", "Αγροτικές")[plot], col=pal, lwd=pallwd, pt.cex=0.6, pch=19, inset=0.03, bg="white", box.col="white")
  } else {
    legend("topleft", c("NUTS1 - Βόρεια Ελλάδα", "NUTS2 - Κεντρική Ελλάδα", "NUTS3 - Αττική", "NUTS4 - Νησιά Αιγαίου & Κρήτη")[plot], col=pal, lwd=pallwd, pt.cex=0.6, pch=19, inset=0.03, bg="white", box.col="white")
  }
  return()
}



fitGastroModel <- function(resG) {
  names(resG)[names(resG)=="gri"] <- "gas"
  names(resG)[names(resG)=="log.gri.sd"] <- "log.gas.sd"
  resG$N <- with(resG, 1 / ((gas/1000) * log.gas.sd^2))
  resG$week <- resG$yearweek %% 100
  resG$t <- 1:nrow(resG)
  resG$knots <- 0; resG$knots[resG$week==1] <- 1
  resG$knots[nrow(resG)-(9:0)] <- 0

  resG$Y <- with(resG, (gas/1000)*N)
  modelG <- glm(Y ~ ns(t, knots=which(knots==1)) + pbs(week, df=4), offset=log(N), data=resG, family="quasipoisson")

  Z <- 2
  od <- max(1,sum(modelG$weights * modelG$residuals^2)/modelG$df.r)

  resG$Pnb <- predict(modelG, type="response")
  resG$stdp <- predict(modelG, se.fit=TRUE)$se.fit
  resG$UPInb <- with(resG, (Pnb^(2/3)+ Z*((4/9)*(Pnb^(1/3))*(od+(stdp^2)*(Pnb)))^(1/2))^(3/2) )
  resG$zscore <- with(resG, (Y^(2/3) - Pnb^(2/3)) / ((4/9)*(Pnb^(1/3))*(od+Pnb*(stdp^2)))^(1/2))

  resG$fitted <- with(resG, Pnb/resG$N*1000)
  resG$UPI <- with(resG, UPInb/resG$N*1000)

  resG$ywtp <- FALSE
  resG$ywtp[with(resG, (week %% 100) %in% c(1,13,26,40))] <- TRUE
  
  attr(resG, "modelG") <- modelG
  attr(resG, "od") <- od
  
  resG
}




gastro_graph <- function(resGastro, back=120, legend=TRUE, ylab="Κρούσματα γαστρεντερίτιδας ανά 1000 επισκέψεις") {
  ymax <-  ceiling(max(c(resGastro$UPI, resGastro$gas), na.rm=TRUE)/10 + 1)*10
  plot(resGastro$gas, type="n", col="brown", xlim=ceiling(nrow(resGastro)/10 + 1)*10 + c(-back,0), ylim=c(0,ymax), xaxt="n", ylab=ylab, xlab=NA, bty="l", yaxs="i")
  abline(v=with(resGastro, t[ywtp]), col="grey", lty="dotted")
  abline(h=seq(0,ymax,by=10), col="grey", lty="dotted")
  axis(1, at=with(resGastro, t[ywtp]), labels=with(resGastro, yearweek[ywtp]), las=2)

  points(resGastro$fitted, type="l", col="green", lwd=2)
  points(resGastro$UPI, type="l", col="red", lwd=2, lty="dashed")
    
  points(resGastro$gas, type="l", col="purple", lwd=2)
  
  if (legend) {
    legend("top", legend=c("Παρατηρούμενα", "Αναμενόμενα", "+2SD"), 
        lty=c("solid","solid","dashed"), lwd=2, horiz=TRUE, 
        col=c("purple","green","red"), bg="white", box.col="white", seg.len=3)
  }
}
