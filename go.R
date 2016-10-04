fitRegionModel <- function(big, NUTSpop, verbose=FALSE, returnModels=FALSE) {
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
        x$prop2 <- x$prop / with(NUTSpop, tapply(prop, nuts, sum))[x$nuts] / table(x$nuts)[x$nuts]
        if (verbose) cat(".")
        x$wgt <- x$prop2/sum(x$prop2)*nrow(x)
        x$nutsf <- factor(x$nuts, levels=1:4)
        m <- glmer(gritot ~ -1 + nutsf + (1|stratum) + (1|codeiat), 
                offset=log(totvis), data=x, family="poisson", weights=wgt)
        # If convergence warnings, try "bobyqa" optimizer
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gritot ~ -1 + nutsf + (1|stratum) + (1|codeiat), 
                offset=log(totvis), data=x, family="poisson", weights=wgt,
                control=glmerControl(optimizer="bobyqa"))
          }
        # If still convergence warnings, set calc.derivs=FALSE
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gritot ~ -1 + nutsf + (1|stratum) + (1|codeiat), 
                offset=log(totvis), data=x, family="poisson", weights=wgt,
                control=glmerControl(calc.derivs=FALSE))
          }
        return(m)
        })
    })
    names(aggrm) <- sort(unique(big$yearweek))
  # Extract point estimate (**exp**, per 1000) and **log**SD
    res <- as.data.frame.matrix(t(sapply(aggrm, function(m) c(
                unname(1000*exp(fixef(m)[paste("nutsf",1:4,sep="")])), 
                unname(coef(summary(m))[,2][paste("nutsf",1:4,sep="")])
    ))))
    names(res) <- c(paste("gri", 1:4, sep=""), paste("log.gri.sd", 1:4, sep=""))
    res$yearweek <- sort(unique(big$yearweek))
    if (returnModels) { attr(res, "aggrm") <- aggrm }
    res
}




#resRegModel <- fitRegionModel(sentinelBig, NUTSpop, verbose=TRUE)
cat("\n\n")



# Γράφημα με το ILI rate κατά NUTS (για μία μόνο χρονιά)
sentinelGraphByNUTS <- function(year, ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις", ygrid=0, ci=FALSE, plot=rep(TRUE,4), alpha=0.1)
{
  if (length(ci)==1) ci <- rep(ci, 4)
  if (length(alpha)==1) alpha <- rep(alpha, 4)
  maxwk <- ifelse(sum(as.integer(isoweek(as.Date(paste(year,"-12-31",sep="")))==53))>0, 53, 52)
  ywk <- as.character(c((year*100+40):(year*100+maxwk),((year+1)*100+1):((year+1)*100+20)))
  
  set <- resRegModel[ywk,match(paste("gri",1:4,sep=""), names(resRegModel))]
  set_logsd <- resRegModel[ywk,match(paste("log.gri.sd",1:4,sep=""), names(resRegModel))]
  
  limrate <- (max(sapply(1:ncol(set), function(i)max(exp(log(set[,i]) + 1.96*set_logsd[,i]*ci[i])*plot[i], na.rm=TRUE))) %/% 10 + 2) * 10
  
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
  pal <- c("magenta", "dodgerblue3", "orange", "green")
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
          polygon(x=c(1:(maxwk-19), (maxwk-19):1),
            y=c(exp(log(set[,i]) - 1.96*set_logsd[,i]), rev(exp(log(set[,i]) + 1.96*set_logsd[,i]))),
            col=addalpha(pal[i], alpha[i]), border=addalpha(pal[i], alpha[i]))
        }
      }
    }
  })
  legend("topleft", c("NUTS1 - Βόρεια Ελλάδα", "NUTS2 - Κεντρική Ελλάδα", "NUTS3 - Αττική", "NUTS4 - Νησιά Αιγαίου & Κρήτη")[plot], col=pal, lwd=pallwd, pt.cex=0.6, pch=19, inset=0.03, bg="white", box.col="white")
  return()
}
 
sentinelGraphByNUTS(2015, ci=TRUE)
