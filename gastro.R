 

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
                offset=log(vis), data=x, family="poisson", weights=wgt)
        }, silent=TRUE)
        if (class(err)=="try-error") return(NA)
        # If convergence warnings, try "bobyqa" optimizer
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gri ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(vis), data=x, family="poisson", weights=wgt,
                control=glmerControl(optimizer="bobyqa"))
          }
        # If still convergence warnings, set calc.derivs=FALSE
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gri ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(vis), data=x, family="poisson", weights=wgt,
                control=glmerControl(calc.derivs=FALSE))
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
              offset=log(totvis), data=x, family="poisson", weights=wgt)
          }, silent=TRUE)
          if (class(err)=="try-error") return(NA)
          attr(m, "fi") <- as.integer(levels(factor(x$grp)))
        } else {
          err <- try({
            m <- glmer(gritot ~ -1 + grp + (1|stratum) + (1|codeiat), 
              offset=log(totvis), data=x, family="poisson", weights=wgt)
          }, silent=TRUE)
          if (class(err)=="try-error") return(NA)
          # If convergence warnings, try "bobyqa" optimizer
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gritot ~ -1 + grp + (1|stratum) + (1|codeiat), 
              offset=log(totvis), data=x, family="poisson", weights=wgt,
              control=glmerControl(optimizer="bobyqa"))
          }
          # If still convergence warnings, set calc.derivs=FALSE
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gritot ~ -1 + grp + (1|stratum) + (1|codeiat), 
              offset=log(totvis), data=x, family="poisson", weights=wgt,
              control=glmerControl(calc.derivs=FALSE))
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


