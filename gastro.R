require(foreign) # Απαιτείται το πακέτο foreign, για την ανάγνωση αρχείων SPSS και EpiData
require(plotrix) # Χρειάζεται για την εκτύπωση 95% CI στο διάγραμμα
require(lme4)
require(splines)
require(pbs)

source("include.R")

path_input = "./data/"
path_output = "./output/"
readerrmsg <- c("Σφάλμα κατά την ανάγνωση του αρχείου","!!! Πιθανά το αρχείο είναι εσφαλμένο ή κατεστραμμένο.\nΕλέγξτε οτι έχετε βάλει το σωστό αρχείο στον κατάλογο")

cat("\nΚαλωσήλθατε στο πρόγραμμα ανάλυσης του συστήματος sentinel για τις γαστρεντερίτιδες.\n")
cat("v3.0 © 2016, Θοδωρής Λύτρας\n")
cat("Διαβάστε το αρχείο README.html για σύντομες οδηγίες χρήσεως.\n")

# ******** ΕΙΣΑΓΩΓΗ ΠΡΟΤΙΜΗΣΕΩΝ ********

curweek<-isoweek(Sys.Date()-7,"both_num")
curyear<-ifelse((curweek%%100)>=40, curweek%/%100, (curweek%/%100)-1)

repeat {
  if(interactive()) {
    input<-readline(paste("\nΓια ποιά εβδομάδα να δείξω αναλυτικά στοιχεία? (YYYYWW) [",curweek,"] ",sep=""))
  } else {
    input <- commandArgs(TRUE)[1]
    if (is.na(input) || input=="-") input <- ""
  }
  if(input=="") { tgtweek<-curweek; break }
  else {
    suppressWarnings(input<-as.integer(input))
    if (!is.na(input) && (input%%100<54) && input>200427 && input<=curweek) { tgtweek<-input; break }
    if (interactive()) {
      cat("\nΕσφαλμένη εισαγωγή - ξαναπροσπαθήστε!\n")
    } else {
      stop("\nΕσφαλμένη εισαγωγή για την ζητούμενη εβδομάδα.\nΔιακοπή του script.\n")
    }
  }
}


formats<-c()
if (capabilities("png")) formats <- append(formats,"png")
if (capabilities("tiff")) formats <- append(formats,"tiff")
if (capabilities("jpeg")) formats <- append(formats,"jpg")
if (capabilities("cairo")) formats <- append(formats,c("ps","pdf","svg"))  # Η βιβλιοθήκη cairo δημιουργεί γραφήματα σε Postscript, PDF και SVG. Συνήθως είναι εγκατεστημένη μόνο στο Linux.
if (length(formats)>0) {
  repeat {
    if(interactive()) {
      graphtype<-readline(paste("\nΣε τι μορφή τα θέλετε τα διαγράμματα? (",paste(formats,collapse=","),") [",formats[1],"] ",sep=""))
    } else {
      graphtype <- commandArgs(TRUE)[5]
      if (is.na(graphtype) || graphtype=="-") graphtype<-""
    }
    if (graphtype=="") graphtype<-formats[1]
    if (graphtype %in% formats) break
    if (interactive()) {
      cat("\nΕσφαλμένη εισαγωγή - ξαναπροσπαθήστε!\n")
    } else {
      stop("\nΕσφαλμένη εισαγωγή για τη μορφή των διαγραμμάτων.\nΔιακοπή του script.\n")
    }
  }
} else {
  cat("\nΣφάλμα! Δεν βρέθηκαν ρουτίνες εξαγωγής γραφικών!")
  cat("\nΔεν δύναμαι να εξάγω γραφήματα!")
  cat("\nΠαρακαλώ εγκαταστήστε τη βιβλιοθήκη cairo, ή μια νεότερη έκδοση του R.\n\n")
  graphtype=NA;
}



# ******** ΕΠΕΞΕΡΓΑΣΙΑ ΣΤΟΙΧΕΙΩΝ ********

#timer<-system.time({   # έναρξη χρονομέτρησης

cat("\nΕπεξεργασία στοιχείων - παρακαλώ περιμένετε...\n\n")


# Φόρτωση αρχείου πληθυσμών ανά νομό
tryCatch({
  nomos_populations <- read.spss(paste(path_input,"nomos_populations.sav",sep=""), reencode="windows-1253", to.data.frame=TRUE, use.value.labels=FALSE)
}, error=function(err){
  stop(paste(readerrmsg[1],"nomos_populations.sav",readerrmsg[2],path_input))
})
names(nomos_populations)<-tolower(names(nomos_populations))
if (FALSE %in% (c("nomadil", "nuts", "ast_p_nu", "agr_p_nu", "as_p_nu1", "ag_p_nu1", "as_p_nu2", "ag_p_nu2") %in% names(nomos_populations))) stop("Το αρχείο nomos_populations.sav δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")   # Ανίχνευση λαθών
nomos_populations<-nomos_populations[order(nomos_populations$nomadil),]


# Create table with total population proportion per stratum (NUTS:astikot)
NUTSpop <- data.frame(nuts=sort(rep(1:4,2)), astikot=rep(1:2,4))
NUTSpop$stratum <- with(NUTSpop, nuts*10 + astikot)
NUTSpop$prop <- mapply(function(n,a){   # prop = Stratum population / Total Greek population
    atxt <- c("ast_p_nu", "agr_p_nu")[a]
    subset(nomos_populations, nuts==n)[1,atxt]
},NUTSpop$nuts, NUTSpop$astikot)



cat("\nΑνάγνωση του αρχείου δηλώσεων...\n")

tryCatch({
  suppressWarnings( sentinelBig <- read.epiinfo(file(paste(path_input,"sent14.rec",sep=""), encoding="windows-1253"), lower.case.names=TRUE) )
}, error=function(err){
  stop(paste(readerrmsg[1],"sent14.rec",readerrmsg[2],path_input))
})
required_fields <- c("nom", "eid", "monada", "etos", "ebdo", "totvis", "gastot", "gritot", "totdays", "hmekat", "arxebd", "telebd")
if (FALSE %in% (required_fields %in% names(sentinelBig))) stop("Το αρχείο sent14.rec δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")   # Ανίχνευση λαθών

dates_to_check <- with(sentinelBig, which((hmekat < arxebd+6) | (as.POSIXlt(arxebd)$wday!=1) | ((telebd-arxebd)>6) | ((telebd-arxebd)<0) | (is.na(arxebd) | is.na(telebd))))
dates_to_check <- sentinelBig$aa[dates_to_check[dates_to_check!=nrow(sentinelBig)]]

if (length(dates_to_check)>0) {
  cat("\nΠΡΟΣΟΧΗ! Προβλήματα στις ημερομηνίες ορισμένων δηλώσεων.\nΠαρακαλώ ελέγξτε τις παρακάτω ημερομηνίες:\n")
  print(dates_to_check)
}

sentinelBig$astikot <- ifelse(sentinelBig$monada=="KEYG" | sentinelBig$monada=="PIAT", 2, 1) # Αγροτικός πληθυσμός(2) αν Κέντρα Υγείας, ειδάλλως αστικός πληθυσμός (1).
sentinelBig$yearweek <- with(sentinelBig, etos*100 + ebdo)
sentinelBig$neweid <- sentinelBig$eid
sentinelBig$eid <- ifelse(sentinelBig$neweid %in% c(2,5,8), 2, 1) # Ειδικότητα με το παλιό σύστημα. 1 = Παθολόγοι, 2 = Παιδίατροι.
sentinelBig$neweid2 <- c(3,1,2,1,3,3,3,1,2,1)[as.integer(sentinelBig$neweid)+1] # Ειδικότητα με NEO σύστημα, όπου οι ιατροί των κέντρων υγείας αντιμετωπίζονται χωριστά (=3).

# Έξτρα χακιά για τους ιατρούς του ΙΚΑ Αμαρουσίου
ika_marousi_pa <- c("7a1228","7a1230","7a1229") # Οι κωδικοί των ιατρών του ΙΚΑ Αμαρουσίου (τους χειριζόμαστε διαφορετικά).
ika_marousi_pd <- c("8a1268") # Οι κωδικοί των ιατρών του ΙΚΑ Αμαρουσίου (τους χειριζόμαστε διαφορετικά).
sentinelBig$neweid2[sentinelBig$codeiat %in% ika_marousi_pa] <- 4 # Παθολόγοι ΙΚΑ Αμαρουσίου
sentinelBig$neweid2[sentinelBig$codeiat %in% ika_marousi_pd] <- 5 # Παιδίατροι ΙΚΑ Αμαρουσίου


sentinelBig$nom <- factor(sentinelBig$nom, levels=nomos_populations$nomadil)
sentinelBig<-merge(sentinelBig[,c(required_fields, "codeiat", "astikot", "yearweek", "neweid", "neweid2")],nomos_populations,by.x="nom", by.y="nomadil",all.x=TRUE)
sentinelBig<-subset(sentinelBig,!is.na(ast_p_nu))   # Έξω οι δηλώσεις που δε ξέρουμε το νομό τους
sentinelBig<-subset(sentinelBig,!is.na(totvis))   # Πετώ έξω όσα (από λάθος) έχουν missing επισκέψεις (δηλαδή δεν έχουν παρονομαστές). Για τους ιδιώτες αυτό έχει ήδη γίνει.
# Σβήνουμε ότι δε χρειάζεται
sentinelBig <- transform(sentinelBig, nomokat=NULL, d_diam=NULL, didia_po=NULL, nomadil2=NULL, arxebd=NULL, telebd=NULL)

sentinelBig <- subset(sentinelBig, !is.na(totvis) & totvis>0)   # Discard if total visits zero or missing
sentinelBig <- subset(sentinelBig, yearweek >= 201440 & yearweek < 210000)


load(sprintf("%s/oldSentinel/big.RData", path_input))
cn <- names(sentinelBig)[names(sentinelBig) %in% names(big)] # Common names
sentinelBig <- rbind(big[,cn], sentinelBig[,cn])




fitGastroModel <- function(big, NUTSpop, verbose=FALSE, returnModels=FALSE) {
  # Match stratum information
    big$stratum <- with(big, paste(nuts, astikot, sep=""))
    big$prop <- NUTSpop$prop[match(big$stratum, NUTSpop$stratum)]
  # Adjust records (if missing ILI cases replace with 0)
    big$gastot[is.na(big$gastot)] <- 0   # If ILI cases missing, set to zero
    big <- subset(big, !is.na(totvis) & totvis>0)   # Discard if total visits zero or missing
  # Fit the models
    suppressWarnings({   # Don't echo glmer warnings
      aggrm <- lapply(sort(unique(big$yearweek)), function(w) {
        x <- subset(big, yearweek==w)
        x$prop2 <- x$prop / table(x$stratum)[x$stratum]
        if (verbose) cat(".")
        x$wgt <- x$prop2/sum(x$prop2)*nrow(x)
        err <- try({
          m <- glmer(gastot ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(totvis), data=x, family="poisson", weights=wgt)
        }, silent=TRUE)
        if (class(err)=="try-error") return(NA)
        # If convergence warnings, try "bobyqa" optimizer
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gastot ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(totvis), data=x, family="poisson", weights=wgt,
                control=glmerControl(optimizer="bobyqa"))
          }
        # If still convergence warnings, set calc.derivs=FALSE
          if (length(summary(m)$optinfo$conv$lme4) > 0) {
            m <- glmer(gastot ~ 1 + (1|stratum) + (1|codeiat), 
                offset=log(totvis), data=x, family="poisson", weights=wgt,
                control=glmerControl(calc.derivs=FALSE))
          }
        return(m)
        })
    })
    names(aggrm) <- sort(unique(big$yearweek))
  # Extract point estimate (**exp**, per 1000) and **log**SD
    res <- as.data.frame.matrix(t(sapply(aggrm, function(m) {
        if (class(m)=="glmerMod") {
            c(gas=unname(1000*exp(fixef(m))), log.gas.sd=unname(sqrt(vcov(m)[1])) )
        } else {
            c(gas=NA, log.gas.sd=NA)
        }
    })))
    res$yearweek <- sort(unique(big$yearweek))
    if (returnModels) { attr(res, "aggrm") <- aggrm }
    res
}


cat("\nΕξαγωγή rate γαστρεντεριτίδων...\n")


resGastro <- fitGastroModel(sentinelBig, NUTSpop, verbose=TRUE)[-1,]
resGastro$N <- with(resGastro, 1 / ((gas/1000) * log.gas.sd^2))
resGastro$week <- resGastro$yearweek %% 100
resGastro$t <- 1:nrow(resGastro)

resGastro$Y <- with(resGastro, (gas/1000)*N)
modelGastro <- glm(Y ~ ns(t, knots=t[week==1]) + pbs(week, df=4), offset=log(N), data=resGastro, family="quasipoisson")

Z <- 2
od <- max(1,sum(modelGastro$weights * modelGastro$residuals^2)/modelGastro$df.r)

resGastro$Pnb <- predict(modelGastro, type="response")
resGastro$stdp <- predict(modelGastro, se.fit=TRUE)$se.fit
resGastro$UPInb <- with(resGastro, (Pnb^(2/3)+ Z*((4/9)*(Pnb^(1/3))*(od+(stdp^2)*(Pnb)))^(1/2))^(3/2) )
resGastro$zscore <- with(resGastro, (Y^(2/3) - Pnb^(2/3)) / ((4/9)*(Pnb^(1/3))*(od+Pnb*(stdp^2)))^(1/2))

resGastro$fitted <- with(resGastro, Pnb/resGastro$N*1000)
resGastro$UPI <- with(resGastro, UPInb/resGastro$N*1000)

resGastro$ywtp <- NA
resGastro$ywtp[with(resGastro, (week %% 100) %in% c(1,26))] <- with(resGastro, yearweek[(week %% 100) %in% c(1,26)])


save(resGastro, file=sprintf("%s/resGastro.RData", path_output))



plot(resGastro$gas, type="l", col="brown", xlim=c(600,660), ylim=c(0,40), xaxt="n", ylab="Rate / 1000 επισκέψεις", xlab=NA, bty="l")
axis(1, at=with(resGastro, t[!is.na(ywtp)]), labels=with(resGastro, yearweek[!is.na(ywtp)]), las=2)

points(resGastro$fitted, type="l", col="green", lwd=3)
points(resGastro$UPI, type="l", col="red", lwd=2)



