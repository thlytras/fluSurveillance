# Εφαρμογή στατιστικής επεξεργασίας των δεδομένων του συστήματος sentinel
# v2.2 © 2015, Θοδωρής Λύτρας
# Βασισμένο σε κώδικα SPSS, του Πάνου Κατερέλου και Σταύρου Πατρινού
# Τελευταία αναθεώρηση: Φεβρουάριος 2015

# **** Read/set global options ****

require(foreign) # Απαιτείται το πακέτο foreign, για την ανάγνωση αρχείων SPSS και EpiData
require(plotrix) # Χρειάζεται για την εκτύπωση 95% CI στο διάγραμμα

# Αν τρέχει ως web, το locale χρειάζεται ρύθμιση...
if(!interactive()) Sys.setlocale(locale="el_GR.UTF-8")

path_input = "./data/"
path_output = "./output/"
opts <- list()
if (file.exists(paste(path_input, "options.txt", sep=""))) {
  tmp <- read.table(paste(path_input, "options.txt", sep=""))
  opts <- as.list(tmp[,2])
  names(opts) <-tmp[,1]
  rm(tmp)
}
# Set default option values, if not already set
if (is.null(opts$calculateOld)) opts$calculateOld <- 1
if (is.null(opts$noDelayed)) opts$noDelayed <- 0
if (is.null(opts$repLimDay)) opts$repLimDay <- 3

required_files_old <- c("IKA5.sav", "Ika5a.rec", "nomos_populations.sav", "sent08.rec", "sent12.rec", "sentinel_doctors_10.2005.sav", "sentKY.rec", "SENTNEWA.sav", "oldSentinel.R")
tmp <- c(paste(path_input, required_files_old[1:8], sep=""), required_files_old[9])
names(tmp) <- required_files_old
required_files_old <- tmp; rm(tmp)
required_files <- c("sent14.rec")
tmp <- paste(path_input, required_files, sep="")
names(tmp) <- required_files
required_files <- tmp; rm(tmp)
readerrmsg <- c("Σφάλμα κατά την ανάγνωση του αρχείου","!!! Πιθανά το αρχείο είναι εσφαλμένο ή κατεστραμμένο.\nΕλέγξτε οτι έχετε βάλει το σωστό αρχείο στον κατάλογο")

# ******** ΠΡΟΕΤΟΙΜΑΣΙΑ ********

if(!interactive()) sink(paste(path_output, "output.txt", sep=""))

cat("\nΚαλωσήλθατε στο πρόγραμμα ανάλυσης του συστήματος sentinel για τις γαστρεντερίτιδες.\n")
cat("v0.1 © 2015, Θοδωρής Λύτρας\n")
cat("Διαβάστε το αρχείο README.html για σύντομες οδηγίες χρήσεως.\n")

# Έλεγχος αν υπάρχουν όλα τα απαραίτητα
if (!file.exists(gsub("/$","",gsub("(\\s).","",path_input)))) stop("Ο κατάλογος εισόδου δεδομένων δεν υπάρχει!")
if (!file.exists(gsub("/$","",gsub("(\\s).","",path_output)))) stop("Ο κατάλογος εξόδου δεδομένων δεν υπάρχει!")
cat("\nΥπάρχουν όλα τα απαραίτητα αρχεία?... ")
if (FALSE %in% file.exists(required_files)) {
  cat("ΌΧΙ! Δε μπορώ να συνεχίσω...\n\n")
  cat("Συγκεκριμένα λείπουν τα αρχεία: ")
  cat(paste(names(required_files)[!file.exists(required_files)], collapse=", "))
  cat("\n")
} else { cat("ΝΑΙ\n") }


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

repeat {
  if(interactive()) {
    input<-readline(paste("\nΠοιάς χρονιάς το διάγραμμα θέλετε? [",curyear,"-",(curyear+1),"] ",sep=""))
  } else {
    input <- commandArgs(TRUE)[2]
    if (is.na(input) || input=="-") input <- ""
  }
  if(input=="") { tgtyear<-curyear; break }
  else {
    input<-strsplit(input,"-")[[1]]
    if (is.na(input[2])) suppressWarnings(input[2]<-as.integer(input[1])+1)
    if (is.na(input[1]) | input[1]=="") suppressWarnings(input[1]<-as.integer(input[2])-1)
    suppressWarnings(input<-as.integer(input[1:2]))
    if ((input[2]-input[1])!=1 | is.na(input[2]-input[1])) input=NA
    if (!is.na(input[1]) && input[1]>2004 && input[1]<=curyear) { tgtyear<-input[1]; break }
    if (interactive()) {
      cat("\nΕσφαλμένη εισαγωγή - ξαναπροσπαθήστε!\n")
    } else {
      stop("\nΕσφαλμένη εισαγωγή για την ζητούμενη χρονια για το διάγραμμα.\nΔιακοπή του script.\n")
    }
  }
}


formats<-c()
if (capabilities("cairo")) formats <- append(formats,c("ps","pdf","svg"))  # Η βιβλιοθήκη cairo δημιουργεί γραφήματα σε Postscript, PDF και SVG. Συνήθως είναι εγκατεστημένη μόνο στο Linux.
if (capabilities("png")) formats <- append(formats,"png")
if (capabilities("tiff")) formats <- append(formats,"tiff")
if (capabilities("jpeg")) formats <- append(formats,"jpg")
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

timer<-system.time({   # έναρξη χρονομέτρησης

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



cat("\nΑνάγνωση του αρχείου δηλώσεων...\n")

tryCatch({
  suppressWarnings( sentinelBig <- read.epiinfo(file(paste(path_input,"sent14.rec",sep=""), encoding="windows-1253"), lower.case.names=TRUE) )
}, error=function(err){
  stop(paste(readerrmsg[1],"sent14.rec",readerrmsg[2],path_input))
})
required_fields <- c("nom", "eid", "monada", "etos", "ebdo", "totvis", "gastot", "totdays", "hmekat", "arxebd", "telebd")
if (FALSE %in% (required_fields %in% names(sentinelBig))) stop("Το αρχείο sent14.rec δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")   # Ανίχνευση λαθών

dates_to_check <- with(sentinelBig, which((hmekat < arxebd+6) | (as.POSIXlt(arxebd)$wday!=1) | ((telebd-arxebd)>6) | ((telebd-arxebd)<0) | (is.na(arxebd) | is.na(telebd))))
if (length(dates_to_check)>0) {
  cat("\nΠΡΟΣΟΧΗ! Προβλήματα στις ημερομηνίες ορισμένων δηλώσεων.\nΠαρακαλώ ελέγξτε τις παρακάτω ημερομηνίες:\n")
  print(dates_to_check)
}

sentinelBig$astikot <- ifelse(sentinelBig$monada=="KEYG", 2, 1) # Αγροτικός πληθυσμός(2) αν Κέντρα Υγείας, ειδάλλως αστικός πληθυσμός (1).
sentinelBig$yearweek <- with(sentinelBig, etos*100 + ebdo)
sentinelBig$oldeid <- ifelse(sentinelBig$eid %in% c(2,5,8), 2, 1) # Ειδικότητα με το παλιό σύστημα. 1 = Παθολόγοι, 2 = Παιδίατροι.
sentinelBig$neweid <- c(1,2,1,3,3,3,1,2,1)[as.integer(sentinelBig$eid)] # Ειδικότητα με NEO σύστημα, όπου οι ιατροί των κέντρων υγείας αντιμετωπίζονται χωριστά (=3).

# Έξτρα χακιά για τους ιατρούς του ΙΚΑ Αμαρουσίου
ika_marousi_pa <- c("7a1228","7a1230","7a1229") # Οι κωδικοί των ιατρών του ΙΚΑ Αμαρουσίου (τους χειριζόμαστε διαφορετικά).
ika_marousi_pd <- c("8a1268") # Οι κωδικοί των ιατρών του ΙΚΑ Αμαρουσίου (τους χειριζόμαστε διαφορετικά).
sentinelBig$neweid[sentinelBig$codeiat %in% ika_marousi_pa] <- 4 # Παθολόγοι ΙΚΑ Αμαρουσίου
sentinelBig$neweid[sentinelBig$codeiat %in% ika_marousi_pd] <- 5 # Παιδίατροι ΙΚΑ Αμαρουσίου


sentinelBig<-merge(sentinelBig[,c(required_fields, "codeiat", "astikot", "yearweek", "oldeid", "neweid")],nomos_populations,by.x="nom", by.y="nomadil",all.x=TRUE)
sentinelBig<-subset(sentinelBig,!is.na(ast_p_nu))   # Έξω οι δηλώσεις που δε ξέρουμε το νομό τους
sentinelBig<-subset(sentinelBig,!is.na(totvis))   # Πετώ έξω όσα (από λάθος) έχουν missing επισκέψεις (δηλαδή δεν έχουν παρονομαστές). Για τους ιδιώτες αυτό έχει ήδη γίνει.
# Σβήνουμε ότι δε χρειάζεται
sentinelBig <- transform(sentinelBig, nomokat=NULL, d_diam=NULL, didia_po=NULL, nomadil2=NULL, arxebd=NULL, telebd=NULL)

# Βγάλε έξω καθυστερημένες δηλώσεις, εφ' όσον έχει ενεργοποιηθεί το σχετικό option
if (opts$noDelayed) sentinelBig <- subset(sentinelBig, isoweek(hmekat-as.integer(opts$repLimDay), "both_num") <= yearweek)

# Υπολογίζουμε την πληρότητα της δήλωσης (αναλυτικά ΜΟΝΟ για το νέο sentinel)
if(tgtweek>=201440) {
  plirotita_eidikotita <- with(subset(sentinelBig, yearweek==tgtweek & !is.na(yearweek)), table(factor(monada, levels=c("IDIO", "KEYG", "PEDY")), factor(oldeid, levels=1:2)))
  plirotita_nuts <- with(subset(sentinelBig, yearweek==tgtweek & !is.na(yearweek)), table(factor(monada, levels=c("IDIO", "KEYG", "PEDY")), factor(nuts, levels=1:4)))  
  msg_marousi <- c(
      if (sum(subset(sentinelBig, yearweek==tgtweek & !is.na(yearweek))$codeiat %in% ika_marousi_pa)>0) "Έχει δηλώσει παθολόγος, " else
	    "ΔΕΝ έχει δηλώσει παθολόγος, ",
      if (sum(subset(sentinelBig, yearweek==tgtweek & !is.na(yearweek))$codeiat %in% ika_marousi_pd)>0) "έχει δηλώσει παιδίατρος" else 
		      "ΔΕΝ έχει δηλώσει παιδίατρος")
}
doc_rep <- with(sentinelBig,table(yearweek,oldeid))
doc_rep_new <- with(sentinelBig,table(yearweek,neweid)) # Με χωριστά τους ιατρούς των ΚΥ


cat("\nΕξαγωγή ILI rate...\n")

aggr1 <- aggregate(sentinelBig[,c("etos", "ebdo", "ast_p_nu", "agr_p_nu", "as_p_nu1", "ag_p_nu1", "as_p_nu2", "ag_p_nu2")], by=list(astikot=sentinelBig$astikot, nuts=sentinelBig$nuts, yearweek=sentinelBig$yearweek), mean, na.rm=TRUE)

aggr2 <- aggregate(sentinelBig[,c("totvis", "gastot")], by=list(astikot=sentinelBig$astikot, nuts=sentinelBig$nuts, yearweek=sentinelBig$yearweek), sum, na.rm=TRUE)

aggr2$gastot[is.na(aggr2$gastot)] = 0
aggr2$gas <- (aggr2$gastot/aggr2$totvis)*1000
aggr2$gas[is.na(aggr2$gas)] = 0

aggr2$gas_w <- ifelse(aggr2$astikot==1,
  (aggr2$gastot/aggr2$totvis)*1000*aggr1$ast_p_nu,
  (aggr2$gastot/aggr2$totvis)*1000*aggr1$agr_p_nu)
aggr2$gas_w[is.na(aggr2$gas_w)] = 0

# Υπολογισμός variance (βάσει διωνυμικής κατανομής και large-sample: p(1-p)/N ) για καθένα από τα ζυγισμένα rate.
aggr2$gas_w_var <- ifelse(aggr2$astikot==1,
  ((aggr2$gastot/aggr2$totvis)*(1-aggr2$gastot/aggr2$totvis)/aggr2$totvis)*(1000*aggr1$ast_p_nu)^2,
  ((aggr2$gastot/aggr2$totvis)*(1-aggr2$gastot/aggr2$totvis)/aggr2$totvis)*(1000*aggr1$agr_p_nu)^2)
aggr2$gas_w_var[is.na(aggr2$gas_w_var)] = 0

aggr2$gas_w1 <- ifelse(aggr2$astikot==1,
  (aggr2$gastot/aggr2$totvis)*1000*aggr1$as_p_nu1,
  (aggr2$gastot/aggr2$totvis)*1000*aggr1$ag_p_nu1)
aggr2$gas_w1[is.na(aggr2$gas_w1)] = 0

aggr2$gas_w2 <- ifelse(aggr2$astikot==1,
  (aggr2$gastot/aggr2$totvis)*1000*aggr1$as_p_nu2,
  (aggr2$gastot/aggr2$totvis)*1000*aggr1$ag_p_nu2)
aggr2$gas_w2[is.na(aggr2$gas_w2)] = 0


aggr3 <- aggregate(aggr2[,c("gas_w","gastot","totvis","gas_w_var")],by=list(yearweek=aggr2$yearweek), sum, na.rm=TRUE)
aggr3 <- subset(aggr3, yearweek>200427 & yearweek<210000)

aggr3 <- merge(aggr3, subset(data.frame(yearweek=as.integer(rownames(doc_rep)), pa=doc_rep[,1], pd=doc_rep[,2]), yearweek>200427 & yearweek<210000), by.x="yearweek")

doc_rep_new <- doc_rep_new[as.character(aggr3$yearweek), , drop=FALSE]


rownames(aggr3) <- aggr3$yearweek
aggr3$yearweek <- NULL
aggr3$gas_w <- round(aggr3$gas_w,2)


ratechart <- aggr3[,1]
names(ratechart) <- rownames(aggr3)
ratechart_var <- aggr3[,4]
names(ratechart_var) <- rownames(aggr3)


# ******** ΕΞΑΓΩΓΗ ΑΠΟΤΕΛΕΣΜΑΤΩΝ ********


cat("\nΔημιουργία γραφημάτων...\n")


# Τα γραφήματα είναι υπό μορφή συναρτήσεων για να είναι επαναχρησιμοποιήσιμα μετά την εκτέλεση του script

# Γράφημα που δείχνει μία περίοδο γρίπης (εβδ. 40 έως εβδ. 20)
# NEA συνάρτηση, με δυνατότητα δεύτερου scaled άξονα y. (Βλέπε README.html για οδηγίες χρήσης).
sentinel_graph <- function(years, col=rainbow(length(years)), 
	yaxis2=NA, mult=1, ygrid=0, lty=rep(1,length(years)), lwd=rep(1,length(years)),
	ylab="Κρούσματα γαστρεντερίτιδας ανά 1000 επισκέψεις",
	ylab2=NA, ylab2rot=TRUE, ci=FALSE)
{
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
      if (ci[i]) plotCI(set[,i], uiw=1.96*sqrt(set_var[,i]), col=col[i], sfrac=0.005, cex=0.01, xpd=TRUE, add=TRUE)
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
      if (ci[i]) plotCI(set[,i], uiw=1.96*sqrt(set_var[,i]), col=col[i], sfrac=0.005, cex=0.01, xpd=TRUE, add=TRUE)
    })
    axis(4)
    # Οι γραμμές του αριστερού y άξονα θέλουμε να βρίσκονται σε πρώτο πλάνο.
    par(new=TRUE)
    plot(0, type="n", axes=FALSE, ylim=c(0,limrate), xlim=c(1,ifelse(maxwk==53,34,33)), ylab=NA, xlab=NA)
    sapply(i1, function(i){
      points(set[,i], type="o", col=col[i], lwd=2*lwd[i], pch=16, cex=lwd[i])
      if (ci[i]) plotCI(set[,i], uiw=1.96*sqrt(set_var[,i]), col=col[i], sfrac=0.005, cex=0.01, xpd=TRUE, add=TRUE)
    })
  } else {
    sapply(1:length(years), function(i){
      points(set[,i], type="o", col=col[i], lwd=2*lwd[i], pch=16, lty=lty[i], cex=lwd[i])
      if (ci[i]) plotCI(set[,i], uiw=1.96*sqrt(set_var[,i]), col=col[i], sfrac=0.005, cex=0.01, xpd=TRUE, add=TRUE)
    })
  }
  return()
}

# Γράφημα με το ILI rate κατά σύστημα (για μία μόνο χρονιά)
sentinelGraphBySystem <- function(year, ylab="Κρούσματα γαστρεντερίτιδας ανά 1000 επισκέψεις", ygrid=0, ci=FALSE, plot=rep(TRUE,5))
{
  astAggr2 <- aggr2
  astAggr2$gas_w <- ifelse(astAggr2$astikot==1,
    (astAggr2$gastot/astAggr2$totvis)*1000*aggr1$as_p_nu2,
    (astAggr2$gastot/astAggr2$totvis)*1000*aggr1$ag_p_nu2)
  astAggr2$gas_w[is.na(astAggr2$gas_w)] = 0
  astAggr2$gas_w_var <- ifelse(astAggr2$astikot==1,
    ((astAggr2$gastot/astAggr2$totvis)*(1-astAggr2$gastot/astAggr2$totvis)/astAggr2$totvis)*(1000*aggr1$as_p_nu2)^2,
    ((astAggr2$gastot/astAggr2$totvis)*(1-astAggr2$gastot/astAggr2$totvis)/astAggr2$totvis)*(1000*aggr1$ag_p_nu2)^2)
  astAggr2$gas_w_var[is.na(astAggr2$gas_w_var)] = 0
  astAggr3 <- aggregate(astAggr2[,c("gas_w","gastot","totvis","gas_w_var")],by=list(yearweek=astAggr2$yearweek, astikot=astAggr2$astikot), sum, na.rm=TRUE)
  astAggr3 <- subset(astAggr3, yearweek>200427 & yearweek<210000)
  spAggr1 <- aggregate(sentinelBig[,c("etos", "ebdo", "ast_p_nu", "agr_p_nu", "as_p_nu1", "ag_p_nu1", "as_p_nu2", "ag_p_nu2")], by=list(monada=sentinelBig$monada, astikot=sentinelBig$astikot, nuts=sentinelBig$nuts, yearweek=sentinelBig$yearweek), mean, na.rm=TRUE)
  spAggr2 <- aggregate(sentinelBig[,c("totvis", "gastot")], by=list(monada=sentinelBig$monada, astikot=sentinelBig$astikot, nuts=sentinelBig$nuts, yearweek=sentinelBig$yearweek), sum, na.rm=TRUE)
  spAggr2$gastot[is.na(spAggr2$gastot)] = 0
  spAggr2$gas <- (spAggr2$gastot/spAggr2$totvis)*1000
  spAggr2$gas[is.na(spAggr2$gas)] = 0
  spAggr2$gas_w <- ifelse(spAggr2$astikot==1,
    (spAggr2$gastot/spAggr2$totvis)*1000*spAggr1$as_p_nu2,
    (spAggr2$gastot/spAggr2$totvis)*1000*spAggr1$ag_p_nu2)
  spAggr2$gas_w[is.na(spAggr2$gas_w)] = 0
  spAggr2$gas_w_var <- ifelse(spAggr2$astikot==1,
    ((spAggr2$gastot/spAggr2$totvis)*(1-spAggr2$gastot/spAggr2$totvis)/spAggr2$totvis)*(1000*spAggr1$as_p_nu2)^2,
    ((spAggr2$gastot/spAggr2$totvis)*(1-spAggr2$gastot/spAggr2$totvis)/spAggr2$totvis)*(1000*spAggr1$ag_p_nu2)^2)
  spAggr2$gas_w_var[is.na(spAggr2$gas_w)] = 0
  spAggr3 <- aggregate(spAggr2[,c("gas_w","gastot","totvis","gas_w_var")],by=list(yearweek=spAggr2$yearweek, monada=spAggr2$monada), sum, na.rm=TRUE)
  spAggr3 <- subset(spAggr3, yearweek>201402)
  astAggr3 <- subset(astAggr3, yearweek>201402)

  if (length(ci)==1) ci <- rep(ci, 5)
  maxwk <- ifelse(sum(as.integer(isoweek(as.Date(paste(year,"-12-31",sep="")))==53))>0, 53, 52)
  ywk <- as.character(c((year*100+40):(year*100+maxwk),((year+1)*100+1):((year+1)*100+20)))
  astAggr3_1 <- subset(astAggr3, astikot==1); rownames(astAggr3_1) <- astAggr3_1$yearweek
  spAggr3_KEYG <- subset(spAggr3, monada=="KEYG"); rownames(spAggr3_KEYG) <- spAggr3_KEYG$yearweek
  spAggr3_PEDY <- subset(spAggr3, monada=="PEDY"); rownames(spAggr3_PEDY) <- spAggr3_PEDY$yearweek
  spAggr3_IDIO <- subset(spAggr3, monada=="IDIO"); rownames(spAggr3_IDIO) <- spAggr3_IDIO$yearweek
  set <- cbind(astAggr3_1[ywk, "gas_w"], spAggr3_KEYG[ywk, "gas_w"], spAggr3_PEDY[ywk, "gas_w"], spAggr3_IDIO[ywk, "gas_w"], aggr3[ywk,1])
  set_var <- cbind(astAggr3_1[ywk, "gas_w_var"], spAggr3_KEYG[ywk, "gas_w_var"], spAggr3_PEDY[ywk, "gas_w_var"], spAggr3_IDIO[ywk, "gas_w_var"], aggr3[ywk,4])
  limrate <- (max(sapply(1:ncol(set), function(i)max((set[,i] + 1.96*sqrt(set_var[,i])*ci[i])*plot[i], na.rm=TRUE))) %/% 10 + 2) * 10
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
  pal <- c("magenta", "dodgerblue3", "orange", "green", "red")
  ltypal <- c("dotted", "solid", "solid", "solid", "solid")
  jitf <- 0.10 * c(0, 0, 1, -1, 0)
  lwdpal <- rep(2,5)
  sapply(1:ncol(set), function(i){
    if (plot[i]) {
      points(1:(maxwk-19) + jitf[i]*ci[i], set[,i], type="l", col=pal[i], lty=ltypal[i], lwd=lwdpal[i])
      if (ci[i]==TRUE) plotCI(1:(maxwk-19) + jitf[i]*ci[i], set[,i], uiw=1.96*sqrt(set_var[,i]), col=pal[i], sfrac=0.005, pch=19, cex=0.6, add=TRUE)
    }
  })
  legend("topleft", c("Αστικός\nπληθυσμός", "ΚΥ", "ΠΕΔΥ", "Ιδιώτες", "Συνολικό")[plot], col=pal, lty=ltypal, lwd=lwdpal, pt.cex=0.6, pch=19, inset=0.03, bg="white", box.col="white")
  return()
}

# Γράφημα με το ILI rate κατά NUTS (για μία μόνο χρονιά)
sentinelGraphByNUTS <- function(year, ylab="Κρούσματα γαστρεντερίτιδας ανά 1000 επισκέψεις", ygrid=0, ci=FALSE, plot=rep(TRUE,5))
{
  spAggr2 <- aggr2
  spAggr2$gastot[is.na(spAggr2$gastot)] = 0
  spAggr2$gas <- (spAggr2$gastot/spAggr2$totvis)*1000
  spAggr2$gas[is.na(spAggr2$gas)] = 0
  spAggr2$gas_w <- ifelse(spAggr2$astikot==1,
    (spAggr2$gastot/spAggr2$totvis)*1000*aggr1$as_p_nu1,
    (spAggr2$gastot/spAggr2$totvis)*1000*aggr1$ag_p_nu1)
  spAggr2$gas_w[is.na(spAggr2$gas_w)] = 0
  spAggr2$gas_w_var <- ifelse(spAggr2$astikot==1,
    ((spAggr2$gastot/spAggr2$totvis)*(1-spAggr2$gastot/spAggr2$totvis)/spAggr2$totvis)*(1000*aggr1$as_p_nu1)^2,
    ((spAggr2$gastot/spAggr2$totvis)*(1-spAggr2$gastot/spAggr2$totvis)/spAggr2$totvis)*(1000*aggr1$ag_p_nu1)^2)
  spAggr2$gas_w_var[is.na(spAggr2$gas_w)] = 0
  spAggr3 <- aggregate(spAggr2[,c("gas_w","gastot","totvis","gas_w_var")],by=list(yearweek=spAggr2$yearweek, nuts=spAggr2$nuts), sum, na.rm=TRUE)
  spAggr3 <- subset(spAggr3, yearweek>201402)

  if (length(ci)==1) ci <- rep(ci, 5)
  maxwk <- ifelse(sum(as.integer(isoweek(as.Date(paste(year,"-12-31",sep="")))==53))>0, 53, 52)
  ywk <- as.character(c((year*100+40):(year*100+maxwk),((year+1)*100+1):((year+1)*100+20)))
  spAggr3_n1 <- subset(spAggr3, nuts==1); rownames(spAggr3_n1) <- spAggr3_n1$yearweek
  spAggr3_n2 <- subset(spAggr3, nuts==2); rownames(spAggr3_n2) <- spAggr3_n2$yearweek
  spAggr3_n3 <- subset(spAggr3, nuts==3); rownames(spAggr3_n3) <- spAggr3_n3$yearweek
  spAggr3_n4 <- subset(spAggr3, nuts==4); rownames(spAggr3_n4) <- spAggr3_n4$yearweek
  set <- cbind(spAggr3_n1[ywk, "gas_w"], spAggr3_n2[ywk, "gas_w"], spAggr3_n3[ywk, "gas_w"], spAggr3_n4[ywk,"gas_w"], aggr3[ywk,1])
  set_var <- cbind(spAggr3_n1[ywk, "gas_w_var"], spAggr3_n2[ywk, "gas_w_var"], spAggr3_n3[ywk, "gas_w_var"], spAggr3_n4[ywk, "gas_w_var"], aggr3[ywk,4])
  limrate <- (max(sapply(1:ncol(set), function(i)max((set[,i] + 1.96*sqrt(set_var[,i])*ci[i])*plot[i], na.rm=TRUE))) %/% 10 + 2) * 10
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
  pal <- c("magenta", "dodgerblue3", "orange", "green", "red")
  jitf <- 0.10 * c(0, 0, 1, -1, 0)
  pallwd <- c(rep(2, 4), 3)
  sapply(1:ncol(set), function(i){
    if (plot[i]) {
      points(1:(maxwk-19) + jitf[i]*ci[i],  set[,i], type="l", col=pal[i], lwd=pallwd)
      if (ci[i]==TRUE) plotCI(1:(maxwk-19) + jitf[i]*ci[i], set[,i], uiw=1.96*sqrt(set_var[,i]), col=pal[i], sfrac=0.005, pch=19, cex=0.6, add=TRUE)
    }
  })
  legend("topleft", c("NUTS1 - Βόρεια Ελλάδα", "NUTS2 - Κεντρική Ελλάδα", "NUTS3 - Αττική", "NUTS4 - Νησιά Αιγαίου & Κρήτη", "Συνολικό")[plot], col=pal, lwd=pallwd, pt.cex=0.6, pch=19, inset=0.03, bg="white", box.col="white")
  return()
}



if (graphtype=="ps" || graphtype=="pdf") {
  graph1 <- call(paste("cairo_", graphtype, sep=""), filename = paste(path_output, "sentinel_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""), width=10, height=6)
  graph3 <- call(paste("cairo_", graphtype, sep=""), filename = paste(path_output, "sentinel_bySystem_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""), width=10, height=6)
  graph4 <- call(paste("cairo_", graphtype, sep=""), filename = paste(path_output, "sentinel_byNUTS_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""), width=10, height=6)
  
} else if (graphtype=="svg") {
  graph1 <- call("svg", filename = paste(path_output,"sentinel_",tgtyear,"-",(tgtyear+1),".svg",sep=""), width=10, height=6)
  graph3 <- call("svg", filename = paste(path_output,"sentinel_bySystem_",tgtyear,"-",(tgtyear+1),".svg",sep=""), width=10, height=6)
  graph4 <- call("svg", filename = paste(path_output,"sentinel_byNUTS_",tgtyear,"-",(tgtyear+1),".svg",sep=""), width=10, height=6)
  
} else if (graphtype=="jpg") {
  graph1 <- call("jpeg", filename = paste(path_output,"sentinel_",tgtyear,"-",(tgtyear+1),".jpg",sep=""), width=2800, height=1680, res=288)
  graph3 <- call("jpeg", filename = paste(path_output,"sentinel_bySystem_",tgtyear,"-",(tgtyear+1),".jpg",sep=""), width=2800, height=1680, res=288)
  graph4 <- call("jpeg", filename = paste(path_output,"sentinel_byNUTS_",tgtyear,"-",(tgtyear+1),".jpg",sep=""), width=2800, height=1680, res=288)
} else if (graphtype=="png") {
  graph1 <- call("png", filename = paste(path_output,"sentinel_",tgtyear,"-",(tgtyear+1),".png",sep=""), width=2800, height=1680, res=288)
  graph3 <- call("png", filename = paste(path_output,"sentinel_bySystem_",tgtyear,"-",(tgtyear+1),".png",sep=""), width=2800, height=1680, res=288)
  graph4 <- call("png", filename = paste(path_output,"sentinel_byNUTS_",tgtyear,"-",(tgtyear+1),".png",sep=""), width=2800, height=1680, res=288)
} else if (graphtype=="tiff") {
  graph1 <- call("tiff", filename = paste(path_output,"sentinel_",tgtyear,"-",(tgtyear+1),".tif",sep=""), width=2800, height=1680, res=288, compression="lzw")
  graph3 <- call("tiff", filename = paste(path_output,"sentinel_bySystem_",tgtyear,"-",(tgtyear+1),".tif",sep=""), width=2800, height=1680, res=288, compression="lzw")
  graph4 <- call("tiff", filename = paste(path_output,"sentinel_byNUTS_",tgtyear,"-",(tgtyear+1),".tif",sep=""), width=2800, height=1680, res=288, compression="lzw")
}

if (is.na(graphtype)) {
  cat("ΔΕΝ δημιουργήθηκαν γραφήματα!\nΑυτή η εγκατάσταση του R δεν έχει δυνατότητα εγγραφής σε κανένα γραφικό file format...\n")
} else {
  eval(graph1)
  sentinel_graph(tgtyear, col="red3", lty=1, lwd=1.5)
  dev.off()
  
  eval(graph3)
  sentinelGraphBySystem(tgtyear, ci=TRUE)
  dev.off()

  eval(graph4)
  sentinelGraphByNUTS(tgtyear, ci=TRUE)
  dev.off()

}

# Μετατρέπω την Postscript σε Encapsulate Postscript με το utility ps2eps, αν είναι εγκατεστημένο (MONO για linux)
if (graphtype=="ps") system(paste("ps2eps -f ", path_output, "*.ps", sep=""), ignore.stderr=TRUE, wait=FALSE)


cat("\nΑποθήκευση ανάλυσης...\n")

write.csv2(aggr3, file = paste(path_output,"ratechart.csv",sep=""))

cat("\nΤέλος!\n")

})   # Τέλος χρονομέτρησης

cat(paste("\nΧρόνος ανάλυσης:",round(timer["elapsed"],1),"δευτερόλεπτα\n\n"))

if(!interactive()) sink()
