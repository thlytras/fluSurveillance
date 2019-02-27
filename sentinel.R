# Εφαρμογή στατιστικής επεξεργασίας των δεδομένων του συστήματος sentinel
# v3.0 © 2016, Θοδωρής Λύτρας
# Βασισμένο σε κώδικα SPSS, του Πάνου Κατερέλου και Σταύρου Πατρινού
# Τελευταία αναθεώρηση: Οκτώβριος 2015

# **** Read/set global options ****

require(foreign) # Απαιτείται το πακέτο foreign, για την ανάγνωση αρχείων SPSS και EpiData
require(plotrix) # Χρειάζεται για την εκτύπωση 95% CI στο διάγραμμα
require(lme4)
require(splines)
require(pbs)
require(odfWeave)

source("include.R")

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
if (is.null(opts$weeksRecalc)) opts$weeksRecalc <- 40

required_files_old <- c("IKA5.sav", "Ika5a.rec", "sent08.rec", "sent12.rec", "sentinel_doctors_10.2005.sav", "sentKY.rec", "SENTNEWA.sav", "nomos_populations.sav", "oldSentinel.R")
tmp <- c(paste(path_input, "oldSentinel/", required_files_old[1:7], sep=""), paste(path_input, required_files_old[8], sep=""), required_files_old[9])
names(tmp) <- required_files_old
required_files_old <- tmp; rm(tmp)
required_files <- c("sent14.rec")
tmp <- paste(path_input, required_files, sep="")
names(tmp) <- required_files
required_files <- tmp; rm(tmp)
readerrmsg <- c("Σφάλμα κατά την ανάγνωση του αρχείου","!!! Πιθανά το αρχείο είναι εσφαλμένο ή κατεστραμμένο.\nΕλέγξτε οτι έχετε βάλει το σωστό αρχείο στον κατάλογο")

# ******** ΠΡΟΕΤΟΙΜΑΣΙΑ ********

if(!interactive()) sink(paste(path_output, "output.txt", sep=""))

cat("\nΚαλωσήλθατε στο πρόγραμμα ανάλυσης του συστήματος sentinel για τη γρίπη.\n")
cat("v2.1 © 2015, Θοδωρής Λύτρας\n")
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
if (opts$calculateOld) {
  opts$weeksRecalc <- 0
  cat("\nΖητήθηκε επανυπολογισμός του rate για το παλιό sentinel (μέχρι 2013-2014).")
  cat("\nΥπάρχουν όλα τα απαραίτητα αρχεία του παλιού sentinel?... ")
  if (FALSE %in% file.exists(required_files_old)) {
    cat("ΌΧΙ! Δε μπορώ να συνεχίσω...\n\n")
    cat("Συγκεκριμένα λείπουν τα αρχεία: ")
    cat(paste(names(required_files_old)[!file.exists(required_files_old)], collapse=", "))
    cat("\n")
  } else { cat("ΝΑΙ\n") }
} else {
  if (!file.exists(paste(path_input, "oldSentinel.RData", sep=""))) {
    stop("Δεν υπάρχει αρχείο αποτελεσμάτων του παλιού sentinel (oldSentinel.RData),\nκαι δε ζητήθηκε επανυπολογισμός του rate για το παλιό sentinel (μέχρι 2013-2014).\nΑδυνατώ να συνεχίσω...")
  }
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

repeat {
  if(interactive()) {
    input<-readline(paste("\nΘέλετε διαστήματα εμπιστοσύνης στις καμπύλες?\n (1=Όχι, 2=Ναι) [2] ",sep=""))
  } else {
    input <- commandArgs(TRUE)[3]
    if (is.na(input) || input=="-") input <- ""
  }
  if(input=="") { ciInPlot <- TRUE; break }
  else {
    suppressWarnings(input<-as.integer(input))
    if (is.na(input)) input <- 0
    if (input==1) { ciInPlot <- FALSE; break }
    if (input==2) { ciInPlot <- TRUE; break }
    if (interactive()) {
      cat("\nΕσφαλμένη εισαγωγή - ξαναπροσπαθήστε!\n")
    } else {
      stop("\nΕσφαλμένη εισαγωγή για το αν θα προσθέσω διαστήματα εμπιστοσύνης στις καμπύλες.\nΔιακοπή του script.\n")
    }
  }
}

repeat {
  if(interactive()) {
    input<-readline(paste("\nΑπό ποιά έως ποιά εβδομάδα να δείξω διαχρονικά την καμπύλη? (YYYYWW) [200426-",(curyear+1)*100+26,"] ",sep=""))
  } else {
    input <- commandArgs(TRUE)[4]
    if (is.na(input) || input=="-") input <- ""
  }
  if(input=="") { diaxyear<-c(200426,(curyear+1)*100+26); break }
  else {
    input<-strsplit(input,"-")[[1]]
    if (is.na(input[2])) input[2]<-(curyear+1)*100+26
    if (is.na(input[1]) | input[1]=="") input[1]<-200426
    suppressWarnings(input<-as.integer(input[1:2]))
    if (!is.na(input[1]) && !is.na(input[2]) && input[1]>=200426 && input[2]<=(curyear+2)*100+26 && input[2]>input[1] && input[1]%%100<54 && input[2]%%100<54) { diaxyear<-input; break }
    cat("\nΕσφαλμένη εισαγωγή - ξαναπροσπαθήστε!\n")
    if (interactive()) {
      cat("\nΕσφαλμένη εισαγωγή - ξαναπροσπαθήστε!\n")
    } else {
      stop("\nΕσφαλμένη εισαγωγή για το ζητούμενο διάστημα της διαχρονικής καμπύλης.\nΔιακοπή του script.\n")
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


# Create table with total population proportion per stratum (NUTS:astikot)
NUTSpop <- data.frame(nuts=sort(rep(1:4,2)), astikot=rep(1:2,4))
NUTSpop$stratum <- with(NUTSpop, nuts*10 + astikot)
NUTSpop$prop <- mapply(function(n,a){   # prop = Stratum population / Total Greek population
    atxt <- c("ast_p_nu", "agr_p_nu")[a]
    subset(nomos_populations, nuts==n)[1,atxt]
},NUTSpop$nuts, NUTSpop$astikot)
# Load population proportions by age group
load(paste(path_input, "NUTSpopAge.RData", sep=""))


if (opts$calculateOld) {
  source("oldSentinel.R", encoding="utf-8") # Υπολογίζει ξανά τα αποτελεσμάτα του παλιού sentinel...
} else { # ...ή τα φορτώνει΄έτοιμα, σωσμένα από προηγούμενη ανάλυση.
  load(paste(path_input, "oldSentinel.RData", sep=""))
}


cat("\nΑνάγνωση του αρχείου δηλώσεων...\n")

tryCatch({
  suppressWarnings( sentinelBig <- read.epiinfo(file(paste(path_input,"sent14.rec",sep=""), encoding="windows-1253"), lower.case.names=TRUE) )
}, error=function(err){
  stop(paste(readerrmsg[1],"sent14.rec",readerrmsg[2],path_input))
})
required_fields <- c("nom", "eid", "monada", "etos", "ebdo", "totvis", "gastot", "gritot", "totdays", "hmekat", "arxebd", "telebd", "gri1", "gri2", "gri3", "gri4", "gas1", "gas2", "gas3", "gas4", "vis1", "vis2", "vis3", "vis4")
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

# Βγάλε έξω καθυστερημένες δηλώσεις, εφ' όσον έχει ενεργοποιηθεί το σχετικό option
if (opts$noDelayed) sentinelBig <- subset(sentinelBig, isoweek(hmekat-as.integer(opts$repLimDay), "both_num") <= yearweek)

sentinelBig <- subset(sentinelBig, !is.na(totvis) & totvis>0)   # Discard if total visits zero or missing
sentinelBig <- subset(sentinelBig, yearweek >= 201440 & yearweek < 210000)

# Υπολογίζουμε την πληρότητα της δήλωσης
olemon <- c("IDIO"="IDIO", "KEYG"="KEYG", "PIAT"="KEYG", "PEDY"="PEDY") # Matrix to merge PIAT into KEDY
plirotites_eidikotita <- with(sentinelBig, table(
        systhma = factor(olemon[monada], levels=c("IDIO", "KEYG", "PEDY")), 
        eidikotita = factor(eid, levels=1:2, labels=c("PA","PD")), 
        yearweek))
plirotites_nuts <- with(sentinelBig, table(
        systhma = factor(olemon[monada], levels=c("IDIO", "KEYG", "PEDY")),
        nuts = factor(nuts, levels=1:4), 
        yearweek))
# Πολλαπλασιάζω τους παρατηρητές των ΚΥ επί 2 (να μη κακοκαρδίσουμε την Κάσσυ)
plirotites_eidikotita[2,,] <- plirotites_eidikotita[2,,]*2
plirotites_nuts[2,,] <- plirotites_nuts[2,,]*2
# Τσεκάρω το Μαρούσι...
if(tgtweek>=201440) {
  msg_marousi <- c(
      if (sum(subset(sentinelBig, yearweek==tgtweek & !is.na(yearweek))$codeiat %in% ika_marousi_pa)>0) "Έχει δηλώσει παθολόγος, " else
	    "ΔΕΝ έχει δηλώσει παθολόγος, ",
      if (sum(subset(sentinelBig, yearweek==tgtweek & !is.na(yearweek))$codeiat %in% ika_marousi_pd)>0) "έχει δηλώσει παιδίατρος" else 
		      "ΔΕΝ έχει δηλώσει παιδίατρος")
}


cat("\nΕξαγωγή ILI rate (βασικό μοντέλο)...\n")

limweek <- 201439
if (opts$weeksRecalc>0 & file.exists(paste(path_output, "res.RData", sep=""))) {
  limweek <- isoweek(isoweekStart(tgtweek)-opts$weeksRecalc*7, "both_num")
  if (limweek<201440) { limweek <- 201439 }
  load(paste(path_output, "res.RData", sep=""))
  limweek <- min(limweek, max(res$yearweek))
  resOld <- rbind(resOld, subset(res, yearweek<=limweek))
  resNutsOld <- rbind(resNutsOld, subset(resNuts, yearweek<=limweek))
  resAstyOld <- rbind(resAstyOld, subset(resAsty, yearweek<=limweek))
  resGastroOld <- rbind(resGastroOld, subset(resGastro, yearweek<=limweek))
  resGri1Old <- subset(resGri1, yearweek<=limweek)
  resGri2Old <- subset(resGri2, yearweek<=limweek)
  resGri3Old <- subset(resGri3, yearweek<=limweek)
  resGri4Old <- subset(resGri4, yearweek<=limweek)
  resGas1Old <- subset(resGas1, yearweek<=limweek)
  resGas2Old <- subset(resGas2, yearweek<=limweek)
  resGas3Old <- subset(resGas3, yearweek<=limweek)
  resGas4Old <- subset(resGas4, yearweek<=limweek)
}
resMainModel <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpop, verbose=TRUE)
descrByWeek <- aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek))
res <- merge(resMainModel, descrByWeek)

doc_rep_new <- with(sentinelBig,table(yearweek,neweid2)) # Με χωριστά τους ιατρούς των ΚΥ
doc_rep_new <- doc_rep_new[as.character(res$yearweek), , drop=FALSE]

# Υπολογισμός "εκτιμώμενου" συνολικού πληθυσμού
# (Πρώην excelάκι Κατερέλου-Καλαμάρα, βάσει του οποίου δηλώνουμε στο TESSy)
abcdland <- read.csv2(paste(path_input,"abcdland.csv",sep=""),header=FALSE)
res$popest <- round(
    (doc_rep_new[,"1"]/abcdland[1,2]*sum(abcdland[5:6,2])) # Παθολόγοι/γεν.ιατροί (πλην ΚΥ)
    + (doc_rep_new[,"2"]/abcdland[2,2]*sum(abcdland[3:4,2])) # Παιδίατροι (πλην ΚΥ)
    + 2*(doc_rep_new[,"3"]/sum(abcdland[1:2,2])*sum(abcdland[3:6,2])) # Ιατροί ΚΥ. Πολλ/ζονται επί 2.
    + 2*(doc_rep_new[,"4"]/abcdland[1,2]*sum(abcdland[5:6,2])) # Παθολόγοι ΙΚΑ Αμαρουσίου. Πολλ/ζονται επί 2.
    + 2*(doc_rep_new[,"5"]/abcdland[1,2]*sum(abcdland[3:4,2]))) # Παιδίατροι ΙΚΑ Αμαρουσίου. Πολλ/ζονται επί 2.

resAll <- rbind(resOld, res) # Συνένωση με τα αποτελέσματα του παλιού sentinel
res <- subset(resAll, yearweek>201439)


cat("\nΕξαγωγή rate γαστρεντερίτιδων (βασικό μοντέλο)...\n")
resGastroModel <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpop, "gastot", "totvis", verbose=TRUE)
descrGastroByWeek <- aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), "gastot", "totvis")
resGastro <- merge(resGastroModel, descrGastroByWeek)
resGastroAll <- rbind(resGastroOld, resGastro) # Συνένωση με τα αποτελέσματα του παλιού sentinel
resGastro <- subset(resGastroAll, yearweek>201439)


cat("\nΕξαγωγή ILI rate (κατά NUTS 1)...\n")
resNuts <- fitFluGroupModel("nuts", subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpop, verbose=TRUE)
resNutsAll <- rbind(resNutsOld, resNuts)
resNuts <- subset(resNutsAll, yearweek>201439)


cat("\nΕξαγωγή ILI rate (κατά αστικότητα)...\n")
resAsty <- fitFluGroupModel("asty", subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpop, verbose=TRUE)
resAstyAll <- rbind(resAstyOld, resAsty)
resAsty <- subset(resAstyAll, yearweek>201439)


cat("\nΕξαγωγή ILI rate (κατά ηλικιακές ομάδες)...\n")
resGri1 <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpopAge[[1]], "gri1", "vis1", verbose=TRUE)
resGri1 <- merge(resGri1, aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), "gri1", "vis1"))
if (exists("resGri1Old")) resGri1 <- rbind(resGri1Old, resGri1)
cat("\n")
resGri2 <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpopAge[[2]], "gri2", "vis2", verbose=TRUE)
resGri2 <- merge(resGri2, aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), "gri2", "vis2"))
if (exists("resGri2Old")) resGri2 <- rbind(resGri2Old, resGri2)
cat("\n")
resGri3 <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpopAge[[3]], "gri3", "vis3", verbose=TRUE)
resGri3 <- merge(resGri3, aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), "gri3", "vis3"))
if (exists("resGri3Old")) resGri3 <- rbind(resGri3Old, resGri3)
cat("\n")
resGri4 <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpopAge[[4]], "gri4", "vis4", verbose=TRUE)
resGri4 <- merge(resGri4, aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), "gri4", "vis4"))
if (exists("resGri4Old")) resGri4 <- rbind(resGri4Old, resGri4)


cat("\nΕξαγωγή rate γαστρεντεριτίδων (κατά ηλικιακές ομάδες)...\n")
resGas1 <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpopAge[[1]], "gas1", "vis1", verbose=TRUE)
resGas1 <- merge(resGas1, aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), "gas1", "vis1"))
if (exists("resGas1Old")) resGas1 <- rbind(resGas1Old, resGas1)
cat("\n")
resGas2 <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpopAge[[2]], "gas2", "vis2", verbose=TRUE)
resGas2 <- merge(resGas2, aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), "gas2", "vis2"))
if (exists("resGas2Old")) resGas2 <- rbind(resGas2Old, resGas2)
cat("\n")
resGas3 <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpopAge[[3]], "gas3", "vis3", verbose=TRUE)
resGas3 <- merge(resGas3, aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), "gas3", "vis3"))
if (exists("resGas3Old")) resGas3 <- rbind(resGas3Old, resGas3)
cat("\n")
resGas4 <- fitFluModel(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), NUTSpopAge[[4]], "gas4", "vis4", verbose=TRUE)
resGas4 <- merge(resGas4, aggrByWeek(subset(sentinelBig, yearweek>limweek & yearweek<=tgtweek), "gas4", "vis4"))
if (exists("resGas4Old")) resGas4 <- rbind(resGas4Old, resGas4)


cat("\nΠροσαρμογή μοντέλου γαστρεντεριτίδων...\n")
resGastroAll <- fitGastroModel(resGastroAll)
resGas1All <- fitGastroModel(resGas1)
resGas2All <- fitGastroModel(resGas2)
resGas3All <- fitGastroModel(resGas3)
resGas4All <- fitGastroModel(resGas4)


cat("\nΑποθήκευση ανάλυσης...\n")

write.csv2(resAll, file = paste(path_output,"ratechart.csv",sep=""))
write.csv2(resGastroAll, file = paste(path_output,"ratechart_gastro.csv",sep=""))
save(res, resNuts, resAsty, resGastro,
        resGri1, resGri2, resGri3, resGri4, 
        resGas1, resGas2, resGas3, resGas4, 
        file = paste(path_output,"res.RData",sep=""))


# ******** ΕΞΑΓΩΓΗ ΑΠΟΤΕΛΕΣΜΑΤΩΝ ********


cat("\nΔημιουργία γραφημάτων...\n")

# Τα γραφήματα είναι υπό μορφή συναρτήσεων για να είναι επαναχρησιμοποιήσιμα 
# μετά την εκτέλεση του script, και βρίσκονται στο αρχείο include.R

makeGraphCalls <- function(graphtype) {
  a1 <- c(ps="cairo_ps", pdf="cairo_pdf", svg="svg", jpg="jpeg", png="png", tiff="tiff")
  a2 <- c(
    paste(path_output, "sentinel_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_allyears.", graphtype, sep=""),
    paste(path_output, "sentinel_bySystem_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_byNUTS_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_AgeGr1_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_AgeGr2_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_AgeGr3_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_AgeGr4_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_Gastro_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_Gas1_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_Gas2_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_Gas3_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_Gas4_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_GasAgeGrAll_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep=""),
    paste(path_output, "sentinel_AgeGrAll_", tgtyear, "-", (tgtyear+1), ".", graphtype, sep="")
  )
  a3 <- c(rep(10, 3), rep(2800, 3)); names(a3) <- names(a1)
  a4 <- c(rep(6, 3), rep(1680, 3)); names(a4) <- names(a1)
  a5 <- c(rep("", 3), rep(", res=288", 3)); names(a5) <- names(a1)
  a6 <- c(rep("", 5), ", compression=lzw"); names(a6) <- names(a1)
  grcalls <- sapply(1:13, function(i) sprintf("%s(\"%s\", width=%s, height=%s%s%s)", a1[graphtype], 
        a2[i], a3[graphtype], a4[graphtype], a5[graphtype], a6[graphtype]) )
  grcalls[14] <- sprintf("%s(\"%s\", width=%s, height=%s%s%s)", a1[graphtype], a2[14], a3[graphtype], a3[graphtype], a5[graphtype], a6[graphtype]) 
  grcalls[15] <- sprintf("%s(\"%s\", width=%s, height=%s%s%s)", a1[graphtype], a2[15], a3[graphtype], a3[graphtype], a5[graphtype], a6[graphtype]) 
  grcalls
}

if (is.na(graphtype)) {
  cat("ΔΕΝ δημιουργήθηκαν γραφήματα!\nΑυτή η εγκατάσταση του R δεν έχει δυνατότητα εγγραφής σε κανένα γραφικό file format...\n")
} else {
  graphCalls <- makeGraphCalls(graphtype)
  eval(parse(text=graphCalls[1]))
  #ytp <- c(tgtyear-2, tgtyear-1,tgtyear)
  ytp <- c(tgtyear-1, tgtyear) # Προσωρινή τροποποίηση για φέτος (θα δείχνουμε μόνο 2 χρονιές)
  if (sum(ytp>=2014)>0 & sum(ytp>=2014)<length(ytp)) {
    sentinel_graph(ytp, col=c("black", "navyblue","red3"), lty=c(3,3,1), lwd=c(1,1,1.5), 
      ci=ciInPlot, alpha=c(0.1,0.15), yaxis2=ytp[ytp<2014], mult=1/5, 
      ylab=paste("Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις\n(Νέο σύστημα επιτήρησης, ",
        paste(paste(ytp[ytp>=2014],"-",ytp[ytp>=2014]+1,sep=""), collapse=", "), ")", sep=""),
      ylab2=paste("Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις\n(Παλιό σύστημα επιτήρησης, ",
        paste(paste(ytp[ytp<2014],"-",ytp[ytp<2014]+1,sep=""), collapse=", "), ")", sep=""))
  } else {
    #sentinel_graph(ytp, col=c("black", "navyblue","red3"), lty=c(3,3,1), lwd=c(1,1,1.5), ci=ciInPlot, alpha=c(0.1,0.15))
    sentinel_graph(ytp, col=c("navyblue","red3"), lty=c(3,1), lwd=c(1,1.5), ci=ciInPlot, alpha=c(0.1,0.15)) # Προσωρινή τροποποίηση για φέτος (θα δείχνουμε μόνο 2 χρονιές)
  }
  dev.off()
  
  eval(parse(text=graphCalls[2]))
  diax_graph(diaxyear, ci=ciInPlot, alpha=0.25)
  dev.off()
  
  eval(parse(text=graphCalls[3]))
  sentinelGraphByGroup(resAstyAll, tgtyear, ci=ciInPlot)
  dev.off()

  eval(parse(text=graphCalls[4]))
  sentinelGraphByGroup(resNutsAll, tgtyear, ci=ciInPlot)
  dev.off()

  eval(parse(text=graphCalls[5]))
  sentinel_graph(ytp, "resGri1", col=c("navyblue","red3"), lty=c(3,1), lwd=c(1,1.5), ci=ciInPlot, alpha=c(0.08,0.15), ylab=NA)
  mtext("Ηλικίες 0-4 ετών", side=3, cex=1.1)
  dev.off()

  eval(parse(text=graphCalls[6]))
  sentinel_graph(ytp, "resGri2", col=c("navyblue","red3"), lty=c(3,1), lwd=c(1,1.5), ci=ciInPlot, alpha=c(0.08,0.15), ylab=NA)
  mtext("Ηλικίες 5-14 ετών", side=3, cex=1.1)
  dev.off()

  eval(parse(text=graphCalls[7]))
  sentinel_graph(ytp, "resGri3", col=c("navyblue","red3"), lty=c(3,1), lwd=c(1,1.5), ci=ciInPlot, alpha=c(0.08,0.15), ylab=NA)
  mtext("Ηλικίες 15-64 ετών", side=3, cex=1.1)
  dev.off()

  eval(parse(text=graphCalls[8]))
  sentinel_graph(ytp, "resGri4", col=c("navyblue","red3"), lty=c(3,1), lwd=c(1,1.5), ci=ciInPlot, alpha=c(0.08,0.15), ylab=NA)
  mtext("Ηλικίες >=65 ετών", side=3, cex=1.1)
  dev.off()

  eval(parse(text=graphCalls[9]))
  gastro_graph(resGastroAll)
  dev.off()

  eval(parse(text=graphCalls[10]))
  gastro_graph(resGas1All)
  mtext("Ηλικίες 0-4 ετών", side=3, cex=1.1, line=1)
  dev.off()

  eval(parse(text=graphCalls[11]))
  gastro_graph(resGas2All)
  mtext("Ηλικίες 5-14 ετών", side=3, cex=1.1, line=1)
  dev.off()

  eval(parse(text=graphCalls[12]))
  gastro_graph(resGas3All)
  mtext("Ηλικίες 15-64 ετών", side=3, cex=1.1, line=1)
  dev.off()

  eval(parse(text=graphCalls[13]))
  gastro_graph(resGas4All)
  mtext("Ηλικίες >=65 ετών", side=3, cex=1.1, line=1)
  dev.off()

  eval(parse(text=graphCalls[14]))
  par(mfrow=c(4,1), oma=c(3,3,1,2))
  gastro_graph(resGas1All, ylab=NA)
  mtext("Ηλικίες 0-4 ετών", side=3, cex=1.1, line=1)
  gastro_graph(resGas2All, ylab=NA)
  mtext("Ηλικίες 5-14 ετών", side=3, cex=1.1, line=1)
  gastro_graph(resGas3All, ylab=NA)
  mtext("Ηλικίες 15-64 ετών", side=3, cex=1.1, line=1)
  gastro_graph(resGas4All, ylab=NA)
  mtext("Ηλικίες >=65 ετών", side=3, cex=1.1, line=1)
  mtext("Κρούσματα γαστρεντερίτιδας ανά 1000 επισκέψεις", side=2, outer=TRUE, line=-1, cex=1.1)
  dev.off()

  eval(parse(text=graphCalls[15]))
  par(mfrow=c(4,1), oma=c(3,3,1,2))
  sentinel_graph(ytp, "resGri1", col=c("navyblue","red3"), lty=c(3,1), lwd=c(1,1.5), ci=ciInPlot, alpha=c(0.08,0.15), ylab=NA)
  mtext("Ηλικίες 0-4 ετών", side=3, cex=1.1)
  sentinel_graph(ytp, "resGri2", col=c("navyblue","red3"), lty=c(3,1), lwd=c(1,1.5), ci=ciInPlot, alpha=c(0.08,0.15), ylab=NA)
  mtext("Ηλικίες 5-14 ετών", side=3, cex=1.1)
  sentinel_graph(ytp, "resGri3", col=c("navyblue","red3"), lty=c(3,1), lwd=c(1,1.5), ci=ciInPlot, alpha=c(0.08,0.15), ylab=NA)
  mtext("Ηλικίες 15-64 ετών", side=3, cex=1.1)
  sentinel_graph(ytp, "resGri4", col=c("navyblue","red3"), lty=c(3,1), lwd=c(1,1.5), ci=ciInPlot, alpha=c(0.08,0.15), ylab=NA)
  mtext("Ηλικίες >=65 ετών", side=3, cex=1.1)
  mtext("Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις", side=2, outer=TRUE, line=-1.3, cex=1.1)
  dev.off()

}

# Μετατρέπω την Postscript σε Encapsulate Postscript με το utility ps2eps, αν είναι εγκατεστημένο (MONO για linux)
if (graphtype=="ps") system(paste("ps2eps -f ", path_output, "*.ps", sep=""), ignore.stderr=TRUE, wait=FALSE)


# Επεξεργασία των πινάκων πληρότητας
if(tgtweek>=201440) {
  plirotita_eidikotita <- plirotites_eidikotita[,,match(as.character(tgtweek), dimnames(plirotites_eidikotita)$yearweek)]
  plirotita_nuts <- plirotites_nuts[,,match(as.character(tgtweek), dimnames(plirotites_nuts)$yearweek)]
} else {
  plirotita_eidikotita <- plirotites_eidikotita_old[,,match(as.character(tgtweek), dimnames(plirotites_eidikotita_old)$yearweek)]
  plirotita_nuts <- plirotites_nuts_old[,,match(as.character(tgtweek), dimnames(plirotites_nuts_old)$yearweek)]
}
# Υπολογίζουμε μερικά αθροίσματα
plirotita_eidikotita<-cbind(plirotita_eidikotita,rowSums(plirotita_eidikotita))
plirotita_nuts<-cbind(plirotita_nuts,rowSums(plirotita_nuts))
plirotita_eidikotita<-rbind(plirotita_eidikotita,colSums(plirotita_eidikotita))
plirotita_nuts<-rbind(plirotita_nuts,colSums(plirotita_nuts))
# Μορφοποίηση περιθωρίων
rownames(plirotita_eidikotita) <- c("Ιδιώτες", "ΚΥ", "ΙΚΑ", "Σύνολο")
colnames(plirotita_eidikotita) <- c("Παθολόγοι / Γεν.Ιατροί", "Παιδίατροι", "Σύνολο")
rownames(plirotita_nuts) <- c("Ιδιώτες", "ΚΥ", "ΙΚΑ", "Σύνολο")
colnames(plirotita_nuts) <- c("Βόρεια Ελλάδα", "Κεντρική Ελλάδα", "Αττική", "Νησιά Αιγαίου & Κρήτη", "Σύνολο")


# Πίνακες γριπωδών συνδρομών & γαστρεντεριτίδων
tbTotGri <- sapply(list(resGri1, resGri2, resGri3, resGri4, res), function(x) {
    a <- unlist(subset(x, yearweek==tgtweek)[,c(5:4,2)])
    if (length(a)==0) a <- c(gritot=0, totvis=0, gri=0)
    prettyNum(round(a, 2), big.mark=".", decimal.mark=",")
})
colnames(tbTotGri) <- c("0-4 ετών", "5-14 ετών", "15-65 ετών", "&gt;65 ετών", "Σύνολο")
rownames(tbTotGri) <- c("Αρ.γριπωδών συνδρομών", "Αρ.επισκέψεων", "Γριπώδεις συνδρομές / 1000 επισκέψεις")

tbTotGas <- sapply(list(resGas1, resGas2, resGas3, resGas4, resGastro), function(x) {
    a <- unlist(subset(x, yearweek==tgtweek)[,c(5:4,2)])
    if (length(a)==0) a <- c(gastot=0, totvis=0, gri=0)
    prettyNum(round(a, 2), big.mark=".", decimal.mark=",")
})
colnames(tbTotGas) <- c("0-4 ετών", "5-14 ετών", "15-65 ετών", "&gt;65 ετών", "Σύνολο")
rownames(tbTotGas) <- c("Αρ.γαστρεντεριτίδων", "Αρ.επισκέψεων", "Γαστρεντερίτιδες / 1000 επισκέψεις")



cat("\nΔημιουργία εκθέσεων sentinel...\n")

styleDefs <- getStyleDefs()

styleDefs$RTable0 <- styleDefs$RTable1
styleDefs$RTable0$marginLeft <- "0"
styleDefs$RTable0$marginRight <- "0"
styleDefs$RTable0$marginTop <- "0"
styleDefs$RTable0$marginBottom <- "0"

styleDefs$MyNormal <- styleDefs$ArialNormal
styleDefs$MyNormal$fontName <- "Fira Sans"
styleDefs$MyNormal$fontSize <- "10pt"

styleDefs$MyCenteredNormal <- styleDefs$MyNormal
styleDefs$MyCenteredNormal$textAlign <- "center"

styleDefs$MyBold <- styleDefs$MyCenteredNormal
styleDefs$MyBold$fontType <- "bold"

setStyleDefs(styleDefs)

styles <- getStyles()
styles$paragraph <- "MyCenteredNormal"
styles$cellText <- "MyCenteredNormal"
styles$headerText <- "MyCenteredNormal"
styles$table <- "RTable0"
setStyles(styles)

options("OutDec"=",")
odfWeave(sprintf("%ssentFluReport.odt", path_input), sprintf("%ssentFluReport-%s.odt", path_output, tgtweek))
odfWeave(sprintf("%ssentGastroReport.odt", path_input), sprintf("%ssentGastroReport-%s.odt", path_output, tgtweek))
options("OutDec"=".")
dev.off()



cat("\nΑποθήκευση συνόλου ανάλυσης...\n")

write.csv2(plirotita_nuts,paste(path_output,"plirotita_nuts_",tgtweek,".csv",sep=""))
write.csv2(plirotita_eidikotita,paste(path_output,"plirotita_eidikotita_",tgtweek,".csv",sep=""))

save.image(file = paste(path_output,"latest_analysis.RData",sep=""))


# Απεικόνιση πληρότητας
cat(paste("\nΙατροί που δήλωσαν κατά την εβδομάδα ",(tgtweek%%100),"/",((tgtweek%/%100)),"\n",sep=""))
print(plirotita_nuts)
print(plirotita_eidikotita)
cat("ΙΚΑ Αμαρουσίου: "); cat(msg_marousi[1]); cat(msg_marousi[2]); cat("\n")


showgri<-function(yweek) {
  result <- t(subset(resAll, yearweek==yweek)[,c("gritot","totvis")])
  rownames(result) <- c("Σύνολο γριπωδών συνδρομών","Σύνολο επισκέψεων")
  if (ncol(result)>1) colnames(result) <- "#"
  result
  }
print(showgri(tgtweek))
cat(paste("Rate για την εβδομάδα ",(tgtweek%%100),"/",((tgtweek%/%100))," :   ",prettyNum(round(subset(resAll, yearweek==tgtweek)$gri,2), decimal.mark=",")," γριπώδεις συνδρομές ανά 1000 επισκέψεις\n",sep=""))

# Απεικόνιση "εκτιμώμενου" συνολικού πληθυσμού
# (Excelάκι Κατερέλου-Καλαμάρα, βάσει του οποίου δηλώνουμε στο TESSy)
cat(paste("\"Εκτιμώμενος\" πληθυσμιακός παρονομαστής (σύνολο): ",prettyNum(subset(resAll, yearweek==tgtweek)$popest, big.mark=".", decimal.mark=","),"\n",sep=""))

cat("\nΤέλος!\n")

})   # Τέλος χρονομέτρησης

cat(paste("\nΧρόνος ανάλυσης:",round(timer["elapsed"],1),"δευτερόλεπτα\n\n"))

if(!interactive()) sink()
