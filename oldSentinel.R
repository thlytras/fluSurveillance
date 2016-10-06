# Εφαρμογή στατιστικής επεξεργασίας των δεδομένων του συστήματος sentinel
# v3.0 © 2015, Θοδωρής Λύτρας
#
# Επεξεργασία των δεδομένων του παλιού sentinel
# Τελευταία αναθεώρηση: Οκτώβριος 2016

cat("\nΕπεξεργασία στοιχείων ΠΑΛΙΟΥ sentinel - παρακαλώ περιμένετε...\n\n")


# 1. Σύστημα sentinel ιδιωτών


cat("\nΑνάγνωση δηλώσεων ιδιωτών...\n\n")

tryCatch({   # Ανίχνευση λαθών
  sentnewa <- read.spss(paste(path_input,"oldSentinel/SENTNEWA.sav",sep=""), reencode="windows-1253", to.data.frame=TRUE, use.value.labels=FALSE)
}, error=function(err){
  stop(paste(readerrmsg[1],"oldSentinel/SENTNEWA.sav",readerrmsg[2],path_input))
})

names(sentnewa)<-tolower(names(sentnewa))

tryCatch({   # Ανίχνευση λαθών
  suppressWarnings( sent08 <- read.epiinfo(file(paste(path_input,"oldSentinel/sent08.rec",sep=""), encoding="windows-1253"), lower.case.names=TRUE) )
}, error=function(err){
  stop(paste(readerrmsg[1],"oldSentinel/sent08.rec",readerrmsg[2],path_input))
})

tryCatch({   # Ανίχνευση λαθών
  suppressWarnings( sent12 <- read.epiinfo(file(paste(path_input,"oldSentinel/sent12.rec",sep=""), encoding="windows-1253"), lower.case.names=TRUE) )
}, error=function(err){
  stop(paste(readerrmsg[1],"oldSentinel/sent12.rec",readerrmsg[2],path_input))
})


required_fields <- c("codeiat", "etos", "ebdo", "monvis", "tuevis", "wedvis", "thuvis", "frivis", "totvis", "gastot", "gritot", "hmekat")

if (FALSE %in% (required_fields %in% names(sentnewa))) stop("Το αρχείο SENTNEWA.sav δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")   # Ανίχνευση λαθών
if (FALSE %in% (required_fields %in% names(sent08))) stop("Το αρχείο sent08.rec δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")   # Ανίχνευση λαθών
if (FALSE %in% (required_fields %in% names(sent12))) stop("Το αρχείο sent12.rec δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")   # Ανίχνευση λαθών

# Διόρθωση για ημερομηνία καταχώρησης
sentnewa$hmekat <- as.Date(ISOdate(1582,10,14)  + sentnewa$hmekat)

sentinel<-rbind(sentnewa[,required_fields],sent08[,required_fields],sent12[,required_fields])
rm(sentnewa,sent08,sent12)

sentinel$astikot = 1
sentinel$yearweek = with(sentinel, etos*100 + ebdo)

sentinel<-subset(sentinel,codeiat!="" & !is.na(codeiat))

# Ανάγνωση της λίστας των ιδιωτών ιατρών
tryCatch({
  sentinel_doctors_10_2005 <- read.spss(paste(path_input,"oldSentinel/sentinel_doctors_10.2005.sav",sep=""), reencode="windows-1253", to.data.frame=TRUE, use.value.labels=FALSE)
}, error=function(err){
  stop(paste(readerrmsg[1],"oldSentinel/sentinel_doctors_10.2005.sav",readerrmsg[2],path_input))
})
names(sentinel_doctors_10_2005) <- tolower(names(sentinel_doctors_10_2005))
if (FALSE %in% (c("codeiat", "oldnew", "nomadil", "nuts", "eidi1") %in% names(sentinel_doctors_10_2005))) stop("Η λίστα των ιδιωτών ιατρών δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")

sentinel_doctors_10_2005 <- sentinel_doctors_10_2005[order(sentinel_doctors_10_2005$codeiat),]

# Ομαδοποιούμε τους γενικούς ιατρούς (eidi1==3), με τους παθολόγους (eidi1==1). Οι παιδίατροι είναι eidi1==2
sentinel_doctors_10_2005$eidi1[sentinel_doctors_10_2005$eidi1==3] = 1

# Συγχώνευση αρχείου δηλώσεων με τα στοιχεία των ιατρών
sentinel <- merge(sentinel, sentinel_doctors_10_2005[,c("codeiat", "oldnew", "nomadil", "nuts", "eidi1")], by="codeiat", all.x=TRUE)

sentinel<-subset(sentinel,!is.na(totvis))   # πετάμε έξω όσους δεν έχουν δηλώσει αριθμό συνολικών επισκέψεων ("αδρή" δήλωση αντί "αναλυτικής" δήλωσης)

sentinel$eid <- factor(sentinel$eidi1, levels=c(1:2))   # Μεταβλητή για την ειδικότητα

# Τριμάρουμε τα πεδία που δεν θα χρειαστούν άλλο...
sentinel <- transform(sentinel, oldnew=NULL, nuts=NULL, eidi1=NULL)
# ...και ορίζουμε από ποιό σύστημα είναι αυτές οι εγγραφές.
sentinel$systhma <- 1


# 2. Σύστημα sentinel Κέντρων Υγείας


cat("\nΑνάγνωση δηλώσεων από ΚΥ (θέλει λίγο χρόνο)...\n")

tryCatch({
  suppressWarnings( sentinel_KY <- read.epiinfo(file(paste(path_input,"oldSentinel/sentKY.rec",sep=""), encoding="windows-1253"), lower.case.names=TRUE) )
}, error=function(err){
  stop(paste(readerrmsg[1],"oldSentinel/sentKY.rec",readerrmsg[2],path_input))
})

required_fields <- c("nom", "eid", "etos", "ebdo", "monvis", "tuevis", "wedvis", "thuvis", "frivis", "totvis", "gastot", "gritot")
if (FALSE %in% (required_fields %in% names(sentinel_KY))) stop("Το αρχείο sentKY.rec δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")   # Ανίχνευση λαθών

sentinel_KY$yearweek = with(sentinel_KY, etos*100 + ebdo)
names(sentinel_KY)[names(sentinel_KY)=="nom"] = "nomadil"
sentinel_KY$astikot = 2
levels(sentinel_KY$nomadil)[levels(sentinel_KY$nomadil)=="Α4"] = "a4"
sentinel_KY$eid <- factor(sentinel_KY$eid, levels=c(6,5), labels=c(1,2))

# Τριμάρουμε τα πεδία που δεν θα χρειαστούν άλλο
sentinel_KY <- transform(sentinel_KY, aa=NULL, syst=NULL, nomos=NULL, eidikot=NULL, iatros=NULL, dilosi=NULL, arxebd=NULL, telebd=NULL, keycode=NULL, satvis=NULL, sunvis=NULL, laatot=NULL, erztot=NULL, koktot=NULL, anetot=NULL, ilatot=NULL, erytot=NULL, partot=NULL, totnos=NULL)
# ...και ορίζουμε από ποιό σύστημα είναι αυτές οι εγγραφές.
sentinel_KY$systhma <- 2


# 3. Σύστημα sentinel ΙΚΑ


cat("\nΑνάγνωση δηλώσεων από ΙΚΑ...\n\n")

tryCatch({   # Ανίχνευση λαθών
  ika5 <- read.spss(paste(path_input,"oldSentinel/IKA5.sav",sep=""), reencode="windows-1253", to.data.frame=TRUE, use.value.labels=FALSE)
}, error=function(err){
  stop(paste(readerrmsg[1],"oldSentinel/IKA5.sav",readerrmsg[2],path_input))
})

names(ika5)<-tolower(names(ika5))

tryCatch({   # Ανίχνευση λαθών
  suppressWarnings( ika5a <- read.epiinfo(file(paste(path_input,"oldSentinel/Ika5a.rec",sep=""), encoding="windows-1253"), lower.case.names=TRUE) )
}, error=function(err){
  stop(paste(readerrmsg[1],"oldSentinel/Ika5a.rec",readerrmsg[2],path_input))
})

ika5a <- transform(ika5a, codeiat=NULL, code1=NULL)
required_fields <- c("year1", "week1", "ika", "monvis", "tuevis", "wedvis", "thuvis", "frivis", "visits", "gripi", "eidik", "hmekat")

ika5$hmekat <- as.Date(NA) # Το IKA5.sav είναι γνωστό οτι δεν περιέχει πεδίο hmekat

if (FALSE %in% (required_fields %in% names(ika5))) stop("Το αρχείο IKA5.sav δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")   # Ανίχνευση λαθών
if (FALSE %in% (required_fields %in% names(ika5a))) stop("Το αρχείο Ika5a.rec δεν περιέχει τα σωστά πεδία! Αδυνατώ να συνεχίσω...")   # Ανίχνευση λαθών

ika<-rbind(ika5[,required_fields],ika5a[,required_fields])
rm(ika5,ika5a)
names(ika)<-recode(names(ika),c("year1","week1","ika","visits","gripi"),c("etos","ebdo","nomadil","totvis","gritot"))

ika$astikot = 1
ika$etos = as.integer(as.character(ika$etos))
ika$ebdo = as.integer(as.character(ika$ebdo))
ika$yearweek = with(ika, etos*100 + ebdo)
ika$codeiat = as.character(ika$nomadil)
levels(ika$nomadil)<-recode(levels(ika$nomadil),
    c("EG","AB","AD","AK","AL","AP","AT","BO","DR","GL","ES","ER","ZA","ZO","HR","IO","KA","KS","KE","KR","KO","MY","NI","PA","PI","PE","PY","RE","SE","TU","TR","SA","LE","VR","IN","HA","ΚΗ"),
    c("a4","a1","54","54","71","a1","a1","43","52","a1","82","a2","21","a1","91","33","17","56","a1","a4","a4","83","a4","13","a4","a4","54","93","62","54","44","84","24","a1","a1","a1","a1"))
# ΠΡΟΣΟΧΗ! Τα τελευταία αρχικά στη λίστα ("ΚΗ") είναι ΕΛΛΗΝΙΚΟΙ χαρακτήρες και όχι αγγλικοί.
ika<-ika[order(ika$nomadil),]

ika$eid <- factor(ika$eidik, levels=c("PA","PD"), labels=1:2)

ika$gastot <- NA

# Διώχνουμε τις μεταβλητές που δε θα χρειαστούμε
ika <- transform(ika, eidik=NULL)
# ...και ορίζουμε από ποιό σύστημα είναι αυτές οι εγγραφές.
ika$systhma <- 3



# Συγχώνευση στοιχείων και εξαγωγή ILI rate


cat("\nΣυγχώνευση στοιχείων από ΙΚΑ, ΚΥ και ιδιώτες (μπορεί να πάρει λίγο χρόνο)...\n")

# Τα αρχεία από τα τρία συστήματα έχουν πλέον τα ίδια πεδία. Τα συγχωνεύουμε...
big<-merge(sentinel,sentinel_KY,all=TRUE)
big<-merge(big,ika,all=TRUE)
# ...και αντιστοιχίζουμε τις πληροφορίες κάθε νομού.
big<-merge(big,nomos_populations,by="nomadil",all.x=TRUE)

# Σβήνουμε ότι δε χρειάζεται
big <- transform(big, nomokat=NULL, d_diam=NULL, didia_po=NULL, nomadil2=NULL)

big<-subset(big,!is.na(ast_p_nu))   # Έξω οι δηλώσεις που δε ξέρουμε το νομό τους
big<-subset(big,!is.na(totvis))   # Πετώ έξω όσα (από λάθος) έχουν missing επισκέψεις (δηλαδή δεν έχουν παρονομαστές). Για τους ιδιώτες αυτό έχει ήδη γίνει.
big$hmekat[big$hmekat<"2004-1-1" | big$hmekat>"2100-1-1"] <- NA  # Διόρθωση για εμφανώς λανθασμένες ημερομηνίες καταχώρησης

big <- subset(big, yearweek>200427 & yearweek<201440)
big <- subset(big, !is.na(totvis) & totvis>0)   # Discard if total visits zero or missing

# Υπολογίζουμε την πληρότητα της δήλωσης
plirotites_eidikotita_old <- with(big, table(
        systhma = factor(systhma, levels=1:3, labels=c("IDIO", "KEYG", "PEDY")), 
        eidikotita = factor(eid, levels=1:2, labels=c("PA","PD")), 
        yearweek))
plirotites_nuts_old <- with(big, table(
        systhma = factor(systhma, levels=1:3, labels=c("IDIO", "KEYG", "PEDY")), 
        nuts = factor(nuts, levels=1:4), 
        yearweek))


cat("\nΕξαγωγή ILI rate (βασικό μοντέλο)...\n")



resMainModelOld <- fitMainModel(big, NUTSpop, verbose=TRUE)
descrByWeekOld <- aggrByWeek(big)

resOld <- merge(resMainModelOld, descrByWeekOld)

cat("\nΕξαγωγή ILI rate (κατά NUTS 1)...\n")
resNutsOld <- fitGroupModel("nuts", big, NUTSpop, verbose=TRUE)
cat("\nΕξαγωγή ILI rate (κατά αστικότητα)...\n")
resAstyOld <- fitGroupModel("asty", big, NUTSpop, verbose=TRUE)

# Υπολογισμός "εκτιμώμενου" συνολικού πληθυσμού
# (Πρώην excelάκι Κατερέλου-Καλαμάρα, βάσει του οποίου δηλώνουμε στο TESSy)
abcland <- read.csv2(paste(path_input,"abcland.csv",sep=""),header=FALSE)
resOld$popest <- round(
        (resOld$pa/abcland[1,2]*sum(abcland[5:6,2])) + 
        (resOld$pd/abcland[2,2]*sum(abcland[3:4,2])))


# Σώζουμε τα αποτελέσματα
save(resOld, resNutsOld, resAstyOld,
        plirotites_eidikotita_old, plirotites_nuts_old, 
        file=paste(path_input, "oldSentinel.RData", sep=""))


cat("\nΤέλος! Συνεχίζουμε με την επεξεργασία των δεδομένων του NEOY sentinel...\n\n")
