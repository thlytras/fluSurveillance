#!/usr/local/bin/Rint

library(stats)
library(utils)
library(FluHMM)

load("output/latest_analysis.RData")

resAll$gri[1:match(201439, resAll$yearweek)] <- resAll$gri[1:match(201439, resAll$yearweek)]*4

myears <- resAll$yearweek[which((resAll$yearweek %% 100) == 42)] %/% 100

rates <- lapply(myears, function(my){
  r <- subset(resAll, (yearweek>=my*100+39) & (yearweek<=(my+1)*100+20))$gri
  names(r) <- subset(resAll, (yearweek>=my*100+39) & (yearweek<=(my+1)*100+20))$yearweek
  r
})
names(rates) <- myears


logsd.rates <- lapply(myears, function(my){
  r <- subset(resAll, (yearweek>=my*100+39) & (yearweek<=(my+1)*100+20))$log.gri.sd
  names(r) <- subset(resAll, (yearweek>=my*100+39) & (yearweek<=(my+1)*100+20))$yearweek
  r
})
names(logsd.rates) <- myears


y <- rev(names(rates))[1]

cat(sprintf("\nRunning FluHMM for season %s-%s", as.integer(y), as.integer(y)+1))

m <- FluHMM(rates[[y]], logSE=logsd.rates[[y]], K=10)

update(m, 100000, thin=10)

cairo_pdf("output/FluHMM-current.pdf", width=10, height=6)
plot(m)
dev.off()

cat("Finished!\n")

