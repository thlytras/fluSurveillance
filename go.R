require(plotrix)

aggr2$gri_w_var <- ifelse(aggr2$astikot==1,
  ((aggr2$gritot/aggr2$totvis)*(1-aggr2$gritot/aggr2$totvis)/aggr2$totvis)*(1000*aggr1$ast_p_nu)^2,
  ((aggr2$gritot/aggr2$totvis)*(1-aggr2$gritot/aggr2$totvis)/aggr2$totvis)*(1000*aggr1$agr_p_nu)^2)
aggr2$gri_w[is.na(aggr2$gri_w)] = 0

astAggr3 <- aggregate(aggr2[,c("gri_w","gritot","totvis","gri_w_var")],by=list(yearweek=aggr2$yearweek, astikot=aggr2$astikot), sum, na.rm=TRUE)
astAggr3 <- subset(astAggr3, yearweek>200427 & yearweek<210000)

#with(subset(astAggr3, astikot==1), plot(gri_w, type="l", ylim=c(0,100), col="red"))
#with(subset(astAggr3, astikot==2), points(gri_w, type="l", col="blue"))

spAggr1 <- aggregate(sentinelBig[,c("etos", "ebdo", "ast_p_nu", "agr_p_nu", "as_p_nu1", "ag_p_nu1", "as_p_nu2", "ag_p_nu2")], by=list(monada=sentinelBig$monada, astikot=sentinelBig$astikot, nuts=sentinelBig$nuts, yearweek=sentinelBig$yearweek), mean, na.rm=TRUE)

spAggr2 <- aggregate(sentinelBig[,c("totvis", "gritot")], by=list(monada=sentinelBig$monada, astikot=sentinelBig$astikot, nuts=sentinelBig$nuts, yearweek=sentinelBig$yearweek), sum, na.rm=TRUE)

spAggr2$gritot[is.na(spAggr2$gritot)] = 0
spAggr2$gri <- (spAggr2$gritot/spAggr2$totvis)*1000
spAggr2$gri[is.na(spAggr2$gri)] = 0

spAggr2$gri_w <- ifelse(spAggr2$astikot==1,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$ast_p_nu,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$agr_p_nu)
spAggr2$gri_w[is.na(spAggr2$gri_w)] = 0

spAggr2$gri_w_var <- ifelse(spAggr2$astikot==1,
  ((spAggr2$gritot/spAggr2$totvis)*(1-spAggr2$gritot/spAggr2$totvis)/spAggr2$totvis)*(1000*spAggr1$ast_p_nu)^2,
  ((spAggr2$gritot/spAggr2$totvis)*(1-spAggr2$gritot/spAggr2$totvis)/spAggr2$totvis)*(1000*spAggr1$agr_p_nu)^2)
spAggr2$gri_w_var[is.na(spAggr2$gri_w)] = 0


spAggr2$gri_w1 <- ifelse(spAggr2$astikot==1,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$as_p_nu1,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$ag_p_nu1)
spAggr2$gri_w1[is.na(spAggr2$gri_w1)] = 0

spAggr2$gri_w2 <- ifelse(spAggr2$astikot==1,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$as_p_nu2,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$ag_p_nu2)
spAggr2$gri_w2[is.na(spAggr2$gri_w2)] = 0

spAggr3 <- aggregate(spAggr2[,c("gri_w","gritot","totvis","gri_w_var")],by=list(yearweek=spAggr2$yearweek, monada=spAggr2$monada), sum, na.rm=TRUE)

spAggr3 <- subset(spAggr3, yearweek>201402)
astAggr3 <- subset(astAggr3, yearweek>201402)

wklen <- nrow(subset(spAggr3, monada=="KEYG"))
jitf <- 0.03

#dev.new()
par(mar=c(6,4,4,2))
with(subset(spAggr3, monada=="KEYG"), plot(1:wklen, gri_w,  type="l", ylim=c(0,130), col="blue", bty="l", ylab="ILI rate", xlab=NA, xaxt="n", main="ILI rate κατά κατηγορία μονάδας"))
with(subset(spAggr3, monada=="KEYG"), plotCI(1:wklen, gri_w, uiw=1.96*sqrt(gri_w_var), col="blue", sfrac=0.005, pch=19, cex=0.5, add=TRUE))
with(subset(spAggr3, monada=="PEDY"), points(1:wklen+1*jitf, gri_w, type="l", col="darkred"))
with(subset(spAggr3, monada=="PEDY"), plotCI(1:wklen+1*jitf, gri_w, uiw=1.96*sqrt(gri_w_var), col="darkred", sfrac=0.005, pch=19, cex=0.5, add=TRUE))
with(subset(spAggr3, monada=="IDIO"), points(1:wklen-1*jitf, gri_w, type="l", col="red"))
with(subset(spAggr3, monada=="IDIO"), plotCI(1:wklen-1*jitf, gri_w, uiw=1.96*sqrt(gri_w_var), col="red", sfrac=0.005, pch=19, cex=0.5, add=TRUE))
with(subset(astAggr3, astikot==1), points(1:wklen, gri_w, type="l", col="deeppink", lty="dotted"))
with(subset(astAggr3, astikot==1), plotCI(1:wklen, gri_w, uiw=1.96*sqrt(gri_w_var), col="deeppink", lty="dotted", sfrac=0.005, pch=19, cex=0.5, add=TRUE))
legend("topleft", c("ΚΥ", "ΠΕΔΥ", "Ιδιώτες", "Αστικός πληθυσμός"), col=c("blue","darkred","red","deeppink"), lty=c(rep("solid",3), "dotted"), pt.cex=0.5, pch=19, bty="n", inset=0.03)
axis(1, at=1:wklen, label=subset(spAggr3, monada=="KEYG")$yearweek, las=2, tick=FALSE, hadj=0.7)
mtext("Εβδομάδα", side=1, line=4)