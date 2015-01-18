
astAggr3 <- aggregate(aggr2[,c("gri_w","gritot","totvis")],by=list(yearweek=aggr2$yearweek, astikot=aggr2$astikot), sum, na.rm=TRUE)
astAggr3 <- subset(astAggr3, yearweek>200427 & yearweek<210000)

with(subset(astAggr3, astikot==1), plot(gri_w, type="l", ylim=c(0,100), col="red"))
with(subset(astAggr3, astikot==2), points(gri_w, type="l", col="blue"))

spAggr1 <- aggregate(sentinelBig[,c("etos", "ebdo", "ast_p_nu", "agr_p_nu", "as_p_nu1", "ag_p_nu1", "as_p_nu2", "ag_p_nu2")], by=list(monada=sentinelBig$monada, astikot=sentinelBig$astikot, nuts=sentinelBig$nuts, yearweek=sentinelBig$yearweek), mean, na.rm=TRUE)

spAggr2 <- aggregate(sentinelBig[,c("totvis", "gritot")], by=list(monada=sentinelBig$monada, astikot=sentinelBig$astikot, nuts=sentinelBig$nuts, yearweek=sentinelBig$yearweek), sum, na.rm=TRUE)

spAggr2$gritot[is.na(spAggr2$gritot)] = 0
spAggr2$gri <- (spAggr2$gritot/spAggr2$totvis)*1000
spAggr2$gri[is.na(spAggr2$gri)] = 0

spAggr2$gri_w <- ifelse(spAggr2$astikot==1,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$ast_p_nu,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$agr_p_nu)
spAggr2$gri_w[is.na(spAggr2$gri_w)] = 0

spAggr2$gri_w1 <- ifelse(spAggr2$astikot==1,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$as_p_nu1,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$ag_p_nu1)
spAggr2$gri_w1[is.na(spAggr2$gri_w1)] = 0

spAggr2$gri_w2 <- ifelse(spAggr2$astikot==1,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$as_p_nu2,
  (spAggr2$gritot/spAggr2$totvis)*1000*spAggr1$ag_p_nu2)
spAggr2$gri_w2[is.na(spAggr2$gri_w2)] = 0

spAggr3 <- aggregate(spAggr2[,c("gri_w","gritot","totvis")],by=list(yearweek=spAggr2$yearweek, monada=spAggr2$monada), sum, na.rm=TRUE)

spAggr3 <- subset(spAggr3, yearweek>201402)
astAggr3 <- subset(astAggr3, yearweek>201402)

dev.new()
with(subset(spAggr3, monada=="KEYG"), plot(gri_w, type="l", ylim=c(0,100), col="blue"))
with(subset(spAggr3, monada=="PEDY"), points(gri_w, type="l", col="darkred"))
with(subset(spAggr3, monada=="IDIO"), points(gri_w, type="l", col="red"))
with(subset(astAggr3, astikot==1), points(gri_w, type="l", col="deeppink", lty="dotted"))

