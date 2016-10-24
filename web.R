library(rCharts)

#path_output = "./output/"

#if (!exists("sentinel_graph")) {
#  load(paste(path_output, "latest_analysis.RData", sep=""))
#}

web_graph <- function(years, col=rainbow(length(years)), 
	yaxis2=NA, mult=1, ygrid=0, lty=rep(1,length(years)), lwd=rep(1,length(years)),
	ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις",
	ylab2=NA, ylab2rot=TRUE, ci=TRUE, alpha=0.15)
{
  ylab <- gsub("\n", "<br/>", ylab); ylab2 <- gsub("\n", "<br/>", ylab2)
  ratechart <- resAll$gri; names(ratechart) <- resAll$yearweek
  ratechart_logsd <- resAll$log.gri.sd; names(ratechart_logsd) <- resAll$yearweek
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

  set <- cbind(set, as.integer(rownames(set))%%100)
  set[,ncol(set)][is.na(set[,ncol(set)])] <- 53
  rownames(set)<-NULL
  set <- as.data.frame(set)
  colnames(set) <- c(paste(years, years+1, sep="-"), "week")
  set$week <- factor(set$week, levels=c(40:maxwk,1:20))
  set <- reshape(set, direction="long", varying=list(1:(ncol(set)-1)), v.names="rate")
  set$year <- paste(years, years+1, sep="-")[set$time]
  set$time <- NULL; set$id <- NULL; rownames(set) <- NULL

#  suppressWarnings({ h <- hPlot(rate ~ week, group="year", data=set) })
#  for(i in 1:length(years)) {
#    h$params$series[[i]]$color <- col[i]
#    h$params$series[[i]]$lineWidth <- lwd[i]
#    h$params$series[[i]]$dashStyle <- lty[i]
#    h$params$series[[i]]$marker$symbol <- "circle"
#  }
  
  h <- Highcharts$new()
  for(i in 1:length(years)) {
      h$series(data=toJSONArray2(subset(set, year==levels(factor(set$year))[i])[,c("week","rate")], 
                names=F, json=F), name=levels(factor(set$year))[i],
          type="line", color=col[i], lineWidth=lwd[i],
          dashStyle=lty[i], zIndex=i, marker=list(symbol="circle"))
  }
  
  if (!is.na(yaxis2[1])) {
    h$yAxis(list(list(title=list(text=ylab, margin=30), min=0, max=limrate), 
      list(title=list(text=ylab2, margin=30), min=0, max=limrate2, opposite="true")))
    for(i in 1:length(years)) {
      h$params$series[[i]]$yAxis <- as.integer(yaxis2[i] %in% years)
    }
  } else {
    h$yAxis(title=list(text=ylab, margin=30), min=0, max=limrate)
  }
  h$params$xAxis[[1]]$title$text <- "Εβδομάδα"
  h$tooltip(formatter="#! function(){ return('<b>' + this.series.name + '</b><br/>Εβδομάδα ' + this.x + ': <b>' + Math.round(this.y*100)/100 + '</b>') } !#")

  if (sum(ci)>0) {
    set_logsd <- cbind(set_logsd, as.integer(rownames(set_logsd))%%100)
    set_logsd[,ncol(set_logsd)][is.na(set_logsd[,ncol(set_logsd)])] <- 53
    rownames(set_logsd)<-NULL
    set_logsd <- as.data.frame(set_logsd)
    colnames(set_logsd) <- c(paste(years, years+1, sep="-"), "week")
    set_logsd <- reshape(set_logsd, direction="long", varying=list(1:(ncol(set_logsd)-1)), v.names="rate")
    set$min <- exp(log(set$rate) - 1.96*set_logsd$rate)
    set$max <- exp(log(set$rate) + 1.96*set_logsd$rate)
    set$col <- col[as.integer(factor(set$year))]
    set$alpha <- alpha[as.integer(factor(set$year))]
    by(set, set$year, function(x){
      h$series(data=toJSONArray2(x[,c("week","min","max")], names=F, json=F), 
          type="arearange", fillOpacity=unique(x$alpha), color=unique(x$col),
          lineWidth=0, zIndex=-1, showInLegend=FALSE, enableMouseTracking=FALSE)
    })
    for(i in length(years)+(1:length(years))) {
      h$params$series[[i]]$yAxis <- as.integer(yaxis2[i-length(years)] %in% years)
    }
    
  }
  
  return(h)
}

#h <- web_graph(c(2013,2014,2015), 
#        col=apply(col2rgb(c("navyblue", "red", "green"))/255, 2, function(x)rgb(x[1],x[2],x[3])), 
#        lty=c("shortdot", "solid"), lwd=c(2), yaxis2=2013, mult=1/5,
#        ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις\n(Νέο σύστημα επιτήρησης)")

#h$save(paste(path_output, "rChart.html", sep=""))

