library(rCharts)

path_output = "./output/"

if (!exists("sentinel_graph")) {
  load(paste(path_output, "latest_analysis.RData", sep=""))
}

web_graph <- function(years, col=rainbow(length(years)), 
	yaxis2=NA, mult=1, ygrid=0, lty=rep(1,length(years)), lwd=rep(1,length(years)),
	ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις",
	ylab2=NA, ylab2rot=TRUE)
{
  ylab <- gsub("\n", "<br/>", ylab); ylab2 <- gsub("\n", "<br/>", ylab2)
  maxwk <- ifelse(sum(as.integer(isoweek(as.Date(paste(years,"-12-31",sep="")))==53))>0, 53, 52)
  set <- sapply(years, function(x){ratechart[as.character(c((x*100+40):(x*100+maxwk),((x+1)*100+1):((x+1)*100+20)))]})
  set <- as.data.frame(set)
  colnames(set) <- paste(years, years+1, sep="-")
  set$week <- factor(as.integer(rownames(set))%%100, levels=c(40:maxwk,1:20))
  limrate <- (max(set[,-ncol(set)],na.rm=TRUE)%/%10+2)*10
  if(!is.na(yaxis2[1])) {
    maxes <- apply(set[,-ncol(set)], 2, max, na.rm=TRUE)
    maxes[match(yaxis2,years)] <- maxes[match(yaxis2,years)]/mult
    limrate <- (max(maxes, na.rm=TRUE)%/%10+2)*10
    limrate2 <- limrate*mult
  }
  set <- reshape(set, direction="long", varying=list(1:3), v.names="rate")
  set$year <- paste(years, years+1, sep="-")[set$time]
  set$time <- NULL; set$id <- NULL; rownames(set) <- NULL
  suppressWarnings({ h <- hPlot(rate ~ week, group="year", data=set) })
  for(i in 1:length(years)) {
    h$params$series[[i]]$color <- col[i]
    h$params$series[[i]]$lineWidth <- lwd[i]
    h$params$series[[i]]$dashStyle <- lty[i]
    h$params$series[[i]]$marker$symbol <- "circle"
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
  h$tooltip(formatter="#! function(){ return('<b>' + this.series.name + '</b><br/>Εβδομάδα ' + this.x + ': <b>' + this.y + '</b>') } !#")
  return(h)
}

h <- web_graph(ytp, col=apply(col2rgb(c("black", "navyblue","red3"))/255, 2, function(x)rgb(x[1],x[2],x[3])), 
      lty=c("shortdot", "shortdot", "solid"), lwd=c(1.5,1.5,3),
      yaxis2=ytp[ytp<2014], mult=1/5, 
      ylab="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις\n(Νέο σύστημα επιτήρησης, 2014-2015)",
      ylab2="Κρούσματα γριπώδους συνδρομής ανά 1000 επισκέψεις\n(Παλιό σύστημα επιτήρησης, 2012-2013, 2013-2014)")

h$save(paste(path_output, "rChart.html", sep=""))

