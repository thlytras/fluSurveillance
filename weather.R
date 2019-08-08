library(FluMoDL)
load("data/map_perif.RData")  # Load region map of Greece (class SpatialPolygonsDataFrame)

if (TRUE) {  # Set to FALSE to NOT update weather
  if (!file.exists("data/gsod.RData")) {
    cat("Downloading all weather stations....\n")

    # Get all Greek weather stations from NOAA, with data after 2013-5-27 
    stations <- NOAA_countryStations("GR", from="2013-5-27")

    # Remove duplicate stations, plus Tymbaki which appears not to have NOAA data
    stations <- stations[-which(stations$station.name %in% 
      c("NAXOS", "KARPATHOS (AIRPORT)", "ZAKINTHOS", "EL VENIZELOS", "TYMBAKION (AIRPORT)")),]

    # Load package sp, and find out the region of each weather station 
    #   using its longitude and latitude
    library(sp)
    coordinates(stations) <- ~lon+lat
    stations$region <- map_perif$nuts2[over(stations, SpatialPolygons(map_perif@polygons))]

    # A few manual corrections in the region code
    stations$region[stations$station.name=="KITHIRA"] <- 25
    stations$region[stations$station.name=="NAFPLION"] <- 25
    stations$region[stations$station.name=="CYCLADES ISLANDS"] <- 42

    # Find out which is the last year
    last_year <- as.integer(format(Sys.Date(), "%Y"))

    # Get all temperatures, by year
    gsod <- list()
    for(y in 2013:last_year) {
      gsod[[as.character(y)]] <- NOAA_getGSOD(stations@data, y, c("station.name","region"))
    }
    gsod <- do.call(rbind, gsod)

    gsod <- subset(gsod, date>="2013-5-27")

    save(stations, gsod, file="data/gsod.RData")

    
  } else {
    
    cat("Updating latest weather... \n")
    load("data/gsod.RData")
    selYear <-  as.integer(format(Sys.Date(), "%Y")) - 1 + as.integer(as.integer(format(Sys.Date(), "%m"))>2):1

    newgsod <- NOAA_getGSOD(stations@data, selYear, c("station.name","region"))
    gsod <- rbind(gsod, newgsod)
    
    # Removing duplicates
    gsod <- gsod[!duplicated(gsod),]
    gsod <- gsod[!rev(duplicated(gsod[nrow(gsod):1,c("date","station.name")])),]
    
    # Sorting
    gsod <- gsod[with(gsod, order(date, region, station.name)),]
    rownames(gsod) <- NULL
    
    save(stations, gsod, file="data/gsod.RData")

  }
  
} else { 
  cat("NOT updating weather, downloading already saved weather tations...\n")
  load("data/gsod.RData") 
}



cat("Aggregating temperatures...\n")
system.time({
  mdTemp <- aggregate(gsod[, "temp", drop=FALSE], gsod[,c("region","date")], mean, na.rm=TRUE)
  mdTemp$pop <- map_perif$pop[match(mdTemp$region, map_perif$nuts2)]
  mdTemp$sum_pop <- with(mdTemp, tapply(pop, date, sum, na.rm=TRUE)[as.character(date)])
  mdTemp$w_temp <- with(mdTemp, temp*pop/sum_pop)
  mdTemp <- aggregate(mdTemp[, "w_temp", drop=FALSE], mdTemp[,"date",drop=FALSE], sum, na.rm=TRUE)
  names(mdTemp)[2] <- "temp"
})

cat("Filling in missing values...\n")
load("data/oldtemp.RData")
names(oldtemp)[1:2] <- names(mdTemp)
mdTemp <- rbind(mdTemp, oldtemp[which(!(oldtemp$date %in% mdTemp$date)),1:2])
mdTemp <- mdTemp[order(mdTemp$date),]

mdTemp <- data.frame(
  date = seq(min(mdTemp$date), max(mdTemp$date), by="days"),
  temp = mdTemp$temp[match(seq(min(mdTemp$date), max(mdTemp$date), by="days"), mdTemp$date)])


cat("Saving mean daily temperatures... ")
save(mdTemp, file="data/mdTemp.RData")

cat("Done!\n\n")
