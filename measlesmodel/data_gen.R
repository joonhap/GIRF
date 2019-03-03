data <- read.csv('measlesUKUS.csv')
UKdata <- subset(x=data, data$country=="UK")

UKcities <- levels(factor(UKdata$loc))
years <- 1944:1964

## BIRTH DATA ##
citybirths <- function(cityname) {
    sapply(years, function(yr) round(sum(subset(x=UKdata, UKdata$loc==cityname & UKdata$year==yr)$rec)) )
}
UKbirths <- cbind(years, sapply(UKcities, citybirths))
colnames(UKbirths) <- c("time", UKcities)
write.csv(UKbirths, file="UKbirths.csv", row.names=FALSE)

## POPULATION DATA ##
citypop <- function(cityname) {
    sapply(years, function(yr) round(subset(x=UKdata, UKdata$loc==cityname & UKdata$year==yr & UKdata$biweek==1)$pop) )
}
UKpop <- cbind(years, sapply(UKcities, citypop))
colnames(UKpop) <- c("time", UKcities)
write.csv(UKpop, file="UKpop.csv", row.names=FALSE)

## MEASLES CASE DATA ##
decimalyear <- subset(x=UKdata, UKdata$loc==UKcities[1] & UKdata$year %in% years)$decimalYear
citycase <- function(cityname) {
    subset(x=UKdata, UKdata$loc==cityname & UKdata$year %in% years)$cases
}
UKcases <- cbind(decimalyear, sapply(UKcities, citycase))
colnames(UKcases) <- c("decimalYear", UKcities)
write.csv(UKcases, file="UKmeasles.csv", row.names=FALSE)

## DISTANCE DATA ##
UKlonlat <- sapply(UKcities, function(cityname) unlist(subset(x=UKdata, UKdata$loc==cityname & UKdata$year==years[1] & UKdata$biweek==1)[,c("lon", "lat")]) )

library(geosphere)
distMat <- matrix(0, nrow=length(UKcities), ncol=length(UKcities))
for (i in 1:length(UKcities))
    for (j in 1:length(UKcities))
        distMat[i,j] <- distHaversine(p1=UKlonlat[,i], p2=UKlonlat[,j])
write.table(x=distMat, file="DistMatrix.txt", row.names=FALSE, col.names=FALSE, sep=" ")
