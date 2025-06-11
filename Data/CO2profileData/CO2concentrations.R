
setwd("~/SOIL-R/FieldWork/2024/Prades/CO2profileData/")

allfiles<-list.files(path=".", pattern=".dat")
pradesFiles<-allfiles[45:115] # These are only the files produced after May 30

allAvgData<-lapply(pradesFiles, FUN=read.csv, skip=4, na.strings = "NAN", header=FALSE)

mergedData<-Reduce(function(...) merge(..., all=TRUE), allAvgData)
cols<-c("TIMESTAMP","RECORD","MeasCycleIDX_Avg","CO2Lvl1_Avg","H2OLvl1_Avg","FlowLvl1_Avg",
        "CO2Lvl1_Std","H2OLvl1_Std","FlowLvl1_Std","CO2Lvl2_Avg","H2OLvl2_Avg","FlowLvl2_Avg",
        "CO2Lvl2_Std","H2OLvl2_Std","FlowLvl2_Std","CO2Lvl3_Avg","H2OLvl3_Avg","FlowLvl3_Avg",
        "CO2Lvl3_Std","H2OLvl3_Std","FlowLvl3_Std","CO2Lvl4_Avg","H2OLvl4_Avg","FlowLvl4_Avg",
        "CO2Lvl4_Std","H2OLvl4_Std","FlowLvl4_Std","CO2Lvl5_Avg","H2OLvl5_Avg","FlowLvl5_Avg","CO2Lvl5_Std",
        "H2OLvl5_Std","FlowLvl5_Std","CO2Lvl6_Avg","H2OLvl6_Avg","FlowLvl6_Avg","CO2Lvl6_Std",
        "H2OLvl6_Std","FlowLvl6_Std","CO2Lvl7_Avg","H2OLvl7_Avg","FlowLvl7_Avg","CO2Lvl7_Std",
        "H2OLvl7_Std","FlowLvl7_Std","BattV_Min","PTemp_Avg")
names(mergedData)<-cols

tm_UTC<-as.POSIXct(mergedData$TIMESTAMP, format="%Y-%m-%d %H:%M:%S", tz="UTC")
tm_CEST<-as.POSIXct(tm_UTC, tz="CET")

CO2cols<-seq(4, 43, by=6)
names(mergedData)[CO2cols] # Confirm these are the correct CO2 columns
heights<-paste("Level", 1:7)

pal<-rainbow(n=length(CO2cols))

pdf("CO2concentrations_withPeak.pdf", width = 7*sqrt(2))
par(mar=c(4,4,1,1))
plot(tm_CEST, mergedData[,4], type="l", ylim=c(420, 800), xlim=c(as.POSIXct("2024-04-30"), as.POSIXct("2024-05-04")), 
     xlab="", ylab="CO2 concentration (ppm)", bty="n")
for(i in 1:length(CO2cols)){
  lines(tm_CEST, mergedData[,CO2cols[i]], col=pal[i])
}
abline(v=as.POSIXct(c("2024-04-30 12:00:00", "2024-05-01 12:00:00", "2024-05-02 12:00:00", "2024-05-03 12:00:00")), lty=2)
abline(v=as.POSIXct(c("2024-05-01 06:00:00", "2024-05-02 06:00:00", "2024-05-03 06:00:00")), lty=3)
legend("topright", legend=c(heights, "Noon", "6:00 am"), lty=c(rep(1, 7), 2, 3), col=c(pal, 1, 1), bty="n")
dev.off()

pdf("CO2concentrations.pdf", width = 7*sqrt(2))
par(mar=c(4,4,1,1))
plot(tm_CEST, mergedData[,4], type="l", ylim=c(420, 460), xlim=c(as.POSIXct("2024-04-30"), as.POSIXct("2024-05-04")), 
     xlab="", ylab="CO2 concentration (ppm)", bty="n")
for(i in 1:length(CO2cols)){
  lines(tm_CEST, mergedData[,CO2cols[i]], col=pal[i])
}
abline(v=as.POSIXct(c("2024-04-30 12:00:00", "2024-05-01 12:00:00", "2024-05-02 12:00:00", "2024-05-03 12:00:00")), lty=2)
abline(v=as.POSIXct(c("2024-05-01 06:00:00", "2024-05-02 06:00:00", "2024-05-03 06:00:00")), lty=3)
legend("topright", legend=c(heights, "Noon", "6:00 am"), lty=c(rep(1, 7), 2, 3), col=c(pal, 1, 1), bty="n")
dev.off()
