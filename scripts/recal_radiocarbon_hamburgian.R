
# Loading the Rcarbon package

library(rcarbon)

# Re-calibrating the weighted averages, of radiocarbon dates from archaeological contexts related
# to the Hamburgian phenomenon, calculated by Pedersen, Maier and Riede (2019), using the IntCal20 calibration curve

x <- calibrate(x = c(12363, 12229, 11719), errors = c(22, 18, 40),
               calCurves = "intcal20",
               ids = c("Pulse 1", "Pulse 2", "Pulse 3"))

# Retrieving the median and one- and two Sigma ranges for the calibrated dates
c14_sum <- summary(x)



# Plotting all three re-calibrated weighted averages and highlightning the 95% higher posterior density interval

par(mfrow=c(2,2))
multiplot(x)
plot(x,1, HPD=TRUE,credMass=0.95, main="Pulse 1")
plot(x,2, HPD=TRUE,credMass=0.95, main="Pulse 2")
plot(x,3, HPD=TRUE,credMass=0.95, main="Pulse 3")


# Exporting summary table
write.csv(c14_sum, file="tables/supplementary/c14_summary.csv",
          row.names = FALSE)
