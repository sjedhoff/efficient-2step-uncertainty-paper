################################################################################
##                                                                            ##
##                    NYC Case Study: Prepare the data                        ##
##                                                                            ##
################################################################################

## load the dataset
data <- readRDS("case_study_multiple_imputation/data/nyc_2016_03_07.rds")

## relevant variables
vars <- c("trip_duration", "osm_distance", "osm_duration", "pickup_latitude",
          "pickup_longitude", "dropoff_latitude", "dropoff_longitude")  

## split by hour
hour <- lapply(0:23, function(i) data[data$pickup_hour == i,])

## Impute the data
source("simulation_multiple_imputation/utils.R")
m <- 20
set.seed(31415)
for(i in 0:23){
  data <- hour[[i+1]]
  imp <- impute_data(m = m, data = data)
  # transform all imputed datasets
  imp_scaled <- list()
  for(j in 1:m){
    dat <- imp[[j]][,vars]
    imp_scaled[[j]] <- data.frame("trip_duration_log" = log(dat$trip_duration),
                                  "osm_distance_scale" = scale(dat$osm_distance),
                                  "osm_duration_scale" = scale(dat$osm_duration),
                                  "pickup_latitude_scale" = scale(dat$pickup_latitude),
                                  "pickup_longitude_scale" = scale(dat$pickup_longitude),
                                  "dropoff_latitude_scale" = scale(dat$dropoff_latitude),
                                  "dropoff_longitude_scale" = scale(dat$dropoff_longitude))
  }
  saveRDS(imp_scaled, file = paste0("case_study_multiple_imputation/data/imp_2016_03_07_",i,".RDS"))
}
