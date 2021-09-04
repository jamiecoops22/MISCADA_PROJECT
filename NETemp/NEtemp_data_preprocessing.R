library(spBayes)
library(dplyr)

data("NETemp.dat")
View(NETemp.dat)

time <- rep(1, 356)
coords_1 <- data.frame(NETemp.dat$UTMX, NETemp.dat$UTMY, time, NETemp.dat$y.1)
View(coords_1)
temp_mat <- as.matrix(NETemp.dat)
View(temp_mat)
dmat <- c()
for (i in 1:356){
  loc_mat <- cbind(rep(temp_mat[i,2], 129), rep(temp_mat[i,3], 129), c(1:129), rep(1, 129), temp_mat[i,4:132])
  dmat <- rbind(dmat, loc_mat)
}
View(dmat)

set.seed(11)
reduced <- sample(1:356, 100, replace = F)
ord <- order(reduced)
reduced <- reduced[ord]

reduced_dmat <- c()
for (i in reduced){
  loc1_mat <- cbind(rep(temp_mat[i,2], 24)/1000, rep(temp_mat[i,3], 24)/1000, c(1:24), rep(1,24), rep(temp_mat[i,1], 24), temp_mat[i,4:27])
  reduced_dmat <- rbind(reduced_dmat, loc1_mat)
}

View(reduced_dmat)

temp_data_full <- data.frame(coords = reduced_dmat[,1:3],
                             x = reduced_dmat[,4:5],
                             y = reduced_dmat[,6])
View(temp_data_full)

temp_data_full <- temp_data_full %>%
  rename(
    Month = coords.3,
    Easting = coords.1,
    Northing = coords.2
  )

# choose model and holdout sets

set.seed(1)
ho <- sort(sample(1:nrow(reduced_dmat), 240))

temp_data_holdout <- data.frame(coords = reduced_dmat[ho,1:3],
                                x = reduced_dmat[ho,4:5],
                                y = reduced_dmat[ho,6])

temp_data_holdout <- temp_data_holdout %>%
  rename(
    Month = coords.3,
    Easting = coords.1,
    Northing = coords.2
  )
View(temp_data_holdout)

temp_data_model <- data.frame(coords = reduced_dmat[-ho,1:3],
                                x = reduced_dmat[-ho,4:5],
                                y = reduced_dmat[-ho,6])

temp_data_model <- temp_data_model %>%
  rename(
    Month = coords.3,
    Easting = coords.1,
    Northing = coords.2
  )
View(temp_data_model)


write.table(reduced_dmat[ho,1:3], "coords.ho", sep="\t", row.names=F, col.names=F)
write.table(reduced_dmat[ho,6], "y.ho", sep="\t", row.names=F, col.names=F)
write.table(reduced_dmat[ho,4:5], "x.ho", sep="\t", row.names=F, col.names=F)

write.table(reduced_dmat[-ho,1:3], "coords.mod", sep="\t", row.names=F, col.names=F)
write.table(reduced_dmat[-ho,6], "y.mod", sep="\t", row.names=F, col.names=F)
write.table(reduced_dmat[-ho,4:5], "x.mod", sep="\t", row.names=F, col.names=F)

temp_coords <- temp_data_full %>%
    select(Easting, Northing)

temp_coords <- temp_coords %>%
  mutate(Easting = Easting*1000,
         Northing = Northing*1000)

write.csv(temp_coords,"temp_coords.csv", row.names = TRUE)

