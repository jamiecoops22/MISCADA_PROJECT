library("dplyr")
library("ggplot2")
library("gstat")
library("maps")

head(temp_data_full)

temp_1 <- subset(temp_data_full, Month %in% c(1,12,24,6,18))

temp_plot <- ggplot(temp_1) + # plot points
  geom_point(aes(x = Easting,y = Northing, # lon and lat
                 colour = y), # attribute color
             size = 5) + # make all points larger
  scale_colour_distiller("Diamond\nclarity") + # attach color scale
  xlab("Easting (km)") + # x-axis label
  ylab("Northing (km)") + # y-axis label
  #geom_path(data = map_data("state"), # add US states map
  #          aes(x = long, y = lat, group = group)) +
  facet_grid(~Month) + # facet by time
  #coord_fixed(xlim = c(-105, -75),
  #            ylim = c(25, 50)) + # zoom in
  theme_bw() # B&W theme
print(temp_plot)

