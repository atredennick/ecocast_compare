library(ncdf4)
library(ncdf4.helpers)
library(PCICt)
library(tidyverse)

climate_filepath <- paste0("~/Desktop/snd_LImon_MIROC-ESM_historical_r1i1p1_185001-200512.nc")
climate_output <- nc_open(climate_filepath)

lon <- ncvar_get(climate_output, varid = "lon")
lat <- ncvar_get(climate_output, varid = "lat")
climate_output$dim$time$units
climate_output$dim$time$calendar

snd_time <- nc.get.time.series(climate_output, v = "snd",
                               time.dim.name = "time")
snd_time[c(1:3, length(snd_time) - 2:0)]
snd <- ncvar_get(climate_output, "snd")

time_index <- which(format(snd_time, "%Y-%m-%d") == "2005-01-16")
snd <- nc.get.var.subset.by.axes(climate_output, "snd",
                                 axis.indices = list(T = time_index))

library(ggmap)
library(viridis)

expand.grid(lon, lat) %>%
  rename(lon = Var1, lat = Var2) %>%
  mutate(lon = ifelse(lon > 180, -(360 - lon), lon),
         snd = as.vector(snd)) %>% 
  ggplot() + 
  geom_raster(aes(x = lon, y = lat, fill = snd)) + 
  borders("world", colour="black", fill=NA) + 
  scale_fill_viridis(name = "Snow Depth", direction = -1) + 
  theme_void() + 
  coord_quickmap() + 
  ggtitle("Modeled snow depth on a day in January 2005") 