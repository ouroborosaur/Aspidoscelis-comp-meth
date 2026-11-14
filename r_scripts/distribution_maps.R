# install packages
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# read KMLs
ang <- st_read("angusticeps.kml")
dep <- st_read("deppii.kml")
coz <- st_read("cozumelus.kml")

# inspect KML
# print(ang)
# plot(st_geometry(ang))

# colors
GROUP_COLORS <- c(
  "A. angusticeps" = "#CC79A7",
  "A. deppii" = "#56B4E9",
  "A. cozumelus" = "#009E73"
)


# angusticeps only
ggplot() +
  
  # geographic background
  geom_sf(
    data = world,
    fill = "grey95",
    color = "grey60",
    linewidth = 0.2
  ) +
  
  # species distributions
  geom_sf(
    data = ang,
    aes(fill = "A. angusticeps"),
    #fill = "#CC79A7",
    alpha = 0.5,
    color = NA
  ) +
  
  coord_sf(
    xlim = c(-110, -80),
    ylim = c(5, 30),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Species"
  ) +
  
  scale_fill_manual(
    values = GROUP_COLORS
  )

# deppii only
ggplot() +
  
  # geographic background
  geom_sf(
    data = world,
    fill = "grey95",
    color = "grey60",
    linewidth = 0.2
  ) +
  
  geom_sf(
    data = dep,
    aes(fill = "A. deppii"),
    #fill = "#56B4E9",
    alpha = 0.5,
    color = NA
  ) +
  
  coord_sf(
    xlim = c(-110, -80),
    ylim = c(5, 30),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Species"
  ) +
  
  scale_fill_manual(
    values = GROUP_COLORS
  )

theme_classic()


# plot together (zoomed out)
ggplot() +
  
  # geographic background
  geom_sf(
    data = world,
    fill = "grey95",
    color = "grey60",
    linewidth = 0.2
  ) +
  
  # species distributions
  geom_sf(
    data = ang,
    aes(fill = "A. angusticeps"),
    #fill = "#CC79A7",
    alpha = 0.5,
    color = NA
  ) +
  
  geom_sf(
    data = dep,
    aes(fill = "A. deppii"),
    #fill = "#56B4E9",
    alpha = 0.5,
    color = NA
  ) +
  
  geom_sf(
    data = coz,
    aes(fill = "A. cozumelus"),
    #fill= "#009E73",
    alpha = 0.7,
    color = NA,
    linewidth = 0.3
  ) +
  
  coord_sf(
    xlim = c(-110, -80),
    ylim = c(5, 30),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Species"
  ) +
  
  scale_fill_manual(
    values = GROUP_COLORS
  )

# plot together (zoomed out)
ggplot() +
  
  # geographic background
  geom_sf(
    data = world,
    fill = "grey95",
    color = "grey60",
    linewidth = 0.2
  ) +
  
  # species distributions
  geom_sf(
    data = ang,
    aes(fill = "A. angusticeps"),
    #fill = "#CC79A7",
    alpha = 0.5,
    color = NA
  ) +
  
  geom_sf(
    data = dep,
    aes(fill = "A. deppii"),
    #fill = "#56B4E9",
    alpha = 0.5,
    color = NA
  ) +
  
  geom_sf(
    data = coz,
    aes(fill = "A. cozumelus"),
    #fill= "#009E73",
    alpha = 0.7,
    color = NA,
    linewidth = 0.3
  ) +
  
  coord_sf(
    xlim = c(-110, -80),
    ylim = c(5, 30),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Species"
  ) +
  
  scale_fill_manual(
    values = GROUP_COLORS
  )
  
theme_classic()


# plot together (zoomed in)
ggplot() +
  
  # geographic background
  geom_sf(
    data = world,
    fill = "grey95",
    color = "grey60",
    linewidth = 0.2
  ) +
  
  # species distributions
  geom_sf(
    data = ang,
    aes(fill = "A. angusticeps"),
    #fill = "#CC79A7",
    alpha = 0.5,
    color = NA
  ) +
  
  geom_sf(
    data = dep,
    aes(fill = "A. deppii"),
    #fill = "#56B4E9",
    alpha = 0.5,
    color = NA
  ) +
  
  geom_sf(
    data = coz,
    aes(fill = "A. cozumelus"),
    #fill= "#009E73",
    alpha = 0.7,
    color = NA,
    linewidth = 0.3
  ) +
  
  coord_sf(
    xlim = c(-93, -86),
    ylim = c(16, 22),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Species"
  ) +
  
  scale_fill_manual(
    values = GROUP_COLORS
  )

