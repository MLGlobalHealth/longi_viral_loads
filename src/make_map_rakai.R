#################
#      AIMS     #
#################

# Plot the communities in Rakai.
# Distinguish .type with pch, pop size with size?

#################
#   Libraries   #
#################

library(data.table)
library(ggplot2)
library(patchwork)
# for maps


gitdir <- here::here()

{
    source(file.path(gitdir, 'R/paths.R'))
    source(file.path(gitdir.functions, "plotting_main_figures.R"))
}




##################
#      Main      #
##################


require(ggmap)

# Get API key required to request the map data. 
stadia_api_key <- Sys.getenv("STADIA_API_KEY")
if( stadia_api_key == "" ){
    stop(
        "Please set the STADIA_API_KEY environment variable to your API key.\n", 
        "For more info: https://docs.stadiamaps.com/authentication/#api-keys"
    )
}


# install the development version of ggmap 

library(ggmap)
try(
    ggmap::register_stadiamaps(stadia_api_key)
)
if(0){
    remove.packages("ggmap")
    devtools::install_github("stadiamaps/ggmap")
}

# us examples
{
    us <- c(left = -125, bottom = 25.75, right = -67, top = 49)
    get_stadiamap(us, zoom = 5, maptype = "stamen_toner") %>% ggmap()

    STADIA_VALID_MAP_TYPES <- c( 
        "stamen_terrain",
        "stamen_toner",
        "stamen_toner_lite",
       "stamen_watercolor",
        "stamen_terrain_background",
        "stamen_toner_background",
        "stamen_terrain_lines",
        "stamen_terrain_labels",
        "stamen_toner_lines",
        "stamen_toner_labels"
    )

    plots <- lapply(
        STADIA_VALID_MAP_TYPES,
        function(x) try({get_stadiamap(us, zoom = 5, maptype = x) |> ggmap()})
    )
    names(plots) <- STADIA_VALID_MAP_TYPES
    ggpubr::ggarrange(plotlist = plots)

    best_ones <- STADIA_VALID_MAP_TYPES[c(1, 5)]
    ggpubr::ggarrange(plotlist = plots[[best_ones]])
}


# now let us try with Rakai and the communities.
fig1a_inset <- plot.all.maps(type="inset")


p_rakai <- plot.rakai.map()
plot_all_maps(delta_inset = .72)

ggsave("tmp.png", width = 7, height = 7)
system("gthumb tmp.png")

# 1. replace get_stamenmap with `get_stadiamap`
# 2. rm source = "stamen"
# 3. maptype <- paste0("stamen_", gsub("-", "_", maptype))

