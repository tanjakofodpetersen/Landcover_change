##------------------------------------------------##
##---  CODE FOR DATA DOWNLOAD AND ANALYSIS FOR ---##
##---  "Urbanisation and land-cover change     ---##
##--- affect functional, but not compositional ---##
##---        turnover of bird communities      ---##
##------------------------------------------------##

##--- PACKAGES ---####
library(ggplot2)
library(sf)
library(sp)
library(rgdal)
library(tibble)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(raster)
library(rasterVis)
library(gridExtra)
library(rgbif)
library(stringr)
library(wicket)
library(rio)
library(purrr)
library(ggmap)
library(ggpubr)
library(mapview)
library(nlme)
library(lme4)
library(MASS)
library(iNEXT)
library(vegan)
library(NbClust)
library(factoextra)
library(scales)
library(dave)
library(car)
library(dunn.test)
library(FSA)
library(betareg)
library(rcompanion)
library(moments)
library(ggmosaic)
library(ade4)
library(adegraphics)
library(betapart)
library(spData)
library(spdep)
library(ncf)
library(mgcv)
library(geoR)
library(glmmTMB)
library(spaMM)
library(mvtnorm)
library(TH.data)
library(multcomp)
library(ggsn)
library(ggtext)
library(mapview)

##-------------####
##--- 1. LAND-COVER DATA  ---####
##--- 1.1 Load and clean up the data ---####
# Load the municipality border polygon, define CRS and convert to sf-object
Trondheim <- readOGR("ko1601admin_omr_f.shp")
proj4string(Trondheim) <- CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +vunits=m +no_defs")
Trondheim <- st_as_sf(Trondheim)
# Remove all unnecessary columns (keep only area)
Trondheim <- Trondheim[,"AREA"]

# Load the landcover maps, check/set/transform CRS:
AR5_2018 <- readOGR(dsn=path.expand("2018"), layer="ar5")  
AR5_2018 <- st_as_sf(AR5_2018)
AR5_2018 <- st_transform(AR5_2018, 32632)

AR5_2012 <- readOGR(dsn=path.expand("2012"), layer="ar5_2012")  
AR5_2012 <- st_as_sf(AR5_2012)
AR5_2012 <- st_transform(AR5_2012, 32632)

# Check/set names of the columns and do some general data cleaning - the most important here is to determine the names of land-cover categories (whether they match),
# the structure of the columns, and the latest verification dates of the data
str(AR5_2018)       # The relevant column is names "artype"
str(AR5_2012)        # The relevant column is names "arealres4" - we have to recode the new data; freshwater and ocean are not differentiated

# Recode new columns in the dataframes to have a matching "landcover"-column - I'll base this on the most conservative classification from AR5_2012
AR5_2018$landcover <- factor(ifelse(AR5_2018$artype=="81" | AR5_2018$artype=="82", paste0("80"), as.character(AR5_2018$artype)))
AR5_2012$landcover <- AR5_2012$arealres_4

# Crop the maps wit the municipality border to ensure similar extents
AR5_2012 <- st_buffer(AR5_2012, dist = 0)   # There are some issues with selfintersection - a hack to solve this is to create a zero-buffer
AR5_2012 <- st_intersection(AR5_2012, Trondheim)      
AR5_2018 <- st_intersection(AR5_2018, Trondheim)      

# Calculate the area(s) of the polygons/dataframes to check of they are somewhat similar (also compare with the Trondheim polygon)
sum(AR5_2012$area_m2 <- st_area(AR5_2012))
sum(AR5_2018$area_m2 <- st_area(AR5_2018))
st_area(Trondheim)


##--- 1.2.1 Combine and split land-cover categories ---####

# Recode new columns in the dataframes to have a matching "landcover_2"-column 
AR5_2018$landcover_3 <- factor(ifelse(AR5_2018$landcover=="11" | AR5_2018$landcover=="12", paste0("11"),
                                      ifelse(AR5_2018$landcover=="21" | AR5_2018$landcover=="22", paste0("20"),
                                             ifelse(AR5_2018$landcover=="30", paste0("30"),
                                                    as.character(AR5_2018$landcover)))))

AR5_2012$landcover_3 <- factor(ifelse(AR5_2012$landcover=="11" | AR5_2012$landcover=="12", paste0("11"),
                                      ifelse(AR5_2012$landcover=="21" | AR5_2012$landcover=="22", paste0("20"),
                                             ifelse(AR5_2012$landcover=="30", paste0("30"),
                                                    as.character(AR5_2012$landcover)))))

##--- 1.2.2 Landcover change in raster-/grid cells ---####
# Calculate and plot the area and changes in area of the new land-cover types in the grid cells - 500*500m2
# Create grid:
# Make a 500*500 grid, intersect it with the Trd-boundary, and calculate the area of each landcover type from each time period within the cells
grid_500_rev <- Trondheim %>%
  st_make_grid(cellsize = 500) %>%   # Define the side length (obs on the units)
  st_intersection(Trondheim) %>%     # Intersect with municipality
  st_cast("MULTIPOLYGON") %>%        # State object type
  st_sf() %>%
  mutate(Pixelnr = row_number())     # Add a pixelnumber/-id

# 2012 grid
{
  grid_2012_rev <- st_intersection(AR5_2012[,c("landcover_3","geometry")], grid_500_rev)       # Join/melt the data
  grid_2012_rev$area_m2 <- st_area(grid_2012_rev)                                             # Calculate the area
  grid_2012_rev <- dcast(grid_2012_rev, Pixelnr~landcover_3, fun.aggregate = sum)      # Calculate the area of each kind of landcover within each rastercell around a point
  grid_2012_rev <- merge(grid_500_rev, grid_2012_rev)      # Merge the grid and the dataframe to get the 'geometry'
  ### Check that the  grid makes sense area-wise
  range(rowSums(st_drop_geometry(grid_2012_rev[,c(2:8)])))
  grid_2012_rev$total <- rowSums(st_drop_geometry(grid_2012_rev[,c(2:8)]))
  # For some reason, fully covered grid cells are not recognised (despite the total area is 250000), and some have some missing area -
  # we will get all the grid cells if I allow for a discrepancy of 100m2 (0.04%)
  grid_2012_complete_rev <- grid_2012_rev[grid_2012_rev$total>249900,]
  }
# 2018 grid
{
  grid_2018_rev <- st_intersection(AR5_2018[,c("landcover_3","geometry")], grid_500_rev)       # Join/melt the data
  grid_2018_rev$area_m2<-st_area(grid_2018_rev)                                             # Calculate the area
  grid_2018_rev <- dcast(grid_2018_rev, Pixelnr~landcover_3, fun.aggregate = sum)      # Calculate the area of each kind of landcover within each rastercell around a point
  grid_2018_rev <- merge(grid_500_rev, grid_2018_rev)      # Merge the grid and the dataframe
  ### Check that the  grid makes sense area-wise
  range(rowSums(st_drop_geometry(grid_2018_rev[,c(2:8)])))
  grid_2018_rev$total <- rowSums(st_drop_geometry(grid_2018_rev[,c(2:8)]))
  grid_2018_complete_rev <- grid_2018_rev[grid_2018_rev$total>249900,]
}

# Calculate the differences within each kind of land-cover within the grid cells, and make sure everything add up
diff_1218_rev <- st_drop_geometry(grid_2018_complete_rev[,c(2:8)]) - st_drop_geometry(grid_2012_complete_rev[,c(2:8)])
diff_1218_rev$Pixelnr <- grid_2018_complete_rev$Pixelnr
diff_1218_rev <- merge(grid_2018_complete_rev[,c("Pixelnr","geometry")],diff_1218_rev, by="Pixelnr")

## Plots with color gradient
{
  # Make a raster with the same offset as the original grid:
  r500_rev <- SpatialGrid(GridTopology(st_bbox(Trondheim.2)[c("xmin", "ymin")], # lower left coordinate
                                       c(500,500), # pixel length and height
                                       c(round(diff(st_bbox(Trondheim.2)[c("xmin", "xmax")])/500),
                                         round(diff(st_bbox(Trondheim.2)[c("ymin", "ymax")])/500)))) # number of pixels in x and y direction
  proj4string(r500) <- CRS(st_crs(Trondheim.2)$proj4string)
  r500_rev <- raster(r500_rev)
  # 2012 - 2018
  {
    # Further, we need to convert the dataframes to rasters (as the functions uses trellis objects generated by rasterVis::levelplot)
    # Make a RasterStack of the individual land cover types - important note: when working with the rasterstack-objects,
    # R fetches them from the indicated location everytime, so do not move or rename anything while you're working!
    for (i in 2:8){
      print(i)
      
      lc <- rasterize(as(diff_1218_rev, 'Spatial'), r500_rev, field=paste("X",names(st_drop_geometry(diff_1218_rev[,i])), sep=""))
      namevec_lc<-paste('raster_1218_complete_rev/', names(st_drop_geometry(diff_1218_rev[,i])))
      writeRaster(lc, filename=namevec_lc, format='GTiff', overwrite=T)
    }
    ras_diff_1218_rev <- list.files('raster_1218_complete_rev/',full.names=T)
    ras_diff_1218_rev <- stack(ras_diff_1218_rev)
    
    # Make the trellis-objects:
    {
      p_11 <- levelplot(ras_diff_1218_rev[["X_11"]],
                        main="Change in 'Bebygd/Samferdsel' (m2), 2012-2018",
                        margin=FALSE, xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                        par.settings = list(axis.line = list(col = "transparent")))
      p_20 <- levelplot(ras_diff_1218_rev[["X_20"]],
                        main="Change in 'Dyrka' (m2), 2012-2018",
                        margin=FALSE, xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                        par.settings = list(axis.line = list(col = "transparent")))
      p_23 <- levelplot(ras_diff_1218_rev[["X_23"]],
                        main="Change in 'Innmarksbeite' (m2), 2012-2018",
                        margin=FALSE, xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                        par.settings = list(axis.line = list(col = "transparent")))
      p_30 <- levelplot(ras_diff_1218_rev[["X_30"]],
                        main="Change in 'Skog' (m2), 2012-2018",
                        margin=FALSE, xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                        par.settings = list(axis.line = list(col = "transparent")))
      p_50 <- levelplot(ras_diff_1218_rev[["X_50"]],
                        main="Change in 'Åpen' (m2), 2012-2018",
                        margin=FALSE, xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                        par.settings = list(axis.line = list(col = "transparent")))
      p_60 <- levelplot(ras_diff_1218_rev[["X_60"]],
                        main="Change in 'Myr' (m2), 2012-2018",
                        margin=FALSE, xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                        par.settings = list(axis.line = list(col = "transparent")))
      p_80 <- levelplot(ras_diff_1218_rev[["X_80"]],
                        main="Change in 'Vann' (m2), 2012-2018",
                        margin=FALSE, xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                        par.settings = list(axis.line = list(col = "transparent")))
    }
    
    # Define a function creating a plot with a diverging colour gradient centered at zero:
    {# The original function was taken from https://stackoverflow.com/questions/33750235/plotting-a-raster-with-the-color-ramp-diverging-around-zero
      # devtools::source_gist('306e4b7e69c87b1826db')
      # I have added an extra 'rev()' to reverse the 'RdBu' color gradient. OBS! For now, it only works on the named palettes from 'brewer.pal'
      
      Mydiverge0 <- function(p, ramp) {
        # p: a trellis object resulting from rasterVis::levelplot
        # ramp: the name of an RColorBrewer palette (as character), a character vector of colour names to interpolate, or a colorRampPalette.
        require(RColorBrewer)
        require(rasterVis)
        if(length(ramp)==1 && is.character(ramp) && ramp %in% 
           row.names(brewer.pal.info)) {
          ramp <- suppressWarnings(colorRampPalette(rev(brewer.pal(11, ramp))))         # Here I added the 'rev()' function to reverse the palette
        } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
          ramp <- colorRampPalette(ramp)      
        } else if(!is.function(ramp)) 
          stop('ramp should be either the name of a RColorBrewer palette, ', 
               'a vector of colours to be interpolated, or a colorRampPalette.')
        rng <- range(p$legend[[1]]$args$key$at)
        s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
        i <- findInterval(rng[which.min(abs(rng))], s)
        zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
        p$legend[[1]]$args$key$at <- s[zlim]
        p[[grep('^legend', names(p))]][[1]]$args$key$col <- ramp(1000)[zlim[-length(zlim)]]
        p$panel.args.common$col.regions <- ramp(1000)[zlim[-length(zlim)]]
        p
      }}
    
    # Make the plots
    {
      grid.arrange(Mydiverge0(p_11, ramp='RdBu'),
                   Mydiverge0(p_20, ramp='RdBu'),
                   Mydiverge0(p_23, ramp='RdBu'),
                   Mydiverge0(p_30, ramp='RdBu'),
                   Mydiverge0(p_50, ramp='RdBu'),
                   Mydiverge0(p_60, ramp='RdBu'),
                   Mydiverge0(p_80, ramp='RdBu'),
                   ncol=3, nrow=3)
    }
    
  }
  
  
}

##-----------------------------------------------####
##--- 2. DOWNLOAD GBIF DATA ---####

# OBS! The following steps are included to show how the original data download was performed - following these steps will no longer result in an
# identical dataset. The original dataset are available through the GBIF permanenet repository:
# "GBIF Occurrence Download 10.15468/dl.nxxuv6 accessed via GBIF.org on 2020-04-26"

##--- 2.1 Original data download steps - DO NOT RUN CODE! ---####

# Make R ask for you login credentials:
options(gbif_user=rstudioapi::askForPassword("my gbif username"))
options(gbif_email=rstudioapi::askForPassword("my registred gbif e-mail"))
options(gbif_pwd=rstudioapi::askForPassword("my gbif password"))

# Create a spatial filter: a bbox. Unfortunately, the municipality border is a bit too complicated for GBIF to handlle
# Transform the Trondheim bbox ("+init=epsg:4326" is the standard GBIF) and retrieve the coordinates for the polygon
st_bbox(st_transform(Trondheim, "+init=epsg:4326"))
my_wkt <- "POLYGON((10.00353 63.30279, 10.0012017 63.51659, 10.72517 63.51659, 10.72517 63.30279, 10.00353 63.30279))" 

# wicket::validate_wkt(my_wkt)    # Check that it is valid
geom_param <- paste("geometry", "within", my_wkt)

# Make a download key. NB! Maximum of 3 download requests handled simultaneously
download_key <- occ_download(
  'hasGeospatialIssue = FALSE',
  'hasCoordinate = TRUE',
  geom_param,
  type = "and"
) %>% 
  occ_download_meta

# Make the "Coffee Break"-function to retrieve the requested data
# define function
download_GBIF_API <- function(download_key,n_try,Sys.sleep_duration,destfile_name){
  start_time <- Sys.time()
  n_try_count <- 1
  
  download_url <- paste("http://api.gbif.org/v1/occurrence/download/request/",
                        download_key[1],sep="")
  
  try_download <- try(download.file(url=download_url,destfile=destfile_name,
                                    quiet=TRUE),silent = TRUE)
  
  while (inherits(try_download, "try-error") & n_try_count < n_try) {   
    Sys.sleep(Sys.sleep_duration)
    n_try_count <- n_try_count+1
    try_download <- try(download.file(url=download_url,destfile=destfile_name,
                                      quiet=TRUE),silent = TRUE)
    print(paste("trying... Download link not ready. Time elapsed (min):",
                round(as.numeric(paste(difftime(Sys.time(),start_time, units = "mins"))),2)))
  }
}


##--- 2.2 Clean and explore the data ---####
occurrence <- fread("occurrence.txt")

# Start out by removing some of the unneccessary columns:
names(occurrence)
occurrence <- occurrence[,c(1,16,57:61,64,68,70:74,80,96:105,108:113,133:136,191:199,218:232)]

# Remove records according to the following criteria:
# (1) Records with the occurrence status “absent” - we do this to ensure that we are assessing actual species richness
occurrence <- occurrence[!occurrence$occurrenceStatus=="absent",]   

# (2) Records with no registered species-level information - make a new column based on genus and species-epithet to compare.
#     We do this to ensure a fixed level of taxonomic resolution across all grid cells
occurrence <- occurrence[!occurrence$specificEpithet=="",]   # 936,176 records

# Create new species column combining the genus and species epithet:
occurrence <- unite(occurrence, species2, c("genus","specificEpithet","infraspecificEpithet"), sep = " ", remove = FALSE, na.rm = TRUE)    
occurrence$species2 <- trimws(occurrence$species2)    # trim trailing whitespace

# (3) Coordinate ucertainty below a threshold. The threshold should be adjusted according the spatial resolution we choose to work with.
#     If we go with the 500*500m, the threshold is as 354m (this also removes NA's). We do this to make sure that the data we keep can be 
#     somewhat expected to actually fall within the given grid cell
GBIF_500 <- occurrence[occurrence$coordinateUncertaintyInMeters<=354,]  # 103,446 records (upd. 2); 683,354 (upd. 3)

# (4) Only records falling within Trondheim Municipality
#     For this, we have to make the dataframe an 'sf'-object, transform the CRS, and crop it with the Trondheim.2 polygon
GBIF_500 <- st_as_sf(GBIF_500, coords = c('decimalLongitude','decimalLatitude'), crs=4326)
GBIF_500 <- st_transform(GBIF_500, crs=32632)
GBIF_500 <- st_intersection(GBIF_500, Trondheim)

##-----------------------------------####
##--- 3. COMBINE SPECIES AND LAND-COVER DATA ---####
##--- 3.1 Assign pixelnr to each record ---####
GBIF_500$ID <- c(1:nrow(GBIF_500))  # Give all records a unique ID for future reference

# Assign Pixelnr to each record for all the point datasets:
GBIF_500$Pixelnr <- as.data.frame(st_drop_geometry(st_join(GBIF_500, grid_500_rev, join = st_intersects)["Pixelnr"]))$Pixelnr

# Remove records not within the defined grid
GBIF_500 <- GBIF_500[!is.na(GBIF_500$Pixelnr),]

# As pointed out in some papers, we should include sampling effort as a covariate somehow, either by using Pixelnr as a random 
# effect or by including number of events as a covariate (or both). Calculate the latter, and save the result for later merging:
pxl <- as.data.frame(table(GBIF_500[GBIF_500$year %in% c(2010,2011,2012,2017,2018,2019),]$Pixelnr))
names(pxl) <- c("Pixelnr","N")
pxl$Pixelnr <- as.integer(as.character(pxl$Pixelnr))

# Calculate the number of records from each time-step to investigate change in sampling effort
# To account for the change in sampling effort rather than juts the sampling effort calculate the difference in N.11 and N.18,
# potentially recalculate as a proportion, and include as a predictor:
pxl.11 <- as.data.frame(table(GBIF_500[GBIF_500$year %in% c(2010,2011,2012),]$Pixelnr))
names(pxl.11) <- c("Pixelnr","N.11")
pxl.11$Pixelnr <- as.integer(as.character(pxl.11$Pixelnr))

pxl.18 <- as.data.frame(table(GBIF_500[GBIF_500$year %in% c(2017,2018,2019),]$Pixelnr))
names(pxl.18) <- c("Pixelnr","N.18")
pxl.18$Pixelnr <- as.integer(as.character(pxl.18$Pixelnr))

pxl <- full_join(full_join(pxl, pxl.11), pxl.18)
pxl[is.na(pxl)] <- 0

# Calculate proportional difference in sampling effort - obs! Anything with no sampling in 2011 will be Inf
pxl$diff <- (pxl$N.18-pxl$N.11)/pxl$N.11

##-------------------------------------------####
##--- 4. SIMILARITY INDICES FOR EACH PIXEL OVER TIME ---####
# Only retain birds and data from 2010-2012 and 2017-2019 - the temporal aspect is used to match the land-cover data:
bird.11 <- GBIF_500[GBIF_500$year>=2010 & GBIF_500$year<=2012 & GBIF_500$class=="Aves",]
bird.18 <- GBIF_500[GBIF_500$year>=2017 & GBIF_500$year<=2019 & GBIF_500$class=="Aves",]

# Plot to assess distribution - we might have to use the large grid cells 
plot(Trondheim.2)
plot(bird.11, pch=".", cex=2, col="black", add=T)
plot(bird.18, pch=".", cex=2, col="red", add=T)

# Data for community matrices
{
  ### Bird
  {
    {
      {
        rar.bird.11 <- st_drop_geometry(bird.11[,c("species2", "Pixelnr")])
        rar.bird.11$species2 <- as.factor(rar.bird.11$species2)
        rar.bird.11$Pixelnr <- as.factor(rar.bird.11$Pixelnr)
        rar.bird.11 <- droplevels(rar.bird.11)
        
        rar.bird.18 <- st_drop_geometry(bird.18[,c("species2", "Pixelnr")])
        rar.bird.18$species2 <- as.factor(rar.bird.18$species2)
        rar.bird.18$Pixelnr <- as.factor(rar.bird.18$Pixelnr)
        rar.bird.18 <- droplevels(rar.bird.18)
      }
      
      # Create empty matrices
      {
        est_bird.11 <- matrix(data=NA, ncol=nlevels(rar.bird.11$species2), nrow=nlevels(rar.bird.11$Pixelnr))
        est_bird.18 <- matrix(data=NA, ncol=nlevels(rar.bird.18$species2), nrow=nlevels(rar.bird.18$Pixelnr))
      }
      # Add column names and row names (species names and Pixelnr)
      {
        colnames(est_bird.11) <- levels(rar.bird.11$species2)
        rownames(est_bird.11) <- levels(rar.bird.11$Pixelnr)
        
        colnames(est_bird.18) <- levels(rar.bird.18$species2)
        rownames(est_bird.18) <- levels(rar.bird.18$Pixelnr)
      }
      # Make tallies of the rar.data
      {
        tallied_bird.11 <- rar.bird.11 %>%
          group_by(species2, Pixelnr) %>%
          tally()
        tallied_bird.18 <- rar.bird.18 %>%
          group_by(species2, Pixelnr) %>%
          tally()
      }
    }
  }
}
# Fill in the community matrices 
{
  ### Bird
  {
    for(r in 1:dim(est_bird.11)[1]){
      for(c in 1:dim(est_bird.11)[2]){
        print(r)
        print(c)
        est_bird.11[r,c]=ntally(i=r, j=c, Tally=tallied_bird.11, Com.matrix=est_bird.11)}
    }
    write.csv(est_bird.11, file="commat_bird_11_final.csv")
    #est_bird.11 <- read.csv("/home/ahomez/t/tanjakp/export/Landcover_change/community_matrices/taxa/commat_bird_11_final.csv", row.names = 1)
    
    for(r in 1:dim(est_bird.18)[1]){
      for(c in 1:dim(est_bird.18)[2]){
        print(r)
        print(c)
        est_bird.18[r,c]=ntally(i=r, j=c, Tally=tallied_bird.18, Com.matrix=est_bird.18)}
    }
    write.csv(est_bird.18, file="commat_bird_18_final.csv")
    #est_bird.18 <- read.csv("/home/ahomez/t/tanjakp/export/Landcover_change/community_matrices/taxa/commat_bird_18_final.csv", row.names = 1)
  }
}

##--- 4.1 Calculate dissimilarity indices ---####
# Make the individual community matrices for the dissimilarity calculations
{
  # Check the matrices - how many species and sampling events are necessary?
  {
    par(mfrow=c(1,2), mar=c(4,4,4,2))
    barplot(table(colSums(est_bird.11 > 0)), las=2, cex.names=0.75,
            xlab="No. grid cells containting species x (2010-2012)")  
    barplot(table(colSums(est_bird.18 > 0)), las=2, cex.names=0.75,
            xlab="No. grid cells containting species x (2017-2019)") 
    
    barplot(table(rowSums(est_bird.11 > 0)), las=2, cex.names=0.75,
            xlab="No. species pr. grid cell (2010-2012)") 
    barplot(table(rowSums(est_bird.18 > 0)), las=2, cex.names=0.75,
            xlab="No. species pr. grid cell (2017-2019)") 
  } 
  
  # Only keep species and pixels found in both matrices
  comm_bird_11 <- est_bird.11[rownames(est_bird.11) %in% rownames(est_bird.18) & rownames(est_bird.11) %in% diff_1218$Pixelnr,
                              colnames(est_bird.11) %in% colnames(est_bird.18)]
  comm_bird_18 <- est_bird.18[rownames(est_bird.18) %in% rownames(est_bird.11) & rownames(est_bird.18) %in% diff_1218$Pixelnr,
                              colnames(est_bird.18) %in% colnames(est_bird.11)]
  
  # Filter according to minimum number of species in each grid cell (3), and minimum number of grid cells containing species x (5)
  # (repeat this step/chunk until no more changes in matrix dimensions occur occur)
  {
    comm_bird_11 <- comm_bird_11[rowSums(comm_bird_11 > 0) >= 3 ,     # Number of species - at least 3 species in each grid cell
                                 colSums(comm_bird_11 > 0) >= 5]     # Number of pixels containing species x
    comm_bird_18 <- comm_bird_18[rowSums(comm_bird_18 > 0) >= 3  ,     # Number of species - at least 3 species in each grid cell
                                 colSums(comm_bird_18 > 0) >= 5]     # Number of pixels containing species x
    
    # Only keep species and pixels found in both matrices
    comm_bird_11 <- comm_bird_11[rownames(comm_bird_11) %in% rownames(comm_bird_18), colnames(comm_bird_11) %in% colnames(comm_bird_18)]
    comm_bird_18 <- comm_bird_18[rownames(comm_bird_18) %in% rownames(comm_bird_11), colnames(comm_bird_18) %in% colnames(comm_bird_11)]
  }
  
  # Combine into one community matrix:
  comm_bird <- rbind(comm_bird_11, comm_bird_18)
  comm_bird[comm_bird > 0] <- 1      # Make it presence/absence 
  rownames(comm_bird) <- c(paste(rownames(comm_bird[c(1:(nrow(comm_bird)/2)),]), "11", sep="."),
                           paste(rownames(comm_bird[c(((nrow(comm_bird)/2)+1):(nrow(comm_bird))),]), "18", sep="."))
}

# Calculate dissimilarity indices
{
  dist.moved_bird_rev <- data.frame(Pixelnr = as.integer(as.character(comm_bird_dt[c(1:(nrow(distdf_bird)/2)),"Pixelnr"])),
                                    dist_11.18 = NA)
  for(i in 1:(nrow(dist.moved_bird_rev))){   
    dist.moved_bird_rev$dist_11.18[i] <- distdf_bird[i+nrow(dist.moved_bird_rev) , i]
  }
  
  df_bird_rev <- left_join(dist.moved_bird_rev, st_drop_geometry(diff_1218_rev))
}


##--- 4.2 Partitioning of beta diversity ---####
# Partioning the beta diversity ((dis-)similarity) into its turnover- and nestedness components:
# Patterns in nestedness are likely caused by the difference in sampling effort, and the effects caused strictly by turnover are likely what we wish
# want to investigate. Based on the papers and methods suggested by Baselga (2010-2012), these can be partitioned.

# Convert the matrices to presence/absence
comm_bird_11.pa <- comm_bird_11
comm_bird_11.pa[comm_bird_11.pa > 0] <- 1
comm_bird_18.pa <- comm_bird_18
comm_bird_18.pa[comm_bird_18.pa > 0] <- 1

beta.bird <- beta.temp(comm_bird_11.pa, comm_bird_18.pa, index.family = "jaccard")
shared_species <- data.frame(Pixelnr = rownames(comm_bird_18.pa), shared_species = NA)
for(i in 1:nrow(shared_species)){
  names1 <- colnames(comm_bird_11.pa[i,colSums(comm_bird_11.pa[i,])>0])  # species present 2011
  names2 <- colnames(comm_bird_18.pa[i,colSums(comm_bird_18.pa[i,])>0])  # species present 2018
  shared_species[i,"shared_species"] <- length(names1[names1 %in% names2])
}

# Plot distributions of components:
par(mfrow=c(1,1))
plot(density(beta.bird$beta.jac), xlim=c(0,1), ylim=c(0, 3.5), xlab=expression(beta), main='', lwd=3, col="gray")
lines(density(beta.bird$beta.jtu),lty=1, lwd=2)
lines(density(beta.bird$beta.jne),lty=2,lwd=2)
legend("topright", legend=c(expression(beta[Jaccard]), expression(beta[turnover]), expression(beta[nestedness])),
       col = c("gray","black","black"), lty=c(1,1,2), lwd=c(3,2,2), )
# Most of the (dis-)similarity is caused by the turnover, only a relatively small part is due to nestedness.
# Another way of assessing this can be percentage-wise - how great a proportion of the total beta-dissimilarity is generally caused by each component:
boxplot(beta.bird$beta.jtu / beta.bird$beta.jac)
boxplot(beta.bird$beta.jne / beta.bird$beta.jac)


##--- 4.3 Visualise and test differences in beta-components ---####
# Add the partitioned beta-diversity components to the dataframe
df_bird_beta_rev <- full_join(beta.bird %>%
                                mutate(Pixelnr = as.integer(rownames(beta.bird))),
                              df_bird_rev[,c("Pixelnr","11","20","23","30","50","60","80")],
                              by='Pixelnr')

# The need for beta-regression (as the response is bounded between 0 and 1) is causing problems, if we also have to take into account spatial autocorrelation.
# The chosen solution to the issue is to logit-transform the variable, and stick to standard lm:
df_bird_beta_rev$beta_logit <- car::logit(df_bird_beta_rev$beta.jtu)


# Plot the different beta diversity components by land-cover change (one plot for each type)
ggarrange(
  ggplot(data=df_bird_beta_rev, aes(x=`11`, y=beta.jac))  +
    geom_point() +
    geom_smooth( method = "loess") +
    labs(x="11", y=expression(beta[total])),
  ggplot(data=df_bird_beta_rev, aes(x=`11`, y=beta.jtu))  +
    geom_point() +
    geom_smooth( method = "loess") +
    labs(x="11", y=expression(beta[turnover])),
  ggplot(data=df_bird_beta_rev, aes(x=`11`, y=beta.jne))  +
    geom_point() +
    geom_smooth( method = "loess") +
    labs(x="11", y=expression(beta[nestedness])),
  ggplot(data=df_bird_beta_rev, aes(x=`11`, y=beta_logit))  +
    geom_point() +
    geom_smooth( method = "loess") +
    labs(x="11", y=expression(beta[logit])),
  
  nrow=2, ncol=2)

# Plot maps of the variables
# Beta-logit
ggplot() + geom_sf(data = diff_1218_rev[!is.na(diff_1218_rev$beta_logit),], aes(fill=beta_logit), color="gray80", size=0.25) +
  labs(x="", y="") + 
  scale_fill_gradient(low="#fff5f0",high="#cb181d", name=expression("logit("*beta[turnover]*")")) +
  theme_minimal() +
  scalebar(data=diff_1218_rev, location = "bottomright", dist = 2, dist_unit = "km", transform=FALSE, st_height = 0.01, st.size = 2, border.size = .1)

##----------------------------------####
##--- 5. MODEL BETA_TURNOVER AS A FUNCTION OF LAND-COVER CHANGE ---####
##--- 5.1 Non-spatial model ---####
# Preliminary model of beta_turnover as a function of change within the different land-covers. First, we need to add the the sampling
# effort as a covariate to the dataframe:
df_bird_beta_rev <- left_join(df_bird_beta_rev, pxl)

# Try making a model including all the land-cover categories, and sampling effort as an covariate (not offset!)
# (we cannot use Pixelnr as a random effect when it only occurs once)
model_betalc_rev <- lm(beta_logit ~ `11` + `20` + `23` + `30` + `60` + `50` + `80` + N + diff,
                       data = df_bird_beta_rev)
# Do model selection to see if any of the land-covers can be removed - in this case, deltaAIC<2, so remove some manually:
step(model_betalc_rev)  

# A problem is that we do not get a hint as to which one to remove - continue with all

# Test the data and the model for Spatial Autocorrelation
# Add the dissimilarity index to the sf-data-frame for easier plotting
diff_1218_plot_rev <- right_join(diff_1218_rev, df_bird_beta_rev[,c("Pixelnr","11","20","50","beta.jtu","N","N.11","N.18","diff")])

# To get the correct coordinates, first get the centroids of each grid cell, then retrieve the xy-coordinates.
xy_beta_rev <- st_coordinates(st_centroid(diff_1218_plot_rev)) 

clust.nb_lc_rev <- dnearneigh(as.matrix(xy_beta_rev[,1:2]), 249, 708) # Find the neighbors - give lower and upper distance class here
# OBS! The classes are in euclidian distance (m), thus we need a reasonable distance to define a neighbouring grid cell. 
# First order neighbours: the upper distance to length(diagonal) = sqrt((500^2)+(500^2))
clust.listw_lc_rev <- nb2listw(clust.nb_lc_rev, zero.policy = T)       # Turns neighbourhood object into a weighted list

# Make a correlogram:
correlog_rev <- correlog(xy_beta_rev[,1], xy_beta_rev[,2], residuals(model_betalc_rev), na.rm = T, increment = 1, resamp = 0, latlon=T)

# Plot the first 20 distance classes
par(mfrow=c(1,1))
par(mar=c(5,5,0.1, 0.1))
plot(correlog_rev$correlation[1:20], type="b", pch=16, lwd=1.5, xaxt="n",
     xlab="distance (m)", ylab="Moran's I");
axis(side=1, at=seq(1,20, by=2), labels=c(500,1000,1415,1582,2000,2122,2500,2693,2916,3042));
abline(h=0); abline(v=2, lty=2, col="red")

# Make a map of the residuals:
plot(xy_beta_rev[,1], xy_beta_rev[,2], col=c("blue", "red")[sign(resid(model_betalc_rev))/2+1.5], pch=19,
     cex=abs(resid(model_betalc_rev))/max(resid(model_betalc_rev))*2, xlab="geographical x- coordinates", ylab="geographical y-coordinates")

# calculate Moran's I values explicitly for a certain distance, and to test for its significance:
GlobMT_rev <- moran.test(residuals(model_betalc_rev), listw=clust.listw_lc_rev, zero.policy = T)
GlobMT_rev     

# Look at it through Monte-Carlo simulation as well:
MC_rev <- moran.mc(residuals(model_betalc_rev), clust.listw_lc_rev,
                   zero.policy = TRUE, nsim=999)   
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
plot(MC_rev, main="")     
abline(v=MC_rev$statistic, lty=2, col="red")

##--- 5.2 Spatially explicit model ---####
### Add longititude and latitude to the dataframe
xy_all_rev <- data.frame(Pixelnr = diff_1218_rev$Pixelnr,
                         longitude = st_coordinates(st_centroid(diff_1218_rev))[,"X"],
                         latitude = st_coordinates(st_centroid(diff_1218_rev))[,"Y"])
df_bird_beta_rev <- left_join(df_bird_beta_rev, xy_all_rev)
# Add a renamed column ("Developed") to make predictions possible
df_bird_beta_rev$Developed <- df_bird_beta_rev$`11`

# Join with the spatial dataframe and plot/map:
outplot <- merge(diff_1218_rev, df_bird_beta_rev[,c(4,13:16,19)], all=F)
# Do some data exploration to see if we need to exclude outliers:
Mydotplot(st_drop_geometry(outplot[,c("beta_logit","N","N.11","N.18","diff","Developed")]))  # we ABOSULTELY have at least one outlier
plot(outplot[,c(9:14)])

mapview(outplot, zcol="N") +
  mapview(outplot, zcol="N.11") +
  mapview(outplot, zcol="N.18") +
  mapview(outplot, zcol="diff")

# The outlier regarding N, N.11 and N.18 is Leangenbukta
# Remove the two outliers, as they might have undue influence on the results
df_bird_beta_rev <- df_bird_beta_rev[df_bird_beta_rev$N < 30000 & df_bird_beta_rev$diff < 100, ]

#---------------------------------####

# Include coordinates as a random effect to account for spatial autocorrelation:
# Only include Developed area:
spaMM_mat11_rev <- step(fitme(beta_logit ~ N + diff + Developed + Matern(1|longitude+latitude),  
                              data = df_bird_beta_rev, verbose = c(trace=TRUE), method = "ML")) 
# According to this, we should not include 'diff' as a predictor
# Based on deltaAIC, we should also exclude Developed, giving a model only including sampling effort (N)

# Only (difference in) sampling effort and SAC
spaMM_mat_null1_rev <- step(fitme(beta_logit ~ N + diff + Matern(1|longitude+latitude), data = df_bird_beta_rev, verbose = c(trace=TRUE), method = "ML"))
# In this case, we should once again take out diff

# Only SAC
spaMM_mat_null2_rev <- fitme(beta_logit ~ 1 + Matern(1|longitude+latitude), data = df_bird_beta_rev, verbose = c(trace=TRUE), method = "ML")

# Compare AICs
AIC(spaMM_mat11_rev)
AIC(spaMM_mat_null1_rev)
AIC(spaMM_mat_null2_rev)

summary(spaMM_mat11_rev)
summary(spaMM_mat_null1_rev)
summary(spaMM_mat_null2_rev)

# Check residuals
par(mfrow=c(3,2))
plot(x = fitted(spaMM_mat11_rev), y = resid(spaMM_mat11_rev, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson residuals", main="Model", cex.lab = 1.5)   
plotNormalHistogram(resid(spaMM_mat11_rev, type = "pearson"))   

plot(x = fitted(spaMM_mat_null1_rev), y = resid(spaMM_mat_null1_rev, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson residuals", main="Null model", cex.lab = 1.5)   
plotNormalHistogram(resid(spaMM_mat_null1_rev, type = "pearson"))  

plot(x = fitted(spaMM_mat_null2_rev), y = resid(spaMM_mat_null2_rev, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson residuals", main="Null model", cex.lab = 1.5)   
plotNormalHistogram(resid(spaMM_mat_null2_rev, type = "pearson")) 

# Check spatial effects/correlation for all models:
par(mfrow=c(2,2))
# Incl. predictors
plot(as.numeric(dist(df_bird_beta_rev[,c("longitude","latitude")])),  # Distance matrix
     as.numeric(MaternCorr(dist(df_bird_beta_rev[,c("longitude","latitude")]), nu = 0.103584319, rho = 0.000847982)),  # Matérn correlation, use nu and rho from the model summary
     xlab = "Distance between pairs of location [in m]", ylab = "Estimated correlation",
     xlim=c(0,5000)); abline(h=0.1, col="red", lty=2)   # ca. 1000 m before the correlaion goes below 0.1

# Incl. one predictor
plot(as.numeric(dist(df_bird_beta_rev[,c("longitude","latitude")])),  # Distance matrix
     as.numeric(MaternCorr(dist(df_bird_beta_rev[,c("longitude","latitude")]), nu = 0.1000900028, rho = 0.0009182647)),  # Matérn correlation, use nu and rho from the model summary
     xlab = "Distance between pairs of location [in m]", ylab = "Estimated correlation",
     xlim=c(0,5000)); abline(h=0.1, col="red", lty=2)   # ca. 1000 m before the correlaion goes below 0.1

# Only sampling effort
plot(as.numeric(dist(df_bird_beta[,c("longitude","latitude")])),  # Distance matrix
     as.numeric(MaternCorr(dist(df_bird_beta[,c("longitude","latitude")]), nu = 0.1058279603, rho = 0.0008758968)),  # Matérn correlation, use nu and rho from the model summary
     xlab = "Distance between pairs of location [in m]", ylab = "Estimated correlation",
     xlim=c(0,5000)); abline(h=0.1, col="red", lty=2)   # ca. 1000 m before the correlaion goes below 0.1

# Only SAC
plot(as.numeric(dist(df_bird_beta[,c("longitude","latitude")])),  # Distance matrix
     as.numeric(MaternCorr(dist(df_bird_beta[,c("longitude","latitude")]), nu = 0.0858025893, rho = 0.0007125543)),  # Matérn correlation, use nu and rho from the model summary
     xlab = "Distance between pairs of location [in m]", ylab = "Estimated correlation",
     xlim=c(0,5000)); abline(h=0.1, col="red", lty=2)   # ca. 1000 m before the correlaion goes below 0.1

##----------------------------------####
##--- 6. TRAITS ---####
##--- 6.1 Functional traits of birds ---####
# Investigate the potential patterns/differences in traits between the birds. To test this, we might have to group the birds according to their
# respnses to land-cover change, e.g. whether the probablity of "lost" increases or decreases with increasing land-cover change.

# Read the trait files - we hve multiple, and they should thus be combined into one:
traits1 <- read.csv2("birdnames_multiple.csv")         # Self-coded; habitat based on descriptions from Birds of the World Online, and conservaton status for Norway
traits2 <- read.csv("/Amniote_Database_Aug_2015.csv")   # Myhrvold et al. (2016)
traits3 <- fread("BirdFuncDat.txt")                    # Wilman et al. (2014)

# Clean-up and combining the dataframes
traits1$species <- as.character(traits1$species)
traits1[traits1 == "Columba.livia.domestica"] <- "Columba.livia"
traits2$species2 <- paste(traits2$genus, traits2$species, sep=".")
traits3$species <- gsub(" ",".",traits3$Scientific)
traits3$ForStratCat <- ifelse((traits3$`ForStrat-watbelowsurf`>=50) | (traits3$`ForStrat-watbelowsurf`>=40 & traits3$`ForStrat-wataroundsurf`>=30 & traits3$`ForStrat-wataroundsurf` < 50), paste0("watbelowsurf"),
                              ifelse(traits3$`ForStrat-wataroundsurf`>=50 | (traits3$`ForStrat-wataroundsurf`>=40 & traits3$`ForStrat-ground`>=40), paste0("wataroundsurf"),
                                     ifelse((traits3$`ForStrat-ground`>=50 & traits3$`ForStrat-wataroundsurf`<50) | (traits3$`ForStrat-ground`>=40 & traits3$`ForStrat-understory`>=20) , paste0("ground"),
                                            ifelse((traits3$`ForStrat-understory`>=50 & traits3$`ForStrat-ground`<50) | (traits3$`ForStrat-understory`>=40 & traits3$`ForStrat-midhigh`<=40 & traits3$`ForStrat-ground`<40), paste0("understory"),
                                                   ifelse((traits3$`ForStrat-midhigh`>=50 & traits3$`ForStrat-understory`<50) | (traits3$`ForStrat-midhigh`>=40 & traits3$`ForStrat-understory`<40) |
                                                            (traits3$`ForStrat-midhigh`==33 & traits3$`ForStrat-understory`==33) | (traits3$`ForStrat-midhigh`==30 & traits3$`ForStrat-canopy`==30) | traits3$`ForStrat-midhigh`==25, paste0("midhigh"),
                                                          ifelse((traits3$`ForStrat-canopy`>=50) | (traits3$`ForStrat-canopy`==33 & traits3$`ForStrat-aerial`==33), paste0("canopy"),
                                                                 ifelse(traits3$`ForStrat-aerial`>=50 | (traits3$`ForStrat-aerial`==20 & traits3$`ForStrat-canopy`==20), paste0("aerial"), NA)))))))

traits_total <- left_join(left_join(traits1[,],
                                    traits2[,c("species2","litter_or_clutch_size_n","litters_or_clutches_per_y","adult_body_mass_g",
                                               "egg_mass_g","incubation_d","longevity_y")], by=c("species" = "species2")),
                          traits3[,c("species","Diet-5Cat","ForStratCat","Nocturnal","BodyMass-Value")], by=c("synonym" = "species"))
traits_total[traits_total == -999] <- NA


##--- 6.2 Specify which species were "lost" and which were "gained" ---####

# Create a list of dataframes with full species lists and -number for each pixel:
spec_bird_rev <- vector(mode = "list", length = nrow(df_bird_beta_rev))     # Empty list with an entry for each pixel
names(spec_bird_rev) <- df_bird_beta_rev$Pixelnr

# Add a community matrix to each list entry (species as columns, years as rows) - only include the species actually present
for(i in 1:nrow(df_bird_beta_rev)){
  print(i)
  spec_bird_rev[[i]] <- comm_bird[c(i, i+nrow(df_bird_beta_rev)), colSums(comm_bird[c(i, i+nrow(df_bird_beta_rev)),]) > 0]
}

# Make lists of the species lost and gained from each pixel:
lost_bird_rev <- list()
for(i in 1:nrow(df_bird_beta_rev)){
  print(i)
  lost_bird_rev[[i]] <- colnames((spec_bird_rev[[i]][, (spec_bird_rev[[i]][1, ] > 0) & (spec_bird_rev[[i]][2, ] == 0), drop=FALSE])) # 'drop=F' is needed - otherwise colnames wll be dropped for pixels with only on appropriate species in either of the periods
}
names(lost_bird_rev) <- df_bird_beta_rev$Pixelnr

gained_bird_rev <- list()
for(i in 1:nrow(df_bird_beta_rev)){
  print(i)
  gained_bird_rev[[i]] <- colnames((spec_bird_rev[[i]][, (spec_bird_rev[[i]][2, ] > 0) & (spec_bird_rev[[i]][1, ] == 0), drop=FALSE]))
}
names(gained_bird_rev) <- df_bird_beta_rev$Pixelnr

remained_bird_rev <- list()  # Make a list of the species for which nothing changed - this list is needed later
for(i in 1:nrow(df_bird_beta_rev)){
  print(i)
  remained_bird_rev[[i]] <- colnames((spec_bird_rev[[i]][, (spec_bird_rev[[i]][1, ] > 0) & (spec_bird_rev[[i]][2, ] > 0), drop=FALSE])) # 'drop=F' is needed - otherwise colnames wll be dropped for pixels with only on appropriate species in either of the periods
}
names(remained_bird_rev) <- df_bird_beta_rev$Pixelnr

# Make lists of the species lost and gained from each pixel. Use the list of 'gained_bird' and 'lost_bird' to get the species lists.
# Create a dataframe (one for 'loss' and one for 'gain') with Pixelnr and land-cover dissimilarity, and a column for each species;
lost_bird.df_rev <- cbind(df_bird_beta_rev[,c("Pixelnr","beta.jtu","beta.jne","beta.jac","beta_logit",
                                              "11","20","23","30","50","60","80","N","N.11","N.18","diff","longitude","latitude")],
                          as.data.frame(matrix(nrow=nrow(df_bird_beta_rev), ncol=length(sort(unique(unlist(lost_bird_rev, use.names = F)))))))
names(lost_bird.df_rev) <- c("Pixelnr","beta.jtu","beta.jne","beta.jac","beta_logit","11","20","23","30","50","60","80","N","N.11","N.18","diff","longitude","latitude",
                             sort(unique(unlist(lost_bird_rev, use.names = F))))
for(i in 1:nrow(lost_bird.df_rev)){
  for(j in 19:ncol(lost_bird.df_rev)){
    lost_bird.df_rev[i,j] <- ifelse(colnames(lost_bird.df_rev[j]) %in% lost_bird_rev[[i]], 1,
                                    ifelse(colnames(lost_bird.df_rev[j]) %in% remained_bird_rev[[i]], 0, NA))
  }
}      # '1' means that the species was lost, '0' means it remained, 'NA' means it was never there

gained_bird.df_rev <- cbind(df_bird_beta_rev[,c("Pixelnr","beta.jtu","beta.jne","beta.jac","beta_logit",
                                                "11","20","23","30","50","60","80","N","N.11","N.18","diff","longitude","latitude")],
                            as.data.frame(matrix(nrow=nrow(df_bird_beta_rev), ncol=length(sort(unique(unlist(gained_bird_rev, use.names = F)))))))
names(gained_bird.df_rev) <- c("Pixelnr","beta.jtu","beta.jne","beta.jac","beta_logit",
                               "11","20","23","30","50","60","80","N","N.11","N.18","diff","longitude","latitude", sort(unique(unlist(gained_bird_rev, use.names = F))))
for(i in 1:nrow(gained_bird.df_rev)){
  for(j in 19:ncol(gained_bird.df_rev)){
    gained_bird.df_rev[i,j] <- ifelse(colnames(gained_bird.df_rev[j]) %in% gained_bird_rev[[i]], 1, 
                                      ifelse(colnames(gained_bird.df_rev[j]) %in% remained_bird_rev[[i]], 0, NA))
  }
}      # '1' means that the species has been gained, '0' means that nothing happened, 'NA' means it was never there

##--- 6.3 Combine species lists with trait data ---####

# Rather than investigating the individual species responses,  group the species into functional groups and/or add their traits as covariates
# in models of (dis-)appearance.  As an initial point, add all of the traits, and do backwards model selection.
# Create long-format dataframes with 'species' as a covariate, and join that with the trait data:
df_bird_trait.lost_rev <- dplyr::left_join(melt(lost_bird.df_rev, id=c("Pixelnr","beta.jtu","beta.jne","beta.jac","beta_logit",
                                                                       "11","20","23","30","50","60","80","N","N.11","N.18","diff","longitude","latitude"),
                                                measure = names(lost_bird.df_rev)[19:169],value.name = "lost", variable.name = "species"),
                                           traits_total, by=c("species"))     
# '1' means that the species was lost, '0' means it remained, 'NA' means it was never there

df_bird_trait.gained_rev <- dplyr::left_join(melt(gained_bird.df_rev, id=c("Pixelnr","beta.jtu","beta.jne","beta.jac","beta_logit",
                                                                           "11","20","23","30","50","60","80","N","N.11","N.18","diff","longitude","latitude"),
                                                  measure = names(gained_bird.df_rev)[19:167],value.name = "gained", variable.name = "species"),
                                             traits_total, by=c("species"))     
# '1' means that the species has been gained, '0' means that nothing happened (remained), 'NA' means it was never there

### OBS! Using the 'gained' dataframe fom above makes the interprtations of the the following models strange - we should rather "swap" the NA and 0 designation;
### thus, species  categorised as '0' means species that remained absent from the grid cells. Species remaining in the grid cells tell us nothing!
df_bird_trait.gained_rev$gained[df_bird_trait.gained_rev$gained==0] <- "remained"  # Extra step to avoid confusion
df_bird_trait.gained_rev$gained[is.na(df_bird_trait.gained_rev$gained)] <- 0
df_bird_trait.gained_rev$gained[df_bird_trait.gained_rev$gained=="remained"] <- NA
df_bird_trait.gained_rev$gained <- as.integer(df_bird_trait.gained_rev$gained)

# We should combine "coastal" and "wetland"
df_bird_trait.lost_rev[df_bird_trait.lost_rev=="coastal"] <- "wetland"
df_bird_trait.gained_rev[df_bird_trait.gained_rev=="coastal"] <- "wetland"
df_bird_trait.lost_rev <- droplevels(df_bird_trait.lost_rev)
df_bird_trait.gained_rev <- droplevels(df_bird_trait.gained_rev)

# Relevel the factorial levels to have a more reasonable baseline
df_bird_trait.lost_rev$habitat.f <- relevel(df_bird_trait.lost_rev$habitat, ref="urban")
df_bird_trait.lost_rev$DietCat.f <- relevel(as.factor(df_bird_trait.lost_rev$`Diet-5Cat`), ref="PlantSeed")
df_bird_trait.lost_rev$ForStratCat.f <- relevel(as.factor(df_bird_trait.lost_rev$ForStratCat), ref="watbelowsurf")

df_bird_trait.gained_rev$habitat.f <- relevel(df_bird_trait.gained_rev$habitat, ref="urban")
df_bird_trait.gained_rev$DietCat.f <- relevel(as.factor(df_bird_trait.gained_rev$`Diet-5Cat`), ref="PlantSeed")
df_bird_trait.gained_rev$ForStratCat.f <- relevel(as.factor(df_bird_trait.gained_rev$ForStratCat), ref="watbelowsurf")

# Leave out Cygnus cygnus as it is an outlier regarding body mass:
df_bird_trait.lost_rev <- df_bird_trait.lost_rev[!(df_bird_trait.lost_rev$species=="Cygnus.cygnus"),]
df_bird_trait.gained_rev <- df_bird_trait.gained_rev[!(df_bird_trait.gained_rev$species=="Cygnus.cygnus"),]

##------------------------------------------####
##--- 7. MODELS OF FUNCTIONAL GROUPS ---####
##--- 7.1 Traits- species loss  ---####

# We have issues with the variable "11" - I here rename it. Also add "Developed" as a single variable
df_lost_rev <- droplevels(df_bird_trait.lost_rev[!(is.na(df_bird_trait.lost_rev$lost) | is.na(df_bird_trait.lost_rev$DietCat.f)),
                                                 c("Pixelnr","11","N","diff","longitude","latitude","adult_body_mass_g","longevity_y","habitat.f","DietCat.f","ForStratCat.f","lost")])
names(df_lost_rev) <- c("Pixelnr","Developed","N","diff","longitude","latitude","adult_body_mass_g","longevity_y","habitat.f","DietCat.f","ForStratCat.f","lost")

# Which Pixelnr in the trait dataframes are not in the beta-dataframes:
df_lost_rev[!(df_lost_rev$Pixelnr %in% df_bird_beta_rev$Pixelnr),]$Pixelnr  # The two outliers
# Remove the outliers:
df_lost_rev <- df_lost_rev[(df_lost_rev$Pixelnr %in% df_bird_beta_rev$Pixelnr),]

# Construct actual models; do stepwise backwards model selection simultaneously
spaMM_lost <- step( fitme(lost ~ Developed + Developed*habitat.f +  Developed*DietCat.f +  Developed*ForStratCat.f +
                            adult_body_mass_g + longevity_y + N + diff +
                            (1 | Pixelnr) + Matern(1 | longitude+latitude),
                          family = binomial, method = "ML", verbose=c(trace=T),
                          data=df_lost_rev))
# Remove extra variable (deltaAIC<2 ; adult_body_mass_g)
spaMM_lost_final <- fitme(lost ~ Developed + Developed*habitat.f +  Developed*DietCat.f +  ForStratCat.f +
                            longevity_y + N + diff +
                            (1 | Pixelnr) + Matern(1 | longitude+latitude),
                          family = binomial, method = "ML", verbose=c(trace=T),
                          data=df_lost_rev)

# Null - only spatial dependency
spaMM_lost_null1 <- fitme(lost ~ 1 + (1 | Pixelnr) + Matern(1 | longitude+latitude),
                          family = binomial, method = "ML", verbose=c(trace=T),
                          data=df_lost_rev)
# Null - only N and diff
spaMM_lost_null2 <- fitme(lost ~ N + diff +
                            (1 | Pixelnr) + Matern(1 | longitude+latitude),
                          family = binomial, method = "ML", verbose=c(trace=T),
                          data=df_lost_rev)

AIC(spaMM_lost)
AIC(spaMM_lost_final)
AIC(spaMM_lost_null1)
AIC(spaMM_lost_null2)

summary(glht(spaMM_lost_final, mcp("ForStratCat.f" = "Tukey"), coef.=fixef.HLfit))

# Check spatial effects/correlation
plot(as.numeric(dist(df_lost_rev[,c("longitude","latitude")])),  # Distance matrix
     as.numeric(MaternCorr(dist(df_lost_rev[,c("longitude","latitude")]), nu = 0.1372550636, rho = 0.0009999514)),  # Matérn correlation, use nu and rho from the model summary
     xlab = "Distance between pairs of location [in m]", ylab = "Estimated correlation",
     xlim=c(0,5000)); abline(h=0.1, col="red", lty=2)   # ca. 1000 m before the correlaion goes below 0.1

# Check numbers more closely:
View(data.frame(dist = c(seq(from=1090, to=1095 , by=0.5)),
                corr = as.numeric(MaternCorr(c(seq(from=1090, to=1095 , by=0.5)), nu = 0.1372550636, rho = 0.0009999514))))

##--- 7.1.1 Model predictions ---####
# To be able to better visualise and interpret the models, create a hypothetical dataframe and make predictions incl.
# confidence intervals. I will exclude the spatial effects from the predictions, and thus only assess the effects of the other
# variables. When not assessed directly, I will leave each variable at the the baseline levels (factor), or the levels perviously
# specified as the mean (numerical). I will insert a constant longitude and latitude (a set of coordinates within Trondheim;
# median of all coordinates). As I am mainly interested in the differences in responses based on the interaction, I will need to
# let both variable fluctuate simultaneously

# Linear effects: longevity, forage stratum, N and diff
{
  ## Rather than only using one baseline level of the other non-focal variables, try and include all other levels as well
  ## to asses the impact of the varying intercepts (and what that will mean for the numbers used later as baselines) 
  # Sampling effort N
  set.seed(123)
  newdata_lost_N <- {data.frame(longitude=569930.1, latitude=7031341,
                                Developed = mean(df_bird_beta_rev$Developed), diff=0, longevity = 15,     # Be aware of the mean(Developed) - this should be from the 'df_bird_beta_rev' dataframe to not have duplicates  
                                N = c(rep(seq(from=min(df_lost_rev$N), to=max(df_lost_rev$N), length.out = 25), 4*7*9)),
                                DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                   levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf"), each=25*9), 4),
                                                       levels=c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf")),
                                habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                   levels=c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland")) )
  }
  names(newdata_lost_N) <- c("longitude","latitude","Developed","diff","longevity_y","N","DietCat.f","ForStratCat.f","habitat.f" )
  newdata_lost_N$lost <- as.numeric(predict(spaMM_lost_final, newdata_lost_N, re.form = NA, type="response"))
  newdata_lost_N <- cbind(newdata_lost_N,
                          get_intervals(spaMM_lost_final, newdata = newdata_lost_N, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_lost_N <- {newdata_lost_N %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Difference in sampling effort, diff
  set.seed(123)
  newdata_lost_diff <- {data.frame(longitude=569930.1, latitude=7031341, Developed = mean(df_bird_beta_rev$Developed),longevity_y = 15, N=10000,
                                   diff = c(rep(seq(from=min(df_lost_rev$diff), to=max(df_lost_rev$diff), length.out = 25), 4*7*9)),
                                   DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                      levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                   ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                    "understory","wataroundsurf"), each=25*9), 4),
                                                          levels=c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                   "understory","wataroundsurf")),
                                   habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet",
                                                                    "open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                      levels=c("urban","generalist","marine","open","open_wet",
                                                               "open_woodland","scrub","water","woodland")) )
  }
  newdata_lost_diff$lost <- as.numeric(predict(spaMM_lost_final, newdata_lost_diff, re.form = NA, type="response"))
  newdata_lost_diff <- cbind(newdata_lost_diff,
                             get_intervals(spaMM_lost_final, newdata = newdata_lost_diff, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_lost_diff <- {newdata_lost_diff %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Longevity
  set.seed(123)
  newdata_lost_longevity <- {data.frame(longitude=569930.1, latitude=7031341, Developed = mean(df_bird_beta_rev$Developed), diff=0, N=10000, 
                                        longevity_y = c(rep(seq(from=min(df_lost_rev$longevity_y),
                                                                to=max(df_lost_rev$longevity_y), length.out = 25), 4*7*9)),
                                        DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                           levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                        ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                         "understory","wataroundsurf"), each=25*9), 4),
                                                               levels=c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                        "understory","wataroundsurf")),
                                        habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet",
                                                                         "open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                           levels=c("urban","generalist","marine","open","open_wet",
                                                                    "open_woodland","scrub","water","woodland")) )
  }
  newdata_lost_longevity$lost <- as.numeric(predict(spaMM_lost_final, newdata_lost_longevity, re.form = NA, type="response"))
  newdata_lost_longevity <- cbind(newdata_lost_longevity,
                                  get_intervals(spaMM_lost_final, newdata = newdata_lost_longevity, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_lost_longevity <- { newdata_lost_longevity %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Plot:
  ggsave("Longevity_lost.pdf",
         ggplot(newdata_lost_longevity, aes(x = longevity_y, y = lost, color=DietCat.f)) +
           geom_path(aes(color=DietCat.f),size=0.5) +
           geom_ribbon(aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975, fill=DietCat.f), alpha = 0.2, linetype=0) +
           geom_point(data=df_lost_rev,
                      size=0.2, shape=4, color="gray60")  +
           theme(legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Longevity (years)", y="Probability of disappearance") +
           scale_fill_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           scale_color_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           guides(fill=guide_legend(nrow=2)) +
           facet_grid(habitat.f ~ ForStratCat.f,
                      labeller = labeller(habitat.f=habitat.labs, ForStratCat.f=forstrat.labs)) +
           coord_cartesian(ylim = c(0,0.4)),   
         width = 16.6, height=20, units="cm", dpi=350)
  
  ggsave("N_lost.pdf",
         ggplot(newdata_lost_N, aes(x = N, y = lost, color=DietCat.f)) +
           geom_path(aes(color=DietCat.f),size=0.5) +
           geom_ribbon(aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975, fill=DietCat.f), alpha = 0.2, linetype=0) +
           geom_point(data=df_lost_rev,
                      size=0.2, shape=4, color="gray60")  +
           theme(legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Total number of observations (N)", y="Probability of disappearance") +
           scale_fill_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           scale_color_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           guides(fill=guide_legend(nrow=2)) +
           facet_grid(habitat.f ~ ForStratCat.f,
                      labeller = labeller(habitat.f=habitat.labs, ForStratCat.f=forstrat.labs)) +
           coord_cartesian(ylim = c(0,1)),  # Based on these graphs, fix N to ca. 10000
         width = 16.6, height=20, units="cm", dpi=350)
  
  ggsave("Diff_lost.pdf",
         ggplot(newdata_lost_diff, aes(x = diff, y = lost, color=DietCat.f)) +
           geom_path(aes(color=DietCat.f),size=0.5) +
           geom_ribbon(aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975, fill=DietCat.f), alpha = 0.2, linetype=0) +
           geom_point(data=df_lost_rev,
                      size=0.2, shape=4, color="gray60")  +
           theme(legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Difference in sampling effort", y="Probability of disappearance") +
           scale_fill_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           scale_color_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           guides(fill=guide_legend(nrow=2)) +
           facet_grid(habitat.f ~ ForStratCat.f,
                      labeller = labeller(habitat.f=habitat.labs, ForStratCat.f=forstrat.labs)) +
           coord_cartesian(ylim = c(0,0.5)),
         width = 16.6, height=20, units="cm", dpi=350)
  
  ### Based on this, it seems that the set sampling effort (and maybe difference in sampling effort) is not appropriate for making
  ### reasonable predictions. I would argue that fixing the change in sampling effort to anything but zero makes intuitive sense.
  ### However, the total number of records should be increased tremendously to get rid of the probability baseline.
  
  # Forage stratum
  set.seed(123)
  newdata_lost_forstrat <- {data.frame(longitude=569930.1, latitude=7031341, Developed = mean(df_bird_beta_rev$Developed), longevity_y=15, diff=0, N=10000, 
                                       DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                          levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                       ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                        "understory","wataroundsurf"), each=25*9), 4),
                                                              levels=c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                       "understory","wataroundsurf")),
                                       habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet",
                                                                        "open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                          levels=c("urban","generalist","marine","open","open_wet",
                                                                   "open_woodland","scrub","water","woodland")) ) }
  newdata_lost_forstrat$lost <- as.numeric(predict(spaMM_lost_final, newdata_lost_forstrat, re.form = NA, type="response")) 
  newdata_lost_forstrat <- cbind(newdata_lost_forstrat,
                                 get_intervals(spaMM_lost_final, newdata = newdata_lost_forstrat, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_lost_forstrat <- {newdata_lost_forstrat %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Plot:
  ggsave("Stratum_lost.pdf",
         ggplot(newdata_lost_forstrat, aes(x = habitat.f, y = lost, color=DietCat.f), alpha=0.5) +
           geom_point(size=1, position=position_dodge(width=0.5)) +
           geom_errorbar(aes(ymin=fixefVar_0.025, ymax=fixefVar_0.975), width=0.5, alpha=0.5,position=position_dodge(width=0.5)) +
           geom_point(data=df_lost_rev,
                      size=0.2, shape=4, alpha=0.5)  +
           theme(axis.text.x = element_text(angle = 90), legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Habitat", y="Probability of disappearance") +
           scale_fill_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           scale_color_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           guides(color=guide_legend(nrow=2)) +
           scale_x_discrete(labels=c("urban" = "Urban", "generalist" = "Generalist",  "marine" = "Marine", "open" = "Open", "woodland"="Woodland",
                                     "open_wet"="Wetland", "open_woodland"="Open \nwoodland","scrub"="Scrub", "water"="Water")) +
           facet_grid(. ~ ForStratCat.f,
                      labeller = labeller(ForStratCat.f=forstrat.labs, DietCat.f=diet.labs)) +
           coord_cartesian(ylim = c(0,0.5)),
         width = 16.6, height=12.5, units="cm", dpi=350)
  
  # Swap x-axis and facet:
  ggsave("Stratum_lost2.pdf",
         ggplot(newdata_lost_forstrat, aes(x = ForStratCat.f, y = lost, color=DietCat.f), alpha=0.5) +
           geom_point(size=1, position=position_dodge(width=0.5)) +
           geom_errorbar(aes(ymin=fixefVar_0.025, ymax=fixefVar_0.975), width=0.5, alpha=0.5,position=position_dodge(width=0.5)) +
           geom_point(data=df_lost_rev,
                      size=0.2, shape=4, alpha=0.5)  +
           theme(axis.text.x = element_text(angle = 90), legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Habitat", y="Probability of disappearance") +
           scale_fill_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           scale_color_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           guides(color=guide_legend(nrow=2)) +
           scale_x_discrete(labels=c("watbelowsurf" = "Water\n(below surf.)", "aerial" = "Aerial",  "canopy" = "Canopy",
                                     "ground" = "Ground", "midhigh"="Midhigh", "understory"="Understorey", "wataroundsurf"="Water\n(around surf.)")) +
           facet_grid(. ~ habitat.f,
                      labeller = labeller(habitat.f=habitat.labs, DietCat.f=diet.labs)) +
           coord_cartesian(ylim = c(0,0.5)),
         width = 16.6, height=12.5, units="cm", dpi=350)
}  

# Interaction effects: change in Developed area, habitat and diet
{
  # Developed area
  set.seed(123)
  newdata_lost_all_rev <- {data.frame(longitude=569930.1, latitude=7031341,
                                      Developed=c(rep(seq(from=min(df_lost_rev$Developed), to=max(df_lost_rev$Developed), length.out = 25), 4*7*9)),
                                      diff=0,
                                      longevity = 15,        
                                      N = 10000,
                                      DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                         levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                      ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf"), each=25*9), 4),
                                                             levels=c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf")),
                                      habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                         levels=c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland")) )
  }
  names(newdata_lost_all_rev) <- c("longitude","latitude","Developed","diff","longevity_y","N","DietCat.f","ForStratCat.f","habitat.f" )
  newdata_lost_all_rev$lost <- as.numeric(predict(spaMM_lost_final, newdata_lost_all_rev, re.form = NA, type="response"))
  newdata_lost_all_rev <- cbind(newdata_lost_all_rev,
                                get_intervals(spaMM_lost_final, newdata = newdata_lost_all_rev, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_lost_all_rev <- {newdata_lost_all_rev %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Plot and export:
  ## Labels for axes:
  {
    habitat.labs <- c("Urban","Generalist","Marine","Open","Wetland","Open\nwoodland","Scrub","Water","Woodland")
    names(habitat.labs) <- c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland")
    diet.labs <- c("Plant/seed","Invertebrate","Omnivore","Vertebrate/\nscavenging")
    names(diet.labs) <- c("PlantSeed","Invertebrate","Omnivore","VertFishScav")
    forstrat.labs <- c("Water\n(below surf.)","Aerial","Canopy","Ground","Midhigh","Understorey","Water\n(around surf.)")
    names(forstrat.labs) <- c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf")
  }
  
  ggsave("Response_lost_int.pdf",
         ggplot(data = newdata_lost_all_rev, aes(x = Developed, y = lost, color=ForStratCat.f, fill=ForStratCat.f)) +
           geom_point(data=df_lost_rev,
                      size=0.2, shape=4, color="gray60")  +
           geom_path(size=0.5) +   # , colour = "#00BFC4"
           geom_ribbon(aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975), alpha = 0.2, linetype=0) +  # , fill = "#00BFC4"
           labs(x="Change in Developed area [m2]", y="Probability of disappearance") +
           facet_grid(habitat.f ~ DietCat.f,
                      labeller = labeller(habitat.f=habitat.labs, DietCat.f=diet.labs)) +
           scale_fill_discrete(name = "Forage stratum", labels = c("Water (below surf.)","Aerial","Canopy","Ground",
                                                                   "Midhigh","Understorey","Water (around surf.)"))  +
           scale_color_discrete(name = "Forage stratum", labels = c("Water (below surf.)","Aerial","Canopy","Ground",
                                                                    "Midhigh","Understorey","Water (around surf.)"))  +
           guides(color=guide_legend(nrow=2)) +
           scale_x_continuous(expand = c(0, 0)) +
           theme(axis.text = element_text(size=5), legend.position = "bottom"),
         width = 16.6, height=20, units="cm", dpi=350)
}

##--- 7.2 Traits- species gain ---####

df_gained_rev <- droplevels(df_bird_trait.gained_rev[!(is.na(df_bird_trait.gained_rev$gained) | is.na(df_bird_trait.gained_rev$DietCat.f)),
                                                     c("Pixelnr","11","N","diff","longitude","latitude","longevity_y","adult_body_mass_g","habitat.f","DietCat.f","ForStratCat.f","gained")])
names(df_gained_rev) <- c("Pixelnr","Developed","N","diff","longitude","latitude","longevity_y","adult_body_mass_g","habitat.f","DietCat.f","ForStratCat.f","gained")

# Which Pixelnr in the trait dataframes are not in the beta-dataframes:
df_gained_rev[!(df_gained_rev$Pixelnr %in% df_bird_beta_rev$Pixelnr),]$Pixelnr  
# Remove the outliers:
df_gained_rev <- df_gained_rev[(df_gained_rev$Pixelnr %in% df_bird_beta_rev$Pixelnr),]

spaMM_gained_rev_upd <- step(fitme(gained ~ Developed + Developed*habitat.f +  Developed*DietCat.f +  Developed*ForStratCat.f +
                                     longevity_y + N + diff +
                                     (1 | Pixelnr) + Matern(1 | longitude+latitude),
                                   family = binomial, method = "ML", verbose=c(trace=T),
                                   init = list(rho=0.05,nu=1,lambda=0.5),  # To try and optimize speed 
                                   control = list(max.iter=50),     # Reduce the number of iterations due to computation time (defaults to 200)
                                   data=df_gained_rev))  

spaMM_gained_rev_upd2 <- fitme(gained ~ Developed + Developed*habitat.f +  DietCat.f +  Developed*ForStratCat.f +
                                 longevity_y + adult_body_mass_g + N + diff +
                                 (1 | Pixelnr) + Matern(1 | longitude+latitude),
                               family = binomial, method = "ML", verbose=c(trace=T),
                               init = list(rho=0.05,nu=1,lambda=0.5),  
                               control = list(max.iter=50),     # 
                               data=df_gained_rev) 

## Null model with only spatial effects
spaMM_gained_rev_null <- fitme(gained ~ 1 + (1 | Pixelnr) + Matern(1 | longitude+latitude),
                               family = binomial, method = "ML", verbose=c(trace=T),
                               init = list(rho=0.05,nu=1,lambda=0.5),  
                               control = list(max.iter=50),     
                               data=df_gained_rev)  
## Model with only spatial effects and sampling effort
spaMM_gained_rev_null2 <- fitme(gained ~ N + diff +
                                  (1 | Pixelnr) + Matern(1 | longitude+latitude),
                                family = binomial, method = "ML", verbose=c(trace=T),
                                init = list(rho=0.05,nu=1,lambda=0.5), 
                                control = list(max.iter=50),    
                                data=df_gained_rev) 

AIC(spaMM_gained_rev_upd)
AIC(spaMM_gained_rev_upd2)
AIC(spaMM_gained_rev_null)
AIC(spaMM_gained_rev_null2)

# Check spatial effects/correlation
plot(as.numeric(dist(df_gained_rev[,c("longitude","latitude")])),  # Distance matrix
     as.numeric(MaternCorr(dist(df_gained_rev[,c("longitude","latitude")]), nu = 0.640838187, rho = 0.003584977)),  # Matérn correlation, use nu and rho from the model summary
     xlab = "Distance between pairs of location [in m]", ylab = "Estimated correlation",
     xlim=c(0,5000)); abline(h=0.1, col="red", lty=2)   # ca. 1000 m before the correlaion goes below 0.1

plot(as.numeric(c(seq(from=0, to=5000 , by=100))),  # Distance matrix
     as.numeric(MaternCorr(c(seq(from=0, to=5000 , by=100)), nu = 0.640838187, rho = 0.003584977)),  # Matérn correlation, use nu and rho from the model summary
     xlab = "Distance between pairs of location [in m]", ylab = "Estimated correlation",
     xlim=c(0,5000)); abline(h=0.1, col="red", lty=2)   # ca. 1000 m before the correlaion goes below 0.1

# Check numbers more closely:
View(data.frame(dist = c(seq(from=720, to=730 , by=0.5)),
                corr = as.numeric(MaternCorr(c(seq(from=720, to=730 , by=0.5)), nu = 0.640838187, rho = 0.003584977))))

##--- 7.2.1 Model predictions ---####
# To be able to better visualise and interpret the models, create a hypothetical dataframe and make predictions incl.
# confidence intervals. I will exclude the spatial effects from the predictions, and thus only assess the effects of the other
# variables. When not assessed directly, I will leave each variable at the the baseline levels (factor), or the levels perviously
# specified as the mean (numerical). I will insert a constant longitude and latitude (a set of coordinates within Trondheim;
# median of all coordinates). As I am mainly interested in the differences in responses based on the interaction, I will need to
# let both variable fluctuate simultaneously

print(summary(spaMM_gained_rev_upd))  # To get full summary

# Linear effects: longevity, diet, N and diff
{
  # Sampling effort N
  set.seed(123)
  newdata_gained_N <- {data.frame(longitude=569930.1, latitude=7031341,
                                  Developed = mean(df_bird_beta_rev$Developed), diff=0, longevity = 15,        
                                  N = c(rep(seq(from=min(df_lost_rev$N), to=max(df_lost_rev$N), length.out = 25), 4*7*9)),
                                  DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                     levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                  ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf"), each=25*9), 4),
                                                         levels=c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf")),
                                  habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                     levels=c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland")) )
  }
  names(newdata_gained_N) <- c("longitude","latitude","Developed","diff","longevity_y","N","DietCat.f","ForStratCat.f","habitat.f" )
  newdata_gained_N$gained <- as.numeric(predict(spaMM_gained_rev_upd, newdata_gained_N, re.form = NA, type="response"))
  newdata_gained_N <- cbind(newdata_gained_N,
                            get_intervals(spaMM_gained_rev_upd, newdata = newdata_gained_N, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_gained_N <- {newdata_gained_N %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Difference in sampling effort, diff
  set.seed(123)
  newdata_gained_diff <- {data.frame(longitude=569930.1, latitude=7031341, Developed = mean(df_bird_beta_rev$Developed),longevity_y = 15, N=10000,
                                     diff = c(rep(seq(from=min(df_gained_rev$diff), to=max(df_gained_rev$diff), length.out = 25), 4*7*9)),
                                     DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                        levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                     ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                      "understory","wataroundsurf"), each=25*9), 4),
                                                            levels=c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                     "understory","wataroundsurf")),
                                     habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet",
                                                                      "open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                        levels=c("urban","generalist","marine","open","open_wet",
                                                                 "open_woodland","scrub","water","woodland")) )
  }
  newdata_gained_diff$gained <- as.numeric(predict(spaMM_gained_rev_upd, newdata_gained_diff, re.form = NA, type="response"))
  newdata_gained_diff <- cbind(newdata_gained_diff,
                               get_intervals(spaMM_gained_rev_upd, newdata = newdata_gained_diff, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_gained_diff <- {newdata_gained_diff %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Longevity
  set.seed(123)
  newdata_gained_longevity <- {data.frame(longitude=569930.1, latitude=7031341, Developed = mean(df_bird_beta_rev$Developed), diff=0, N=10000, 
                                          longevity_y = c(rep(seq(from=min(df_gained_rev$longevity_y),
                                                                  to=max(df_gained_rev$longevity_y), length.out = 25), 4*7*9)),
                                          DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                             levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                          ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                           "understory","wataroundsurf"), each=25*9), 4),
                                                                 levels=c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                          "understory","wataroundsurf")),
                                          habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet",
                                                                           "open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                             levels=c("urban","generalist","marine","open","open_wet",
                                                                      "open_woodland","scrub","water","woodland")) )
  }
  newdata_gained_longevity$gained <- as.numeric(predict(spaMM_gained_rev_upd, newdata_gained_longevity, re.form = NA, type="response"))
  newdata_gained_longevity <- cbind(newdata_gained_longevity,
                                    get_intervals(spaMM_gained_rev_upd, newdata = newdata_gained_longevity, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_gained_longevity <- { newdata_gained_longevity %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Plot:
  ggsave("Longevity_gained.pdf",
         ggplot(newdata_gained_longevity, aes(x = longevity_y, y = gained, color=DietCat.f)) +
           geom_path(aes(color=DietCat.f),size=0.5) +
           geom_ribbon(aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975, fill=DietCat.f), alpha = 0.2, linetype=0) +
           geom_point(data=df_gained_rev,
                      size=0.2, shape=4, color="gray60")  +
           theme(legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Longevity (years)", y="Probability of appearance") +
           scale_fill_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           scale_color_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           guides(fill=guide_legend(nrow=2)) +
           facet_grid(habitat.f ~ ForStratCat.f,
                      labeller = labeller(habitat.f=habitat.labs, ForStratCat.f=forstrat.labs)) +
           coord_cartesian(ylim = c(0,1)),   
         width = 16.6, height=20, units="cm", dpi=350)
  
  ggsave("N_gained.pdf",
         ggplot(newdata_gained_N, aes(x = N, y = gained, color=DietCat.f)) +
           geom_path(aes(color=DietCat.f),size=0.5) +
           geom_ribbon(aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975, fill=DietCat.f), alpha = 0.2, linetype=0) +
           geom_point(data=df_gained_rev,
                      size=0.2, shape=4, color="gray60")  +
           theme(legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Total number of observations (N)", y="Probability of appearance") +
           scale_fill_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           scale_color_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           guides(fill=guide_legend(nrow=2)) +
           facet_grid(habitat.f ~ ForStratCat.f,
                      labeller = labeller(habitat.f=habitat.labs, ForStratCat.f=forstrat.labs)) +
           coord_cartesian(ylim = c(0,1)),  
         width = 16.6, height=20, units="cm", dpi=350)
  
  ggsave("Diff_gained.pdf",
         ggplot(newdata_gained_diff, aes(x = diff, y = gained, color=DietCat.f)) +
           geom_path(aes(color=DietCat.f),size=0.5) +
           geom_ribbon(aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975, fill=DietCat.f), alpha = 0.2, linetype=0) +
           geom_point(data=df_gained_rev,
                      size=0.2, shape=4, color="gray60")  +
           theme(legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Difference in sampling effort", y="Probability of appearance") +
           scale_fill_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           scale_color_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/Scavenging"))  +
           guides(fill=guide_legend(nrow=2)) +
           facet_grid(habitat.f ~ ForStratCat.f,
                      labeller = labeller(habitat.f=habitat.labs, ForStratCat.f=forstrat.labs)) +
           coord_cartesian(ylim = c(0,1)),
         width = 16.6, height=20, units="cm", dpi=350)
  
  # Diet
  set.seed(123)
  newdata_gained_diet <- {data.frame(longitude=569930.1, latitude=7031341, Developed = mean(df_bird_beta_rev$Developed), longevity_y=15, diff=0, N=10000, 
                                     DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                        levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                     ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                      "understory","wataroundsurf"), each=25*9), 4),
                                                            levels=c("watbelowsurf","aerial","canopy","ground","midhigh",
                                                                     "understory","wataroundsurf")),
                                     habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet",
                                                                      "open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                        levels=c("urban","generalist","marine","open","open_wet",
                                                                 "open_woodland","scrub","water","woodland")) ) }
  newdata_gained_diet$gained <- as.numeric(predict(spaMM_gained_rev_upd, newdata_gained_diet, re.form = NA, type="response")) 
  newdata_gained_diet <- cbind(newdata_gained_diet,
                               get_intervals(spaMM_gained_rev_upd, newdata = newdata_gained_diet, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_gained_gained <- {newdata_gained_diet %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Plot:
  ggsave("Diet_gained.pdf",
         ggplot(newdata_gained_diet, aes(x = habitat.f, y = gained, color=ForStratCat.f), alpha=0.5) +
           geom_point(size=1, position=position_dodge(width=0.5)) +
           geom_errorbar(aes(ymin=fixefVar_0.025, ymax=fixefVar_0.975), width=0.5, alpha=0.5,position=position_dodge(width=0.5)) +
           geom_point(data=df_gained_rev,
                      size=0.2, shape=4, alpha=0.5)  +
           theme(axis.text.x = element_text(angle = 90), legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Habitat", y="Probability of appearance") +
           scale_fill_discrete(name = "Forage stratum", labels = c("Water /below surf.)","Aerial","Canopy","Ground","Midhigh",
                                                                   "Understorey","Water (around surf.)"))  +
           scale_color_discrete(name = "Forage stratum", labels = c("Water /below surf.)","Aerial","Canopy","Ground","Midhigh",
                                                                    "Understorey","Water (around surf.)"))  +
           guides(color=guide_legend(nrow=2)) +
           scale_x_discrete(labels=c("urban" = "Urban", "generalist" = "Generalist",  "marine" = "Marine", "open" = "Open", "woodland"="Woodland",
                                     "open_wet"="Wetland", "open_woodland"="Open \nwoodland","scrub"="Scrub", "water"="Water")) +
           facet_grid(. ~ DietCat.f,
                      labeller = labeller(ForStratCat.f=forstrat.labs, DietCat.f=diet.labs)) +
           coord_cartesian(ylim = c(0,0.6)),
         width = 16.6, height=12.5, units="cm", dpi=350)
  
  # Swap x-axis and facet:
  ggsave("Diet_gained2.pdf",
         ggplot(newdata_gained_diet, aes(x = DietCat.f, y = gained, color=ForStratCat.f), alpha=0.5) +
           geom_point(size=1, position=position_dodge(width=0.5)) +
           geom_errorbar(aes(ymin=fixefVar_0.025, ymax=fixefVar_0.975), width=0.5, alpha=0.5,position=position_dodge(width=0.5)) +
           geom_point(data=df_gained_rev,
                      size=0.2, shape=4, alpha=0.5)  +
           theme(axis.text.x = element_text(angle = 90), legend.position = "bottom", axis.text = element_text(size=5)) +
           labs(x="Main dietary component", y="Probability of appearance") +
           scale_fill_discrete(name = "Forage stratum", labels = c("Water /below surf.)","Aerial","Canopy","Ground","Midhigh",
                                                                   "Understorey","Water (around surf.)"))  +
           scale_color_discrete(name = "Forage stratum", labels = c("Water /below surf.)","Aerial","Canopy","Ground","Midhigh",
                                                                    "Understorey","Water (around surf.)"))  +
           guides(color=guide_legend(nrow=2)) +
           scale_x_discrete(labels=c("PlantSeed" = "Plant/seed", "Invertebrate" = "Invertebrate",
                                     "Omnivore" = "Omnivore", "VertFishScav" = "Vertebrate/\nscavenging")) +
           facet_grid(. ~ habitat.f,
                      labeller = labeller(ForStratCat.f=forstrat.labs, habitat.f=habitat.labs)) +
           coord_cartesian(ylim = c(0,0.6)),
         width = 16.6, height=12.5, units="cm", dpi=350)
}  

# Interaction effects: change in Developed area, habitat and forage stratum
{
  # Developed area
  set.seed(123)
  newdata_gained_all_rev <- {data.frame(longitude=569930.1, latitude=7031341,
                                        Developed=c(rep(seq(from=min(df_gained_rev$Developed), to=max(df_gained_rev$Developed), length.out = 25), 4*7*9)),
                                        diff=0,
                                        longevity = 15,        
                                        N = 10000,
                                        DietCat.f = factor(rep(c("PlantSeed","Invertebrate","Omnivore","VertFishScav"), each=25*7*9),
                                                           levels = c("PlantSeed","Invertebrate","Omnivore","VertFishScav")),
                                        ForStratCat.f = factor(rep(rep(c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf"), each=25*9), 4),
                                                               levels=c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf")),
                                        habitat.f = factor(rep(rep(rep(c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland"), each=25), 7), 4),
                                                           levels=c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland")) )
  }
  names(newdata_gained_all_rev) <- c("longitude","latitude","Developed","diff","longevity_y","N","DietCat.f","ForStratCat.f","habitat.f" )
  newdata_gained_all_rev$gained <- as.numeric(predict(spaMM_gained_rev_upd, newdata_gained_all_rev, re.form = NA, type="response"))
  newdata_gained_all_rev <- cbind(newdata_gained_all_rev,
                                  get_intervals(spaMM_gained_rev_upd, newdata = newdata_gained_all_rev, intervals = "fixefVar", re.form = NA, type="response") )
  # Remove predictions for combinations we do not have data for in the observed data:
  newdata_gained_all_rev <- {newdata_gained_all_rev %>%
      filter((ForStratCat.f=="watbelowsurf" & habitat.f=="marine" & (DietCat.f=="Invertebrate" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="watbelowsurf" & habitat.f=="water" & DietCat.f!="PlantSeed") |
               (ForStratCat.f=="aerial" & habitat.f=="generalist" & DietCat.f=="VertFishScav") |
               (ForStratCat.f=="aerial" & habitat.f=="open" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="aerial" & habitat.f=="open_woodland" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="canopy" & (habitat.f=="open" | habitat.f=="woodland") & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="ground" & habitat.f=="urban" & DietCat.f=="PlantSeed") |
               (ForStratCat.f=="ground" & (habitat.f=="generalist" | habitat.f=="marine") & (DietCat.f=="Omnivore" | DietCat.f=="VertFishScav")) |
               (ForStratCat.f=="ground" & (habitat.f=="open" | habitat.f=="open_woodland" | habitat.f=="woodland")) |
               (ForStratCat.f=="ground" & habitat.f=="open_wet" & DietCat.f!="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="scrub" & DietCat.f=="PlantSeed" ) |
               (ForStratCat.f=="ground" & habitat.f=="water" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="generalist" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="midhigh" & habitat.f=="open" & (DietCat.f=="PlantSeed" | DietCat.f=="Invertebrate")) |
               (ForStratCat.f=="midhigh" & habitat.f=="open_wet" & DietCat.f=="Invertebrate") |
               (ForStratCat.f=="midhigh" & habitat.f=="open_woodland" & DietCat.f=="Omnivore") |
               (ForStratCat.f=="midhigh" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="open" & DietCat.f=="Invertebrate" ) |
               (ForStratCat.f=="understory" & habitat.f=="open_woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="understory" & habitat.f=="woodland" & DietCat.f!="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="" & (DietCat.f=="" | DietCat.f=="")) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="generalist" & DietCat.f=="Omnivore" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="marine" & DietCat.f=="VertFishScav" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="water" ) |
               (ForStratCat.f=="wataroundsurf" & habitat.f=="woodland" & DietCat.f=="Invertebrate" ))}
  
  # Plot and export:
  ## Labels for axes:
  {
    habitat.labs <- c("Urban","Generalist","Marine","Open","Wetland","Open\nwoodland","Scrub","Water","Woodland")
    names(habitat.labs) <- c("urban","generalist","marine","open","open_wet","open_woodland","scrub","water","woodland")
    diet.labs <- c("Plant/seed","Invertebrate","Omnivore","Vertebrate/\nscavenging")
    names(diet.labs) <- c("PlantSeed","Invertebrate","Omnivore","VertFishScav")
    forstrat.labs <- c("Water\n(below surf.)","Aerial","Canopy","Ground","Midhigh","Understorey","Water\n(around surf.)")
    names(forstrat.labs) <- c("watbelowsurf","aerial","canopy","ground","midhigh","understory","wataroundsurf")
  }
  
  ggsave("Response_gained_int.pdf",
         ggplot(data = newdata_gained_all_rev, aes(x = Developed, y = gained, color=DietCat.f, fill=DietCat.f)) +
           geom_point(data=df_gained_rev,
                      size=0.2, shape=4, color="gray60")  +
           geom_path(size=0.5) +   # , colour = "#00BFC4"
           geom_ribbon(aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975), alpha = 0.2, linetype=0) +  # , fill = "#00BFC4"
           labs(x="Change in Developed area [m2]", y="Probability of appearance") +
           facet_grid(habitat.f ~ ForStratCat.f,
                      labeller = labeller(habitat.f=habitat.labs, ForStratCat.f=forstrat.labs)) +
           scale_fill_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/scavenging"))  +
           scale_color_discrete(name = "Main dietary component", labels = c("Plant/Seed", "Invertebrate", "Omnivore", "Vertebrate/scavenging"))  +
           guides(color=guide_legend(nrow=2)) +
           scale_x_continuous(expand = c(0, 0)) +
           theme(axis.text = element_text(size=5), legend.position = "bottom"),
         width = 16.6, height=20, units="cm", dpi=350)
}


##----------------------------####
##--- ADDITIONAL PLOTS ---####
# Change in developed area and beta_turnover combined in a single figure:
diff_plot <- diff_1218_rev
diff_plot$beta_logit <- ifelse(diff_plot$Pixelnr %in% df_bird_beta_rev$Pixelnr, diff_plot$beta_logit, NA)  # Subset to only pixels incl. in analyses

ras_pxl <- rasterize(as(diff_1218_rev, 'Spatial'), r500_rev, field="Pixelnr")  # Rasterize Pixelnr
ras_beta <- rasterize(as(diff_1218_rev, 'Spatial'), r500_rev, field="beta_logit")  # Rasterize beta
ras_dev <- ras_diff_1218_rev[["X_11"]]  # Subset to only "Developed"

ras_plot <- addLayer(addLayer(ras_pxl, ras_beta), ras_dev) # Combine with area change in a RasterStack
names(ras_plot) <- c("Pixelnr","beta_logit","Developed")

p_dev <- levelplot(ras_plot[["Developed"]],
                   main="",
                   margin=FALSE, xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                   par.settings = list(axis.line = list(col = "transparent"), layout.heights=list(xlab.key.padding=2)),
                   colorkey=list(space="bottom", width=0.75, height=0.75, labels=list(cex=0.5))) +
  layer(sp.text(loc = c(568180.1,7020250), txt="Change in Developed area (m2)", cex=0.6)) +
  layer(sp.text(loc = c(552500.1,7042591), txt="a)", cex=1))

p_beta <- levelplot(ras_plot[["beta_logit"]],
                    main="",
                    margin=FALSE, xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                    par.settings = list(axis.line = list(col = "transparent"), layout.heights=list(xlab.key.padding=2)),
                    colorkey=list(space="bottom", width=0.75, height=0.75, labels=list(cex=0.5))) +
  layer(sp.text(loc = c(568180.1,7020250), txt=expression("logit("*beta[turnover]*")"), cex=0.6)) +
  layer(sp.text(loc = c(552500.1,7042591), txt="b)", cex=1))

grid.arrange(Mydiverge0(p_dev, ramp='RdBu'),
             Mydiverge0(p_beta, ramp='YlOrRd'),
             ncol=2, nrow=1)

ggsave("maps.pdf",
       grid.arrange(Mydiverge0(p_dev, ramp='RdBu'),
                    Mydiverge0(p_beta, ramp='YlOrRd'),
                    ncol=2, nrow=1),
       width=16.6, height=10, units="cm", dpi=350)
