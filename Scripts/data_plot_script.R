######################
#Script to plot nest numbers, teselation, and chicks consensus on one plot
#
#
######################


# Setup -------------------------------------------------------------------

rm(list=ls())


# Load packages -----------------------------------------------------------

#install.packages("deldir", dependencies = TRUE)
#install.packages("SDMTools", dependencies = TRUE)
#install.packages("dplyr", dependencies = TRUE)
##install.packages("tidyr", dependencies = TRUE)
##install.packages("jpeg", dependencies = TRUE)

require(dplyr)
library(jpeg)
library(stringr)



# function ----------------------------------------------------------------


#nest_coords is nest locations with y transformed for upper left origin (not lower left)
#consensus is consensus chick locations with y transformed for upper left origin (not lower left)
#jpeg_dir is jpeg directory


pt_img_fun <- function(nest_coords, consensus, jpeg_dir, output_dir)
{
  poly_fun <- function(KM_REV_ORTHO)
  {
    #KM_REV_ORTHO <- km_rev_ortho
    width <- 1000
    height <- 750
    
    #Voronoi tesselation using specified nest sites - deldir makes tesselation
    vt <- suppressWarnings(deldir(KM_REV_ORTHO[,1], KM_REV_ORTHO[,2], 
                                  rw= c(0, width, 0, height)))
    
    w <- tile.list(vt) #polygon coordinates
    
    polys <- vector(mode= 'list', length= length(w))
    for (i in seq(along= polys))
    {
      pcrds <- cbind(w[[i]]$x, w[[i]]$y)
      pcrds <- rbind(pcrds, pcrds[1,])
      polys[[i]] <- pcrds #arrange polygon coordinates
    }
    return(polys)
  }
  
  #determine polygons from nest coordinates
  polys <- poly_fun(nest_coords)
  
  #get jpeg names
  jpeg_files <- list.files(path = jpeg_dir)
  
  #ggplot colors function
  gg_color_hue <- function(n, OUT = 'HEX')
  {
    #n = 6
    hues = seq(15, 375, length=n+1)
    tmp <- hcl(h=hues, l=65, c=100)[1:n]
    gg_cols <- col2rgb(tmp)/255
    
    if (OUT == 'HEX')
    {
      return(tmp)
    }
    if (OUT == 'RGB')
    {
      return(gg_cols)
    }
    if (OUT != 'HEX' & OUT != 'RGB')
    {
      stop('"OUT" argument must be "HEX" or "RGB"')
    }
  }
  
  #plot jpeg function
  plot_jpeg = function(path, add=FALSE)
  {
    par(mar=  c(0,0,0,0))
    jpg = readJPEG(path, native=T) # read the file
    res = dim(jpg)[1:2] # get the resolution
    if (!add) # initialize an empty plot area if add==FALSE
      plot(1,1,xlim=c(1,res[2]),ylim=c(1,res[1]),asp=1,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
    rasterImage(jpg,1,1,res[2],res[1])
    par(mar= c(5, 4, 4, 2))
  }
  
  #determine colors to use (as many as there are nests)
  cols <- gg_color_hue(NROW(nest_coords))
  
  #vector of numbers to plot for nests
  let <- str_pad(1:99, 2, pad = '0')
  
  #loop to plot camera image with nest 'zones', nest numbers, and consensus chick clicks
  dir <- substr(jpeg_files[1], start = 1, stop = 10)
  dir.create(paste0(output_dir, '/', dir))
  for (i in 1:length(jpeg_files))
  {
    #i <- 1
    #create jpeg
    jpeg(filename = paste0(output_dir, '/', dir, '/', jpeg_files[i]), width = 1000, height = 750)
    #plot jpeg camera image
    plot_jpeg(paste0(jpeg_dir, '/', jpeg_files[i]))
    
    #filter for consensus clicks from a single image
    jpg_name <- strsplit(jpeg_files[i], split = '.', fixed = TRUE)[[1]][1]
    
    filt_clicks <- filter(consensus, name == jpg_name)
    for (j in 1:length(polys))
    {
      #plot polygons
      lines(polys[[j]], lwd = 5)
      #nests numbers on plot
      text(nest_coords$x[j], nest_coords$y[j], labels = let[j], col = cols[j], lwd = 12)
    }
    #consensus clicks
    points(filt_clicks$x, filt_clicks$y, pch = 19, col = 'red', lwd = 4)
    dev.off()
  }
}




# Load data ---------------------------------------------------------------

#data <- read.csv("C:/Users/lady3793/Dropbox/Penguin_DPhil/Survival_paper/LOCKb2014/LOCKb2014_nestcoords.csv", header = TRUE, sep = ",")
#data_user <- read.csv("C:/Users/lady3793/Dropbox/Penguin_DPhil/Survival_paper/Filtered_clusters/LOCKb_filtered2.csv", header = TRUE, sep = ",")

dir <- '~/Google_Drive/R/Project_archive/Old_penguin_watch_camera_processing/Scripts/'
#NEST COORDINATES
nest_coords_p <- read.csv(paste0(dir, 'LOCKb2014_nestcoords.csv'))
#CONSENSUS CLICKS
consensus_p <- read.csv(paste0(dir, 'LOCKb2014_data_user.csv'))
#IMAGE DIR
path <- paste0(dir, 'LOCKb2014b_jpeg/')



# transform y coordinates -------------------------------------------------


#transform y coordinates as image origin is top left, rather than bottom left
nest_coords <- data.frame(x = nest_coords_p$x, 
                          y = 750 - nest_coords_p$y)
consensus <- data.frame(name = consensus_p$name, 
                        x = consensus_p$x, 
                        y = 750 - consensus_p$y)


# run function ------------------------------------------------------------

jpeg_dir <- "/Users/caseyyoungflesh/Google_Drive/R/Project_archive/Old_penguin_watch_camera_processing/Scripts/LOCKb2014b_jpeg"
output_dir <- "/Users/caseyyoungflesh/Google_Drive/R/Project_archive/Old_penguin_watch_camera_processing/Scripts/Output_jpeg/"



pt_img_fun(nest_coords, consensus, jpeg_dir, output_dir)
