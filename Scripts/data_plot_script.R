######################
#Script to plot nest numbers, teselation, and chicks consensus on one plot
#
#INSTRUCTIONS:
#STEP 1 - run code to L140
#STEP 2 - set dir
#STEP 3 - copy code from previous site run and change object names as needed
#STEP 4 - create output dir in `images_with_polys/` (specified with `output_dir` object in R)
#NOTE: make sure dir objects in R end with '/'
#STEP 5 - run code to specify objects
#STEP 6 - run function to create images
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

library(SDMTools)
library(deldir)
library(dplyr)
library(jpeg)
library(stringr)



# function ----------------------------------------------------------------


#nest_coords is nest locations - UNTRANSFORMED Y (bottom left is origin)
#consensus is consensus chick locations - UNTRANSFORMED Y (bottom left is origin)
#jpeg_dir is jpeg directory
#output_dir is output directory


pt_img_fun <- function(nest_coords, 
                       consensus, 
                       jpeg_dir, 
                       output_dir, 
                       dim = c(2048, 1536))
{
  # transform y coords
  trans_fun <- function(INPUT, TYPE = 'COORDS', DIM = dim)
  {
    #scale factor
    x_scale <- 2048/1000
    y_scale <- 1536/750
    
    #scale to original image size
    x_input_sc <- INPUT$x*x_scale
    y_input_sc <- INPUT$y*y_scale
    
    if (TYPE == 'COORDS')
    {
      OUTPUT <- data.frame(x = x_input_sc, y = DIM[2] - y_input_sc)
    }
    if (TYPE == 'CONSENSUS')
    {
      OUTPUT <- data.frame(name = INPUT$name, x = x_input_sc, y = DIM[2] - y_input_sc)
    }
    return(OUTPUT)
  }
  
  #polygon function
  poly_fun <- function(KM_REV_ORTHO, DIM = dim)
  {
    width <- DIM[1]
    height <- DIM[2]
    
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
  
  #-----------#
  #test data
  # dir <- '~/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/'
  # AITCd2014_nc <- read.csv(paste0(dir, 'Nest_coords/AITCd2014a_nestcoords.csv'))
  # AITCd2014_con <- read.csv(paste0(dir, 'Consensus_data/AITCd2014_consensus.csv'))
  # dim = c(2048, 1536)
  # nest_coords <- trans_fun(AITCd2014_nc, TYPE = 'COORDS', DIM = dim)
  # consensus <- trans_fun(AITCd2014_con, TYPE = 'CONSENSUS', DIM = dim)
  # jpeg_dir <- paste0(dir, 'Full_res_images/AITCd2014/')
  # output_dir <- paste0(dir, 'Images_with_polys/AITCd2014/')
  #-----------#
  
  #determine polygons from nest coordinates
  polys <- poly_fun(nest_coords, DIM = dim)
  
  #get jpeg names
  jf <- list.files(path = jpeg_dir)
  jpeg_files <- jf[grep('.JPG', jf)]

  #ggplot colors function
  gg_color_hue <- function(n, ALPHA = 1)
  {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100, alpha = ALPHA)[1:n]
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
  cols <- gg_color_hue(NROW(nest_coords), ALPHA = 0.7)
  
  #vector of numbers to plot for nests
  let <- str_pad(1:99, 2, pad = '0')
  
  #loop to plot camera image with nest 'zones', nest numbers, and consensus chick clicks
  for (i in 1:length(jpeg_files))
  {
    #i <- 301
    #create jpeg
    jpeg(filename = paste0(output_dir, jpeg_files[i]), width = dim[1], height = dim[2])
    #plot jpeg camera image
    plot_jpeg(paste0(jpeg_dir, jpeg_files[i]))
    
    #filter for consensus clicks from a single image
    jpg_name <- strsplit(jpeg_files[i], split = '.', fixed = TRUE)[[1]][1]
    
    filt_clicks <- filter(consensus, name == jpg_name)
    for (j in 1:length(polys))
    {
      #j <- 1
      #plot polygons
      lines(polys[[j]], lwd = 3, col = rgb(1,0,0,0.15))
      #nests numbers on plot
      #if(i < 5)
      #{
        text(nest_coords$x[j], nest_coords$y[j]-10, labels = let[j], col = cols[j], cex = 1.5)
      #}
    }
    #consensus clicks
    points(filt_clicks$x, filt_clicks$y, pch = 19, col = rgb(1,0,0,0.2), lwd = 3)
    dev.off()
  }
}







# set user dirs ---------------------------------------------------------------

# Fiona
#NEST COORDINATES
#nest_coords_p <- read.csv("C:/Users/lady3793/Dropbox/Penguin_DPhil/Survival_paper/LOCKb2014/LOCKb2014_nestcoords.csv", header = TRUE, sep = ",")
#CONSENSUS CLICKS
#data_user <- read.csv("C:/Users/lady3793/Dropbox/Penguin_DPhil/Survival_paper/Filtered_clusters/LOCKb_filtered2.csv", header = TRUE, sep = ",")

# Casey

dir <- '~/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/'





# AITCd2014 --------------------------------------------------------------

#NEST COORDINATES
AITCd2014_nc <- read.csv(paste0(dir, 'Nest_coords/AITCd2014a_nestcoords.csv'))
#CONSENSUS CLICKS
AITCd2014_con <- read.csv(paste0(dir, 'Consensus_data/AITCd2014_consensus.csv'))

# set input/output
jpeg_dir <- paste0(dir, 'Full_res_images/AITCd2014/')
output_dir <- paste0(dir, 'Images_with_polys/AITCd2014/')

# Run function
pt_img_fun(AITCd2014_nc, AITCd2014_con, jpeg_dir, output_dir, dim = c(2048, 1536))




# AITCd2014 --------------------------------------------------------------

#NEST COORDINATES
AITCd2014_nc <- read.csv(paste0(dir, 'Nest_coords/AITCd2014a_nestcoords.csv'))
#CONSENSUS CLICKS
AITCd2014_con <- read.csv(paste0(dir, 'Consensus_data/AITCd2014_consensus.csv'))

# set input/output
jpeg_dir <- paste0(dir, 'Full_res_images/AITCd2014/')
output_dir <- paste0(dir, 'Images_with_polys/AITCd2014/')

# Run function
pt_img_fun(AITCd2014_nc, AITCd2014_con, jpeg_dir, output_dir, dim = c(2048, 1536))
