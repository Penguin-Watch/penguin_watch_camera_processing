######################
#Script to plot nest numbers, tessellation, and chicks consensus clicks on camera images
#
#Used for data QC
#
#INSTRUCTIONS:
#STEP 1 - set user dir
#STEP 2 - read in consensus data and nest coordinates
#STEP 3 - specify input and output dir
    #- MAKE SURE R OBJECTS END WITH '/'
    #- MAKE SURE DIRECTORIES EXIST (CREATE IF THEY DO NOT)
#STEP 5 - run function to create images with polygons
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


#nest_coords is nest locations - UNTRANSFORMED Y (bottom left is origin) - always based on 1000 x 750
#consensus is consensus chick locations - UNTRANSFORMED Y (bottom left is origin) - always based on 1000 x 750
#jpeg_dir is jpeg directory
#output_dir is output directory
#dim is camera image dimesions
#poly_tr is transparency of polygon lines and red dots


pt_img_fun <- function(nest_coords, 
                       consensus, 
                       jpeg_dir, 
                       output_dir, 
                       dim = c(2048, 1536),
                       poly_tr = 0.2)
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
  # poly_tr = 0.3
  # dir <- '~/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/'
  # nest_coords <- read.csv(paste0(dir, 'Nest_coords/BAILa2013_nestcoords.csv'))
  # consensus <- read.csv(paste0(dir, 'Consensus_data/BAILa2013_consensus.csv'))
  # consensus$name <- paste0('BAILa2013', substring(BAILa2013_con$name, 10, 17))
  # dim = c(2048, 1536)
  # jpeg_dir <- paste0(dir, 'Full_res_images/BAILa2013/')
  # output_dir <- paste0(dir, 'Images_with_polys/BAILa2013/')
  #-----------#
   
  NEST_COORDS <- trans_fun(nest_coords, TYPE = 'COORDS', DIM = dim)
  CONSENSUS <- trans_fun(consensus, TYPE = 'CONSENSUS', DIM = dim)
  
  #determine polygons from nest coordinates
  polys <- poly_fun(NEST_COORDS, DIM = dim)
  
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
  cols <- gg_color_hue(NROW(NEST_COORDS), ALPHA = 0.7)
  
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
    
    filt_clicks <- filter(CONSENSUS, name == jpg_name)
    for (j in 1:length(polys))
    {
      #j <- 1
      #plot polygons
      lines(polys[[j]], lwd = 3, col = rgb(1,0,0,poly_tr))
      #nests numbers on plot
      #if(i < 5)
      #{
        text(NEST_COORDS$x[j], NEST_COORDS$y[j]-20, labels = let[j], 
             col = rgb(1,0,0, 0.8), #cols[j], 
             cex = 1.5)
      #}
    }
    #consensus clicks
    points(filt_clicks$x, filt_clicks$y, pch = 19, col = rgb(1,0,0,poly_tr), lwd = 3)
    dev.off()
  }
}







# set user dirs ---------------------------------------------------------------

# Fiona
#dir <- ''

# Casey
dir <- '~/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/'




# AITCd2014 --------------------------------------------------------------

# #NEST COORDINATES
# AITCd2014_nc <- read.csv(paste0(dir, 'Nest_coords/AITCd2014a_nestcoords.csv'))
# #CONSENSUS CLICKS
# AITCd2014_con <- read.csv(paste0(dir, 'Consensus_data/AITCd2014_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/AITCd2014/')
# output_dir <- paste0(dir, 'Images_with_polys/AITCd2014/')
# 
# # Run function
# pt_img_fun(AITCd2014_nc, AITCd2014_con, jpeg_dir, output_dir, dim = c(2048, 1536))




# BAILa2013 --------------------------------------------------------------


# #NEST COORDINATES
# BAILa2013_nc <- read.csv(paste0(dir, 'Nest_coords/BAILa2013_nestcoords.csv'))
# #CONSENSUS CLICKS
# BAILa2013_con <- read.csv(paste0(dir, 'Consensus_data/BAILa2013_consensus.csv'))
# #images names were changed because wrong year was there - changed conensus names bc of this
# BAILa2013_con$name <- paste0('BAILa2013', substring(BAILa2013_con$name, 10, 17))
# 
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BAILa2013/')
# output_dir <- paste0(dir, 'Images_with_polys/BAILa2013/')
# 
# # Run function
# pt_img_fun(BAILa2013_nc, BAILa2013_con, jpeg_dir, output_dir, 
#            dim = c(2048, 1536), poly_tr = 0.3)




# BAILa2014 --------------------------------------------------------------

# #NEST COORDINATES
# BAILa2014_nc <- read.csv(paste0(dir, 'Nest_coords/BAILa2014_nestcoords.csv'))
# #CONSENSUS CLICKS
# BAILa2014_con <- read.csv(paste0(dir, 'Consensus_data/BAILa2014_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BAILa2014/')
# output_dir <- paste0(dir, 'Images_with_polys/BAILa2014/')
# 
# # Run function
# pt_img_fun(BAILa2014_nc, BAILa2014_con, jpeg_dir, output_dir, dim = c(2048, 1536))




# BOOTb2013 --------------------------------------------------------------

# #NEST COORDINATES
# BOOTb2013_nc <- read.csv(paste0(dir, 'Nest_coords/BOOTb2013_nestcoords.csv'))
# #CONSENSUS CLICKS
# BOOTb2013_con <- read.csv(paste0(dir, 'Consensus_data/BOOTb2013_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BOOTb2013/')
# output_dir <- paste0(dir, 'Images_with_polys/BOOTb2013/')
# 
# # Run function
# pt_img_fun(BOOTb2013_nc, BOOTb2013_con, jpeg_dir, output_dir, dim = c(2048, 1536))




# GEORa2013 --------------------------------------------------------------

#NEST COORDINATES
GEORa2013_nc <- read.csv(paste0(dir, 'Nest_coords/GEORa2013_nestcoords.csv'))
#CONSENSUS CLICKS
GEORa2013_con <- read.csv(paste0(dir, 'Consensus_data/GEORa2013_consensus.csv'))

# set input/output
jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2013/')
output_dir <- paste0(dir, 'Images_with_polys/GEORa2013/')

# Run function
pt_img_fun(GEORa2013_nc, GEORa2013_con, jpeg_dir, output_dir, 
           dim = c(2048, 1536), poly_tr = 0.3)




# GEORa2014 --------------------------------------------------------------

# #NEST COORDINATES
# GEORa2014_nc <- read.csv(paste0(dir, 'Nest_coords/GEORa2014_nestcoords.csv'))
# #CONSENSUS CLICKS
# GEORa2014_con <- read.csv(paste0(dir, 'Consensus_data/GEORa2014_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2014/')
# output_dir <- paste0(dir, 'Images_with_polys/GEORa2014/')
# 
# # Run function
# pt_img_fun(GEORa2014_nc, GEORa2014_con, jpeg_dir, output_dir, dim = c(2048, 1536))




# HALFb2013 --------------------------------------------------------------

# #NEST COORDINATES
# HALFb2013_nc <- read.csv(paste0(dir, 'Nest_coords/HALFb2013_nestcoords.csv'))
# #CONSENSUS CLICKS
# HALFb2013_con <- read.csv(paste0(dir, 'Consensus_data/HALFb2013_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/HALFb2013/')
# output_dir <- paste0(dir, 'Images_with_polys/HALFb2013/')
# 
# # Run function
# pt_img_fun(HALFb2013_nc, HALFb2013_con, jpeg_dir, output_dir, dim = c(2048, 1536))




# LOCKb2013 --------------------------------------------------------------

# #NEST COORDINATES
# LOCKb2013_nc <- read.csv(paste0(dir, 'Nest_coords/LOCKb2013_nestcoords.csv'))
# #CONSENSUS CLICKS
# LOCKb2013_con <- read.csv(paste0(dir, 'Consensus_data/LOCKb2013_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/LOCKb2013/')
# output_dir <- paste0(dir, 'Images_with_polys/LOCKb2013/')
# 
# # Run function
# pt_img_fun(LOCKb2013_nc, LOCKb2013_con, jpeg_dir, output_dir, dim = c(2048, 1536))




# LOCKb2014 --------------------------------------------------------------

# #NEST COORDINATES
# LOCKb2014_nc <- read.csv(paste0(dir, 'Nest_coords/LOCKb2014_nestcoords.csv'))
# #CONSENSUS CLICKS
# LOCKb2014_con <- read.csv(paste0(dir, 'Consensus_data/LOCKb2014_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/LOCKb2014/')
# output_dir <- paste0(dir, 'Images_with_polys/LOCKb2014/')
# 
# # Run function
# pt_img_fun(LOCKb2014_nc, LOCKb2014_con, jpeg_dir, output_dir, dim = c(2048, 1536))




# NEKOc2013 --------------------------------------------------------------

# #NEST COORDINATES
# NEKOc2013_nc <- read.csv(paste0(dir, 'Nest_coords/NEKOc2013_nestcoords.csv'))
# #CONSENSUS CLICKS
# NEKOc2013_con <- read.csv(paste0(dir, 'Consensus_data/NEKOc2013_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/NEKOc2013/')
# output_dir <- paste0(dir, 'Images_with_polys/NEKOc2013/')
# 
# # Run function
# pt_img_fun(NEKOc2013_nc, NEKOc2013_con, jpeg_dir, output_dir, dim = c(2048, 1536))




# ORNEa2014 --------------------------------------------------------------
# 
# #NEST COORDINATES
# ORNEa2014_nc <- read.csv(paste0(dir, 'Nest_coords/ORNEa2014_nestcoords.csv'))
# #CONSENSUS CLICKS
# ORNEa2014_con <- read.csv(paste0(dir, 'Consensus_data/ORNEa2014_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/ORNEa2014/')
# output_dir <- paste0(dir, 'Images_with_polys/ORNEa2014/')
# 
# # Run function
# pt_img_fun(ORNEa2014_nc, ORNEa2014_con, jpeg_dir, output_dir, dim = c(2048, 1536))


