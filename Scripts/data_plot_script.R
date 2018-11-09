######################
#Script to plot nest numbers, polygons, and clicks on camera images
#
#
#INSTRUCTIONS FOR POLYGON IMAGES:
#STEP 1 - set user dir
#STEP 2 - read in nest coordinate data
#STEP 3 - specify input and output dir
#- MAKE SURE R OBJECTS END WITH '/'
#- MAKE SURE DIRECTORIES EXIST (CREATE IF THEY DO NOT)
#STEP 5 - run function with `TYPE = 'POLY'` to create images with polygons
#
#
#INSTRUCTIONS FOR QC IMAGES (black and white images, with polygons, nest numbers, and clicks):
#STEP 1 - set user dir
#STEP 2 - read in nest coordinate data
#STEP 3 - read in click data
#STEP 3 - specify input and output dir
#- MAKE SURE R OBJECTS END WITH '/'
#- MAKE SURE DIRECTORIES EXIST (CREATE IF THEY DO NOT)
#STEP 5 - run function with `TYPE = 'BOTH'` to create QC images
#
#
#NOTES:
#-Imagemagick required
#-Scripts designed to run on macOS -> calls imagemagick through bash with `system()`
######################



# Setup -------------------------------------------------------------------

rm(list=ls())



# set user dirs ---------------------------------------------------------------

# Fiona
#dir <- ''

# Casey
dir <- '~/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/'

#dir <- 'C:/Users/Lynch Lab 7/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/'

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



# polygon function ----------------------------------------------------------------


#nest_coords is nest locations - UNTRANSFORMED Y (bottom left is origin) - always based on 1000 x 750
#consensus is consensus chick locations - UNTRANSFORMED Y (bottom left is origin) - always based on 1000 x 750
#jpeg_dir is jpeg directory
#output_dir is output directory
#dim is camera image dimesions
#poly_tr is transparency of polygon lines and red dots
#TYPE is whether to plot 1) JUST POLYGONS or 2) polygons, nest numbers, and consensus clicks ('POLY' or 'BOTH') - POLY does not require consensus clicks, just nest coords
#NEST_IMG_SZ is whether nest classifications were done on a full or scaled image - see google doc sheet for this information - PARITAL or FULL
#keep is whether to process all images or half the images - 'all' or 'half
#keep_oe is to keep odd or even images (with the first image being 1, not the actual image number) with keep = 'half' - SHOULD BE CHECKED TO MAKE SURE IT PROCESSES CORRECTLY - images processed earlier in time used 'even' (though this option wasn't yet available) - images processed later in time used 'odd' to accomodate keeping a quarter rather than a half - this option applies both to the polygon image plotting and both polygons and consensus clicks

pt_img_fun <- function(nest_coords, 
                       consensus, 
                       jpeg_dir, 
                       output_dir, 
                       dim = c(2048, 1536),
                       poly_tr = 0.2,
                       TYPE = 'POLY',
                       NEST_IMG_SZ = 'PARTIAL',
                       img_st = NULL,
                       img_end = NULL,
                       keep = 'all',
                       keep_oe = 'odd')
{
  #change permissions
  system(paste0('chmod -R 755 ', output_dir))
  
  # transform x and y coords
  trans_fun <- function(INPUT, TYPE = 'COORDS', NEST_IMG_SZ = 'PARTIAL', DIM = dim)
  {
    if (NEST_IMG_SZ == 'PARTIAL')
    {
      #scale factor
      if (DIM[1] == 2048)
      {
        x_scale <- 2048/1000
        y_scale <- 1536/750
      }
      if (DIM[1] == 1920)
      {
        x_scale <- 1920/1000
        y_scale <- 1080/562.5
      }
    
      #scale to original image size
      x_input_sc <- INPUT$x*x_scale
      y_input_sc <- INPUT$y*y_scale
    }
    if (NEST_IMG_SZ == 'FULL')
    {
      x_input_sc <- INPUT$x
      y_input_sc <- INPUT$y
    }
    if (NEST_IMG_SZ != 'PARTIAL' & NEST_IMG_SZ != 'FULL')
    {
      stop("Invalid arg for 'NEST_IMG_SZ'")
    }
    
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
  
  #ggplot colors function
  gg_color_hue <- function(n, ALPHA = 1)
  {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100, alpha = ALPHA)[1:n]
  }
  
  #-----------#
  # #test data
  # dir <- '~/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/'
  # jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2017')
  # output_dir <- paste0(dir, 'QC_images/test/')
  # nest_coords <- read.csv(paste0(dir, 'Nest_coords/GEORa2017_V1_nestcoords.csv'))
  # consensus <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/CUVEb2014_classifications.csv'))
  # jpeg_dir = jpeg_dir
  # output_dir = output_dir
  # #dim = c(1920, 1080)
  # dim = c(2048, 1536)
  # poly_tr = 0.6
  # TYPE = ''
  # img_st = 'GEORa2017a_000341'
  # img_end = 'GEORa2017a_000700'
  # plot_jpeg('~/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/Full_res_images/NEKOc2013/NEKOc2013c_000468.JPG')
  # NEST_COORDS <- trans_fun(nest_coords, TYPE = 'COORDS', NEST_IMG_SZ = 'FULL', DIM = dim)
  # CONSENSUS <- trans_fun(consensus, TYPE = 'CONSENSUS', NEST_IMG_SZ = 'FULL', DIM = dim)
  # #points(NEST_COORDS, col = rgb(1,0,0,0.8), pch = 19)
  # polys <- poly_fun(NEST_COORDS, DIM = dim)
  # for (j in 1:length(polys))
  # {
  #   #j <- 1
  #   #plot polygons
  #   lines(polys[[j]], lwd = 4, col = rgb(0,1,0,poly_tr))
  # }
  # points(CONSENSUS$x, CONSENSUS$y, col = rgb(0,1,0,0.8), pch = '.')
  # GEORa2017_nc_V1 <- read.csv(paste0(dir, 'Nest_coords/GEORa2017_V1_nestcoords.csv'))
  # GEORa2017_nc_V2 <- read.csv(paste0(dir, 'Nest_coords/GEORa2017_V2_nestcoords.csv'))
  # 
  # # set input/output
  # jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2017/')
  # output_dir <- paste0(dir, 'Images_with_polys/GEORa2017/')
  # 
  # # Run function - start end must be even numbers (if starting image is even) if keeping half
  # nest_coords = GEORa2017_nc_V1
  #            jpeg_dir = jpeg_dir
  #            output_dir = output_dir
  #            dim = c(2048, 1536)
  #            poly_tr = 0.6
  #            TYPE = 'POLY'
  #            NEST_IMG_SZ = 'FULL'
  #            img_st = 'GEORa2017a_000180.JPG'
  #            img_end = 'GEORa2017a_000340.JPG'
  #            keep = 'half'
  #-----------#
  
  NEST_COORDS <- trans_fun(nest_coords, NEST_IMG_SZ = NEST_IMG_SZ, TYPE = 'COORDS', DIM = dim)
  
  if (TYPE == 'BOTH')
  {
    CONSENSUS <- trans_fun(consensus, TYPE = 'CONSENSUS', NEST_IMG_SZ = 'FULL', DIM = dim)
  }
  
  #determine polygons from nest coordinates
  polys <- poly_fun(NEST_COORDS, DIM = dim)
  
  #get jpeg names
  jf <- list.files(path = jpeg_dir)
  jpeg_files_p <- jf[grep('.JPG', jf)]
  
  if (keep == 'half')
  {
    td <- seq(1, length(jpeg_files_p), 2)
    
    if (keep_oe == 'odd')
    {
      jpeg_files <- jpeg_files_p[td]
    }
    if (keep_oe == 'even')
    {
      jpeg_files <- jpeg_files_p[-td]
    }
  }
  if (keep == 'quarter')
  {
    td <- seq(1, length(jpeg_files_p), 4)
    jpeg_files <- jpeg_files_p[td]
  }
  if (keep == 'all')
  {
    jpeg_files <- jpeg_files_p
  }
  
  if (keep != 'half' & keep != 'quarter' & keep != 'all')
  {
    stop('valid args for keep are "half" and "all"')
  }
  
  #determine colors to use (as many as there are nests)
  cols <- gg_color_hue(NROW(NEST_COORDS), ALPHA = 0.7)
  
  #vector of numbers to plot for nests
  let <- str_pad(1:99, 2, pad = '0')
  
  if (!is.null(img_st) & !is.null(img_end))
  {
    ind_st <- grep(img_st, jpeg_files)
    ind_end <- grep(img_end, jpeg_files)
    jpeg_files_list <- jpeg_files[ind_st:ind_end]
  } else {
    jpeg_files_list <- jpeg_files
  }
    
  #loop to plot camera image with nest 'zones', nest numbers, and consensus chick clicks
  for (i in 1:length(jpeg_files_list))
  {
    #i <- 301
    #create jpeg
    jpeg(filename = paste0(output_dir, jpeg_files_list[i]), width = dim[1], height = dim[2])
    
    if (TYPE != 'BOTH')
    {
      #plot jpeg camera image
      plot_jpeg(paste0(jpeg_dir, jpeg_files_list[i]))
    } else {
      system(paste0('convert ',  jpeg_dir, jpeg_files_list[i], 
                    ' -set colorspace Gray -separate -average ', jpeg_dir, jpeg_files_list[i], 'tmp'))
      
      plot_jpeg(paste0(jpeg_dir, jpeg_files_list[i], 'tmp'))
      system(paste0('rm ', jpeg_dir, jpeg_files_list[i], 'tmp'))
    }
    
    #filter for consensus clicks from a single image
    jpg_name <- strsplit(jpeg_files_list[i], split = '.', fixed = TRUE)[[1]][1]
    
    if (TYPE == 'BOTH')
    {
      filt_clicks <- filter(CONSENSUS, name == jpg_name)
    }
    for (j in 1:length(polys))
    {
      #j <- 1
      if (TYPE != 'BOTH')
      {
      #plot polygons
      lines(polys[[j]], lwd = 3, col = rgb(1,0,0,poly_tr))
      }
      #nests numbers on plot
      if (TYPE == 'BOTH')
      {
        lines(polys[[j]], lwd = 3, col = rgb(0,0,1,poly_tr))
        text(NEST_COORDS$x[j], NEST_COORDS$y[j]-20, labels = let[j], 
             col = rgb(0,1,0, 0.8), #cols[j], 
             cex = 3)
      }
    }
    if (TYPE == 'BOTH')
    {
      #consensus clicks
      points(filt_clicks$x, filt_clicks$y, pch = 19, col = rgb(1,0,0,poly_tr), cex = 2.5)
    }
    dev.off()
  }
  #change permissions back
  system(paste0('chmod -R 555 ', output_dir))
}



# reduce polygon image size function ----------------------------------------------------------------

#determine which images are greater than 1MB -> reduce quality of these images to make file smaller (PW Pro will only take images under 1MB)
#jpeg_dir is directory that jpeg images are in


rd_img_fun <- function(jpeg_dir)
{
  #jpeg_dir <- '~/Desktop/test/'
  #change permissions
  system(paste0('chmod -R 755 ', jpeg_dir))
  
  #list files
  jf <- list.files(path = jpeg_dir)
  jpeg_files <- paste0(jpeg_dir, jf[grep('.JPG', jf)])
  #determine file sizes
  jpeg_sizes <- file.size(jpeg_files)
  #larger than 1MB
  large_files <- jpeg_files[which(jpeg_sizes >= 1000000)]
  
  if (length(large_files) > 0)
  {
    #quality level
    PER <- 85
    for (i in 1:length(large_files))
    {
      #determine number of characters in filename
      num_char <- nchar(large_files[i])
      #imagemagick to reduce file size
      system(paste0('convert -strip -interlace Plane -sampling-factor 4:2:0 -quality ', 
                    PER, '% ', large_files[i], ' ', substring(large_files[i], first = 1, last = num_char-4), '.JPG'))
    }
  }
  
  #recheck file sizes
  #determine file sizes
  jpeg_sizes <- file.size(jpeg_files)
  #larger than 1MB
  large_files <- jpeg_files[which(jpeg_sizes >= 1000000)]
  if (length(large_files) < 1)
  {
    print('SUCCESS!')
  } else {
    print('Looks like there are still some large files! Running through again.')
  
    PER <- 75
    for (i in 1:length(large_files))
    {
      #determine number of characters in filename
      num_char <- nchar(large_files[i])
      #imagemagick to reduce file size
      system(paste0('convert -strip -interlace Plane -sampling-factor 4:2:0 -quality ', PER, '% ', large_files[i], ' ', substring(large_files[i], first = 1, last = num_char-4), '.JPG'))
    }
  
    #recheck file sizes
    #determine file sizes
    jpeg_sizes <- file.size(jpeg_files)
    #larger than 1MB
    large_files <- jpeg_files[which(jpeg_sizes >= 1000000)]
    if (length(large_files) < 1)
    {
      print('SUCCESS!')
    } else {
      print('Looks like there are still some large files! Running through again.')
    }
  }
  
  #change permissions back
  system(paste0('chmod -R 555 ', jpeg_dir))
}




# AITCd2014 --------------------------------------------------------------

# #NEST COORDINATES
# AITCd2014_nc <- read.csv(paste0(dir, 'Nest_coords/AITCd2014a_nestcoords.csv'))
# #CONSENSUS CLICKS
# AITCd2014_con <- read.csv(paste0(dir, 'Consensus_data/PW_clicks/AITCd2014_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/AITCd2014/')
# output_dir <- paste0(dir, 'QC_images/AITCd2014/')
# 
# # Run function
# pt_img_fun(nest_coords = AITCd2014_nc,
#            consensus = AITCd2014_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            TYPE = 'BOTH',
#            poly_tr = 0.6,
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'all')





# BAILa2013 --------------------------------------------------------------


# #NEST COORDINATES
# BAILa2013_nc <- read.csv(paste0(dir, 'Nest_coords/BAILa2013_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BAILa2013/')
# output_dir <- paste0(dir, 'Images_with_polys/BAILa2013/')
# 
# # Run function
# pt_img_fun(nest_coords = BAILa2013_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'PARTIAL')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# BAILa2014 --------------------------------------------------------------

# #NEST COORDINATES
# BAILa2014_nc <- read.csv(paste0(dir, 'Nest_coords/BAILa2014_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BAILa2014/')
# output_dir <- paste0(dir, 'Images_with_polys/BAILa2014/')
# 
# # Run function
# pt_img_fun(nest_coords = BAILa2014_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'PARTIAL')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# BOOTb2013 --------------------------------------------------------------

# #NEST COORDINATES - three different sets
# BOOTb2013_nc_V1 <- read.csv(paste0(dir, 'Nest_coords/BOOTb2013_V1_nestcoords.csv'))
# BOOTb2013_nc_V2 <- read.csv(paste0(dir, 'Nest_coords/BOOTb2013_V2_nestcoords.csv'))
# BOOTb2013_nc_V3 <- read.csv(paste0(dir, 'Nest_coords/BOOTb2013_V3_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BOOTb2013/')
# output_dir <- paste0(dir, 'Images_with_polys/BOOTb2013/')
# 
# # Run function
# pt_img_fun(nest_coords = BOOTb2013_nc_V1,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            img_st = 'BOOTb2013a_000052',
#            img_end = 'BOOTb2014a_000454',
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'half')
# 
# pt_img_fun(nest_coords = BOOTb2013_nc_V2,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            img_st = 'BOOTb2014a_000455',
#            img_end = 'BOOTb2014a_000606',
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'half')
# 
# pt_img_fun(nest_coords = BOOTb2013_nc_V3,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            img_st = 'BOOTb2014a_000607',
#            img_end = 'BOOTb2014a_001080',
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# GEORa2013 --------------------------------------------------------------

# #NEST COORDINATES
# GEORa2013_nc <- read.csv(paste0(dir, 'Nest_coords/GEORa2013_nestcoords.csv'))
# #CONSENSUS CLICKS
# GEORa2013_con <- read.csv(paste0(dir, 'Consensus_data/PW_clicks/GEORa2013_consensus.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2013/')
# output_dir <- paste0(dir, 'QC_images/GEORa2013/')
# 
# # Run function
# pt_img_fun(nest_coords = GEORa2013_nc,
#            consensus = GEORa2013_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.3,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'PARTIAL',
#            img_st = 'GEORa2013b_000001.JPG',
#            img_end = 'GEORa2013b_000508.JPG',
#            keep = 'all')




# GEORa2014 --------------------------------------------------------------

# #NEST COORDINATES
# GEORa2014_nc <- read.csv(paste0(dir, 'Nest_coords/GEORa2014_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2014/')
# output_dir <- paste0(dir, 'Images_with_polys/GEORa2014/')
# 
# # Run function
# pt_img_fun(nest_coords = GEORa2014_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'all')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# GEORa2014_nc <- read.csv(paste0(dir, 'Nest_coords/GEORa2014_nestcoords.csv'))
# 
# #PW PRO CLICKS
# GEORa2014_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/GEORa2014_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2014/')
# output_dir <- paste0(dir, 'QC_images/GEORa2014/')
# 
# # Run function
# pt_img_fun(nest_coords = GEORa2014_nc,
#            consensus = GEORa2014_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'all')



# HALFb2013 --------------------------------------------------------------

# #NEST COORDINATES
# HALFb2013_nc <- read.csv(paste0(dir, 'Nest_coords/HALFb2013_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/HALFb2013/')
# output_dir <- paste0(dir, 'Images_with_polys/HALFb2013/')
# 
# # Run function
# pt_img_fun(nest_coords = HALFb2013_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'all')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# LOCKb2013 --------------------------------------------------------------

# #NEST COORDINATES
# LOCKb2013_nc_V1 <- read.csv(paste0(dir, 'Nest_coords/LOCKb2013_V1_nestcoords.csv'))
# LOCKb2013_nc_V2 <- read.csv(paste0(dir, 'Nest_coords/LOCKb2013_V2_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/LOCKb2013/')
# output_dir <- paste0(dir, 'Images_with_polys/LOCKb2013/')
# 
# #Run function
# pt_img_fun(nest_coords = LOCKb2013_nc_V1,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'PARTIAL',
#            img_st = 'LOCKb2013b_000241',
#            img_end = 'LOCKb2013b_000388',
#            keep = 'all')
# 
# pt_img_fun(nest_coords = LOCKb2013_nc_V2,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'PARTIAL',
#            img_st = 'LOCKb2013b_000389',
#            img_end = 'LOCKb2013b_000770',
#            keep = 'all')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# #BOTH
# #NEST COORDINATES
# LOCKb2013_nc_V1 <- read.csv(paste0(dir, 'Nest_coords/LOCKb2013_V1_nestcoords.csv'))
# LOCKb2013_nc_V2 <- read.csv(paste0(dir, 'Nest_coords/LOCKb2013_V2_nestcoords.csv'))
# 
# #PW PRO CLICKS
# LOCKb2013_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/LOCKb2013_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/LOCKb2013/')
# output_dir <- paste0(dir, 'QC_images/LOCKb2013/')
# 
# # Run function
# pt_img_fun(nest_coords = LOCKb2013_nc_V1,
#            consensus = LOCKb2013_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'PARTIAL',
#            img_st = 'LOCKb2013b_000241',
#            img_end = 'LOCKb2013b_000388',
#            keep = 'all')
# 
# pt_img_fun(nest_coords = LOCKb2013_nc_V2,
#            consensus = LOCKb2013_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'PARTIAL',
#            img_st = 'LOCKb2013b_000389',
#            img_end = 'LOCKb2013b_000770',
#            keep = 'all')







# LOCKb2014 --------------------------------------------------------------

# #NEST COORDINATES
# LOCKb2014_nc <- read.csv(paste0(dir, 'Nest_coords/LOCKb2014_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/LOCKb2014/')
# output_dir <- paste0(dir, 'Images_with_polys/LOCKb2014/')
# 
# # Run function
# pt_img_fun(nest_coords = LOCKb2014_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# #BOTH
# #NEST COORDINATES
# LOCKb2014_nc <- read.csv(paste0(dir, 'Nest_coords/LOCKb2014_nestcoords.csv'))
# 
# #PW PRO CLICKS
# LOCKb2014_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/LOCKb2014_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/LOCKb2014/')
# output_dir <- paste0(dir, 'QC_images/LOCKb2014/')
# 
# # Run function
# pt_img_fun(nest_coords = LOCKb2014_nc,
#            consensus = LOCKb2014_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'all')




# NEKOc2013 --------------------------------------------------------------

# #NEST COORDINATES
# NEKOc2013_nc <- read.csv(paste0(dir, 'Nest_coords/NEKOc2013_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/NEKOc2013/')
# output_dir <- paste0(dir, 'Images_with_polys/NEKOc2013/')
# 
# #Run function
# pt_img_fun(nest_coords = NEKOc2013_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            keep = 'half',
#            NEST_IMG_SZ = 'PARTIAL')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# NEKOc2013_nc <- read.csv(paste0(dir, 'Nest_coords/NEKOc2013_nestcoords.csv'))
# 
# #PW PRO CLICKS
# NEKOc2013_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/NEKOc2013_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/NEKOc2013/')
# output_dir <- paste0(dir, 'QC_images/NEKOc2013/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = NEKOc2013_nc,
#            consensus = NEKOc2013_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'PARTIAL',
#            keep = 'half',
#            keep_oe = 'even')



# ORNEa2014 --------------------------------------------------------------

# #NEST COORDINATES
# ORNEa2014_nc_V1 <- read.csv(paste0(dir, 'Nest_coords/ORNEa2014_V1_nestcoords.csv'))
# ORNEa2014_nc_V2 <- read.csv(paste0(dir, 'Nest_coords/ORNEa2014_V2_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/ORNEa2014/')
# output_dir <- paste0(dir, 'Images_with_polys/ORNEa2014/')
# 
# # Run function
# pt_img_fun(nest_coords = ORNEa2014_nc_V1,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            img_st = 'SPIGa2014a_000207',
#            img_end = 'SPIGa2014b_000163',
#            keep = 'all',
#            NEST_IMG_SZ = 'PARTIAL')
# 
# 
# pt_img_fun(nest_coords = ORNEa2014_nc_V2,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            img_st = 'ORNEa2014c_000001',
#            img_end = 'ORNEa2014c_000450',
#            keep = 'all',
#            NEST_IMG_SZ = 'PARTIAL')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# GEORa2015 ---------------------------------------------------------------


# #NEST COORDINATES
# GEORa2015_nc_V1 <- read.csv(paste0(dir, 'Nest_coords/GEORa2015_V1_nestcoords.csv'))
# GEORa2015_nc_V2 <- read.csv(paste0(dir, 'Nest_coords/GEORa2015_V2_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2015/')
# output_dir <- paste0(dir, 'Images_with_polys/GEORa2015/')
# 
# # Run function
# pt_img_fun(nest_coords = GEORa2015_nc_V1,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'GEORa2015a_000001',
#            img_end = 'GEORa2015a_000624',
#            keep = 'half')
# 
# 
# pt_img_fun(nest_coords = GEORa2015_nc_V2,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'GEORa2015a_000625',
#            img_end = 'GEORa2015a_001219',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# GEORa2015_nc_V1 <- read.csv(paste0(dir, 'Nest_coords/GEORa2015_V1_nestcoords.csv'))
# GEORa2015_nc_V2 <- read.csv(paste0(dir, 'Nest_coords/GEORa2015_V2_nestcoords.csv'))
# 
# #PW PRO CLICKS
# GEORa2015_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/GEORa2015_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2015/')
# output_dir <- paste0(dir, 'QC_images/GEORa2015/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = GEORa2015_nc_V1,
#            consensus = GEORa2015_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'GEORa2015a_000001',
#            img_end = 'GEORa2015a_000624',
#            keep = 'all')
# 
# pt_img_fun(nest_coords = GEORa2015_nc_V2,
#            consensus = GEORa2015_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'GEORa2015a_000625',
#            img_end = 'GEORa2015a_001219',
#            keep = 'all')



# CUVEb2014 ---------------------------------------------------------------


# #NEST COORDINATES
# CUVEb2014_nc <- read.csv(paste0(dir, 'Nest_coords/CUVEb2014a_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/CUVEb2014/')
# output_dir <- paste0(dir, 'Images_with_polys/CUVEb2014/')
# 
# # Run function
# pt_img_fun(nest_coords = CUVEb2014_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'all')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# CUVEb2014_nc <- read.csv(paste0(dir, 'Nest_coords/CUVEb2014a_nestcoords.csv'))
# 
# #PW PRO CLICKS
# CUVEb2014_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/CUVEb2014_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/CUVEb2014/')
# output_dir <- paste0(dir, 'QC_images/CUVEb2014/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = CUVEb2014_nc,
#            consensus = CUVEb2014_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'all')



# DAMOa2014 ---------------------------------------------------------------


# #NEST COORDINATES
# DAMOa2014_nc <- read.csv(paste0(dir, 'Nest_coords/DAMOa2014_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/DAMOa2014/')
# output_dir <- paste0(dir, 'Images_with_polys/DAMOa2014/')
# 
# # Run function
# pt_img_fun(nest_coords = DAMOa2014_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# DAMOa2014_nc <- read.csv(paste0(dir, 'Nest_coords/DAMOa2014_nestcoords.csv'))
# 
# #PW PRO CLICKS
# DAMOa2014_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/DAMOa2014_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/DAMOa2014/')
# output_dir <- paste0(dir, 'QC_images/DAMOa2014/')
# 
# # Run function
# pt_img_fun(nest_coords = DAMOa2014_nc,
#            consensus = DAMOa2014_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'all')



# DANCb2013 ---------------------------------------------------------------


# #NEST COORDINATES
# DANCb2013_nc <- read.csv(paste0(dir, 'Nest_coords/DANCb2013_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/DANCb2013/')
# output_dir <- paste0(dir, 'Images_with_polys/DANCb2013/')
# 
# # Run function
# pt_img_fun(nest_coords = DANCb2013_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# DANCb2013_nc <- read.csv(paste0(dir, 'Nest_coords/DANCb2013_nestcoords.csv'))
# 
# #PW PRO CLICKS
# DANCb2013_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/DANCb2013_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/DANCb2013/')
# output_dir <- paste0(dir, 'QC_images/DANCb2013/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = DANCb2013_nc,
#            consensus = DANCb2013_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half',
#            keep_oe = 'even')



# DANCb2014 ---------------------------------------------------------------


# #NEST COORDINATES
# DANCb2014_nc <- read.csv(paste0(dir, 'Nest_coords/DANCb2014_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/DANCb2014/')
# output_dir <- paste0(dir, 'Images_with_polys/DANCb2014/')
# 
# # Run function
# pt_img_fun(nest_coords = DANCb2014_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# DANCb2014_nc <- read.csv(paste0(dir, 'Nest_coords/DANCb2014_nestcoords.csv'))
# 
# #PW PRO CLICKS
# DANCb2014_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/DANCb2014_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/DANCb2014/')
# output_dir <- paste0(dir, 'QC_images/DANCb2014/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = DANCb2014_nc,
#            consensus = DANCb2014_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half',
#            keep_oe = 'even')



# CUVEb2015 ---------------------------------------------------------------

#deleted everything before CUVEb2015a_000372 as there was significant cam movement

# #NEST COORDINATES
# CUVEb2015_nc <- read.csv(paste0(dir, 'Nest_coords/CUVEb2015_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/CUVEb2015/')
# output_dir <- paste0(dir, 'Images_with_polys/CUVEb2015/')
# 
# # Run function
# pt_img_fun(nest_coords = CUVEb2015_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# CUVEb2015_nc <- read.csv(paste0(dir, 'Nest_coords/CUVEb2015_nestcoords.csv'))
# 
# #PW PRO CLICKS
# CUVEb2015_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/CUVEb2015_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/CUVEb2015/')
# output_dir <- paste0(dir, 'QC_images/CUVEb2015/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = CUVEb2015_nc,
#            consensus = CUVEb2015_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half',
#            keep_oe = 'even')



# MAIVb2013 ---------------------------------------------------------------


# #NEST COORDINATES
# MAIVb2013_nc <- read.csv(paste0(dir, 'Nest_coords/MAIVb2013a_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/MAIVb2013/')
# output_dir <- paste0(dir, 'Images_with_polys/MAIVb2013/')
# 
# # Run function
# pt_img_fun(nest_coords = MAIVb2013_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# MAIVb2013_nc <- read.csv(paste0(dir, 'Nest_coords/MAIVb2013a_nestcoords.csv'))
# 
# #PW PRO CLICKS
# MAIVb2013_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/MAIVb2013_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/MAIVb2013/')
# output_dir <- paste0(dir, 'QC_images/MAIVb2013/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = MAIVb2013_nc,
#            consensus = MAIVb2013_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half',
#            keep_oe = 'even')


# MAIVc2013 ---------------------------------------------------------------


# #NEST COORDINATES
# MAIVc2013_nc <- read.csv(paste0(dir, 'Nest_coords/MAIVc2013_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/MAIVc2013/')
# output_dir <- paste0(dir, 'Images_with_polys/MAIVc2013/')
# 
# # Run function
# pt_img_fun(nest_coords = MAIVc2013_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


#BOTH
#NEST COORDINATES
MAIVc2013_nc <- read.csv(paste0(dir, 'Nest_coords/MAIVc2013_nestcoords.csv'))

#PW PRO CLICKS
MAIVc2013_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/MAIVc2013_classifications.csv'))

# set input/output
jpeg_dir <- paste0(dir, 'Full_res_images/MAIVc2013/')
output_dir <- paste0(dir, 'QC_images/MAIVc2013/')


# Run function
pt_img_fun(nest_coords = MAIVc2013_nc,
           consensus = MAIVc2013_con,
           jpeg_dir = jpeg_dir,
           output_dir = output_dir,
           dim = c(1920, 1080),
           poly_tr = 0.6,
           TYPE = 'BOTH',
           NEST_IMG_SZ = 'FULL',
           keep = 'half',
           keep_oe = 'even')


# GEORa2017 ---------------------------------------------------------------

# #NEST COORDINATES
# GEORa2017_nc_V1 <- read.csv(paste0(dir, 'Nest_coords/GEORa2017_V1_nestcoords.csv'))
# GEORa2017_nc_V2 <- read.csv(paste0(dir, 'Nest_coords/GEORa2017_V2_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2017/')
# output_dir <- paste0(dir, 'Images_with_polys/GEORa2017/')
# 
# # Run function - start end must be odd numbers (if starting image is even) if keeping half
# pt_img_fun(nest_coords = GEORa2017_nc_V1,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'GEORa2017a_000181.JPG',
#            img_end = 'GEORa2017a_000341.JPG',
#            keep = 'half')
# 
# pt_img_fun(nest_coords = GEORa2017_nc_V2,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'GEORa2017a_000343',
#            img_end = 'GEORa2017a_000699',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# #BOTH
# #NEST COORDINATES
# GEORa2017_nc_V1 <- read.csv(paste0(dir, 'Nest_coords/GEORa2017_V1_nestcoords.csv'))
# GEORa2017_nc_V2 <- read.csv(paste0(dir, 'Nest_coords/GEORa2017_V2_nestcoords.csv'))
# 
# #PW PRO CLICKS
# GEORa2017_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/GEORa2017_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/GEORa2017/')
# output_dir <- paste0(dir, 'QC_images/GEORa2017/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = GEORa2017_nc_V1,
#            consensus = GEORa2017_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'GEORa2017a_000181.JPG',
#            img_end = 'GEORa2017a_000341.JPG',
#            keep = 'half',
#            keep_oe = 'even')
# 
# pt_img_fun(nest_coords = GEORa2017_nc_V2,
#            consensus = GEORa2017_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'GEORa2017a_000343',
#            img_end = 'GEORa2017a_000699',
#            keep = 'half',
#            keep_oe = 'even')


# LOCKb2015 ---------------------------------------------------------------


# #NEST COORDINATES
# LOCKb2015_nc <- read.csv(paste0(dir, 'Nest_coords/LOCKb2015a_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/LOCKb2015/')
# output_dir <- paste0(dir, 'Images_with_polys/LOCKb2015/')
# 
# # Run function
# pt_img_fun(nest_coords = LOCKb2015_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'quarter')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# #BOTH
# #NEST COORDINATES
# LOCKb2015_nc <- read.csv(paste0(dir, 'Nest_coords/LOCKb2015a_nestcoords.csv'))
# 
# #PW PRO CLICKS
# LOCKb2015_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/LOCKb2015a_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/LOCKb2015/')
# output_dir <- paste0(dir, 'QC_images/LOCKb2015/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = LOCKb2015_nc,
#            consensus = LOCKb2015_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'quarter',
#            keep_oe = 'even')




# AITCb2015 ---------------------------------------------------------------


# #NEST COORDINATES
# AITCb2015_nc <- read.csv(paste0(dir, 'Nest_coords/AITCb2015_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/AITCb2015/')
# output_dir <- paste0(dir, 'Images_with_polys/AITCb2015/')
# 
# # Run function
# pt_img_fun(nest_coords = AITCb2015_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# DANCb2015 ---------------------------------------------------------------


# #NEST COORDINATES
# DANCb2015_v1_nc <- read.csv(paste0(dir, 'Nest_coords/DANCb2015_V1_nestcoords.csv'))
# DANCb2015_v2_nc <- read.csv(paste0(dir, 'Nest_coords/DANCb2015_V2_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/DANCb2015/')
# output_dir <- paste0(dir, 'Images_with_polys/DANCb2015/')
# 
# # Run function
# pt_img_fun(nest_coords = DANCb2015_v1_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'DANCb2015a_000001.JPG',
#            img_end = 'DANCb2015a_000551.JPG',
#            keep = 'half')
# 
# pt_img_fun(nest_coords = DANCb2015_v2_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'DANCb2015b_000001.JPG',
#            img_end = 'DANCb2015b_000209.JPG',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# NEKOc2014 ---------------------------------------------------------------


# #NEST COORDINATES
# NEKOc2014_nc <- read.csv(paste0(dir, 'Nest_coords/NEKOc2014_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/NEKOc2014/')
# output_dir <- paste0(dir, 'Images_with_polys/NEKOc2014/')
# 
# # Run function
# pt_img_fun(nest_coords = NEKOc2014_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# NEKOc2015 ---------------------------------------------------------------


# #NEST COORDINATES
# NEKOc2015_nc <- read.csv(paste0(dir, 'Nest_coords/NEKOc2015_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/NEKOc2015/')
# output_dir <- paste0(dir, 'Images_with_polys/NEKOc2015/')
# 
# # Run function
# pt_img_fun(nest_coords = NEKOc2015_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# NEKOc2016 ---------------------------------------------------------------


# #NEST COORDINATES
# NEKOc2016_nc <- read.csv(paste0(dir, 'Nest_coords/NEKOc2016_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/NEKOc2016/')
# output_dir <- paste0(dir, 'Images_with_polys/NEKOc2016/')
# 
# # Run function
# pt_img_fun(nest_coords = NEKOc2016_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# BEAVa2016 ---------------------------------------------------------------


# #NEST COORDINATES
# BEAVa2016_nc <- read.csv(paste0(dir, 'Nest_coords/BEAVa2016_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BEAVa2016/')
# output_dir <- paste0(dir, 'Images_with_polys/BEAVa2016/')
# 
# # Run function
# pt_img_fun(nest_coords = BEAVa2016_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# COOPa2015 ---------------------------------------------------------------

# #NEST COORDINATES
# COOPa2015_nc <- read.csv(paste0(dir, 'Nest_coords/COOPa2015_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/COOPa2015/')
# output_dir <- paste0(dir, 'Images_with_polys/COOPa2015/')
# 
# # Run function
# pt_img_fun(nest_coords = COOPa2015_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)

# #BOTH
# #NEST COORDINATES
# COOPa2015_nc <- read.csv(paste0(dir, 'Nest_coords/COOPa2015_nestcoords.csv'))
# 
# #PW PRO CLICKS
# COOPa2015_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/COOPa2015_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/COOPa2015/')
# output_dir <- paste0(dir, 'QC_images/COOPa2015/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = COOPa2015_nc,
#            consensus = COOPa2015_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half',
#            keep_oe = 'odd')



# COOPx2017 ---------------------------------------------------------------

# #NEST COORDINATES
# COOPx2017_nc <- read.csv(paste0(dir, 'Nest_coords/COOPx2017a_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/COOPx2017/')
# output_dir <- paste0(dir, 'Images_with_polys/COOPx2017/')
# 
# # Run function
# pt_img_fun(nest_coords = COOPx2017_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# #BOTH
# #NEST COORDINATES
# COOPx2017_nc <- read.csv(paste0(dir, 'Nest_coords/COOPx2017a_nestcoords.csv'))
# 
# #PW PRO CLICKS
# COOPx2017_con <- read.csv(paste0(dir, 'Consensus_data/PW_Pro_clicks/COOPx2017_classifications.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/COOPx2017/')
# output_dir <- paste0(dir, 'QC_images/COOPx2017/')
# 
# 
# # Run function
# pt_img_fun(nest_coords = COOPx2017_nc,
#            consensus = COOPx2017_con,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'BOTH',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half',
#            keep_oe = 'odd')





# OCEAa2015 ---------------------------------------------------------------

# #NEST COORDINATES
# OCEAa2015_nc <- read.csv(paste0(dir, 'Nest_coords/OCEAa2015_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/OCEAa2015/')
# output_dir <- paste0(dir, 'Images_with_polys/OCEAa2015/')
# 
# # Run function
# pt_img_fun(nest_coords = OCEAa2015_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.9,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)



# BEAVa2015 ---------------------------------------------------------------

# #NEST COORDINATES
# BEAVa2015_nc <- read.csv(paste0(dir, 'Nest_coords/BEAVa2015_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BEAVa2015/')
# output_dir <- paste0(dir, 'Images_with_polys/BEAVa2015/')
# 
# # Run function
# pt_img_fun(nest_coords = BEAVa2015_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.9,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# BEAVa2016 ---------------------------------------------------------------
# 
# #NEST COORDINATES
# BEAVa2016_nc <- read.csv(paste0(dir, 'Nest_coords/BEAVa2016a_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BEAVa2016/')
# output_dir <- paste0(dir, 'Images_with_polys/BEAVa2016/')
# 
# # Run function
# pt_img_fun(nest_coords = BEAVa2016_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.9,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# BEAVa2017 ---------------------------------------------------------------

# #NEST COORDINATES
# BEAVa2017_nc <- read.csv(paste0(dir, 'Nest_coords/BEAVa2017a_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BEAVa2017/')
# output_dir <- paste0(dir, 'Images_with_polys/BEAVa2017/')
# 
# # Run function
# pt_img_fun(nest_coords = BEAVa2017_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.9,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# BLUFb2016 ---------------------------------------------------------------

# #NEST COORDINATES
# BLUFb2016_nc <- read.csv(paste0(dir, 'Nest_coords/BLUFb2016_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/BLUFb2016/')
# output_dir <- paste0(dir, 'Images_with_polys/BLUFb2016/')
# 
# # Run function
# pt_img_fun(nest_coords = BLUFb2016_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(1920, 1080),
#            poly_tr = 0.9,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# # GODHb2018 ---------------------------------------------------------------
# 
# #NEST COORDINATES
# GODHb2018_nc <- read.csv(paste0(dir, 'Nest_coords/GODHb2018_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/GODHb2018/')
# output_dir <- paste0(dir, 'Images_with_polys/GODHb2018/')
# 
# # Run function
# pt_img_fun(nest_coords = GODHb2018_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.9,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)


# # JOUGa2017 ---------------------------------------------------------------
# 
# #NEST COORDINATES
# JOUGa2017_v1_nc <- read.csv(paste0(dir, 'Nest_coords/JOUGa2017_V1_nestcoords.csv'))
# JOUGa2017_v2_nc <- read.csv(paste0(dir, 'Nest_coords/JOUGa2017_V2_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/JOUGa2017/')
# output_dir <- paste0(dir, 'Images_with_polys/JOUGa2017/')
# 
# #run function
# pt_img_fun(nest_coords = JOUGa2017_v1_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'JOUGa2017a_000001.JPG',
#            img_end = 'JOUGa2017a_000197.JPG',
#            keep = 'half')
# 
# pt_img_fun(nest_coords = JOUGa2017_v2_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.6,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            img_st = 'JOUGa2017a_000199.JPG',
#            img_end = 'JOUGa2017a_000799.JPG',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)




# # MIKKa2017 ---------------------------------------------------------------
# 
# #NEST COORDINATES
# MIKKa2017_nc <- read.csv(paste0(dir, 'Nest_coords/MIKKa2017a_nestcoords.csv'))
# 
# # set input/output
# jpeg_dir <- paste0(dir, 'Full_res_images/MIKKa2017/')
# output_dir <- paste0(dir, 'Images_with_polys/MIKKa2017/')
# 
# # Run function
# pt_img_fun(nest_coords = MIKKa2017_nc,
#            jpeg_dir = jpeg_dir,
#            output_dir = output_dir,
#            dim = c(2048, 1536),
#            poly_tr = 0.9,
#            TYPE = 'POLY',
#            NEST_IMG_SZ = 'FULL',
#            keep = 'half')
# 
# #reduce size of large images for PW Pro
# rd_img_fun(output_dir)




