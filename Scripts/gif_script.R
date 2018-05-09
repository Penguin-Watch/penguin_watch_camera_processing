#####################
#Script to produce gifs from camera images
#
#Casey Youngflesh - Feb 5, 2018
#####################

###Note - set directory carefully (i.e. not on shared network).


#NOTES:
#requires .jpg images named with the following convention:
#SITEaYEARa_NUMBER.JPG [e.g., HALFc2013a_000001.JPG] - letter after year is not strictly necessary
#FUNCTION REQUIRES INSTALLATION OF IMAGEMAGICK on machine: https://www.imagemagick.org/script/index.php
#BE SURE TO CHECK 'INSTALL LEGACY UTILITIES' WHEN INSTALLING ON WINDOWS

#USE:
#navigate to directory with all camera images (should only contain camera images), then simply run function:
#>gif_fun()
#time delay between images can be specified with DELAY arg; size can be specified with SIZE arg

#FUNCTION:
#generates .gifs for all unique camera placements and puts them in a directory called '/gifs', inside the current directory with all images to be processed

rm(list=ls())
ls()

gif_fun <- function(DELAY = 10, SIZE = '384X316')
{
  #NOT ALL CAMS HAVE THE SAME NAMING CONVENTION - SOME HAVE LETTERS AFTER YEAR, SOME DO NOT
  
  #names of all images in folder
  images <- list.files()
  
  #test whether there is a letter after the year
  end_letter <- substr(images, 10, 10)
  lo <- function(x) !grepl("[^A-Za-z]", x)
  el_res <- which(lo(end_letter) == TRUE)
  
  #images with a letter after year
  img_letter <- images[el_res]
  gen_letter <- substr(img_letter, 1, 10) #just stem names
  
  #images with no letter after year
  img_no_letter <- images[-el_res]
  gen_no_letter <- substr(img_no_letter, 1, 9) #just stem names
  
  #combine all stem names
  stem <- c(gen_letter, gen_no_letter)
  #unique stem names
  u_stem <- unique(stem)
  
  #determine os
  if(Sys.info()[['sysname']] == 'Windows')
  {
    #resize and convert to gif using imagemagick
    system('cmd.exe', input = 'mkdir gifs')
    for (i in 1:length(u_stem))
    {
      system('cmd.exe ', input = paste0('convert -delay ', DELAY, 
                                        ' -resize ',  SIZE, ' ',
                                        u_stem[i], '*.JPG ',
                                        'gifs/',
                                        u_stem[i], '.gif'))
    }
  }
  if(Sys.info()[['sysname']] == 'Darwin')
  {
    #resize and convert to gif using imagemagick
    system('mkdir gifs')
    for (i in 1:length(u_stem))
    {
      system(paste0('convert -delay ', DELAY,
                    ' -resize ',  SIZE, ' ',
                    u_stem[i], '*.JPG ',
                    'gifs/',
                    u_stem[i], '.gif'))
    }
  }
}




# example use -------------------------------------------------------------

##setwd('~/Google_Drive/R/MISC/test_gif/')
gif_fun()
#gifs will be located in '/test_gif/gifs'
