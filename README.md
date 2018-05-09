# penguin_watch_camera_processing

Use `Consensus_data/`, `Nest_coords/`, and `Full_res_images/` with `data_plot_script.R` to create `Images_with_polys/`. Use those images to correct `Raw_data_files/` into `QC_data_files`.

## Scripts
* `/Scripts/data_plot_script.R` - plots nest polygons and chick 'consensus clicks' onto camera images
* `/Scripts/gif_script.R` - creates gif from camera images

## Shared folder structure
* `/Consensus_data/` - chick 'consensus clicks'
* `/Nest_coords/` - coordinates of nests for each site (may be more than one set per season)
* `/Full_res_images/` - full resolution jpg images from time-lapse cameras
* `/Low_res_images/` - low resultion jpg images from time-lapse cameras - **UNUSED**
* `/Images_with_polys/` - camera images that have nest polygons and chick consensus clicks plotted on top
* `/Raw_data_files/` - raw data files (number of chicks in each nest zone for each time step)
* `/QC_data_files/` - manually correct data files (number of chicks in each nest zone for each time step)

