# bacteria_segmentation

## Description

This script allows for batch segmentation of bacteria in 3D images using [omnipose](https://www.nature.com/articles/s41592-022-01639-4), filter out some by volume or position in the image (touching borders) and then returns the mean and median intensity of the bacteria in predefined channels.

## Requirements

The script is made in Python and requires [Fiji](https://doi.org/10.1038/nmeth.2019) to be ran. On top of this, multiple update sites need to be activated following [this guide](https://imagej.net/update-sites/#following-an-update-site): 
* BIOP
* MorpholibJ

On top of that, a python environment containing a modified version of [omnipose](https://github.com/kevinjohncutler/omnipose) to handle specific `cudatoolkit` version corresponding to the hardware available in-house. To recreate the environment, a YAML file is available which will put all the packages to the required versions.

Once the environment is created and the update sites activated, just drag and drop the script in the main Fiji window and click on the RUN button.3

As Fiji is operating system independant, this should be runnable on Windows, Mac and Linux. As Deep Learning is involved, a machine with a Nvidia graphics card is suggested, to make the process faster. 

## Run script

### Input

The script will prompt for a folder containing the images to be analysed as well as a file extension to only analyze the targeted files. A calibrated size for the bacteria is required as well as the omnipose environment path to run the python command and find the bacteria.

### Runtime

For all selected files, the images will be opened using [Bio-Formats](https://doi.org/10.1083/jcb.201004104), a median filter is applied on all slices to smooth and help the segmentation and the images are then segmented using omnipose and the predefined model `bact_fluor_omni` and a label image is outputed.

Using the [3D ROI Manager](http://dx.doi.org/10.1093/bioinformatics/btt276), measurements are done on the individual bacteria to remove objects that are below a volume threshold or objects touching the X and Y borders of the image. Mean and median intensities are then measured on the remaining bacteria for both GFP and mCherry and saved in a CSV. 

### Output

A CSV containing all the intensity measurements is saved and the remaining bacteria are exported as 3D ROIs allowing the user to double check the segmentation and correct some wrongly detected bacteria.
