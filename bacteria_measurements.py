# ─── SCRIPT PARAMETERS ──────────────────────────────────────────────────────────

#@ File(label="Folder with your images", style="directory", description="Folder with the images") src_dir
#@ String(label="Extension for the images to look for", value="lif") filename_filter
#@ Double(label="Bacteria diameter", value=1) bact_diameter_calibrated
#@ File(label="Path to omnipose environment", style="directory", description="Path to omnipose env") omnipose_env
#@ CommandService command

# ─── IMPORTS ────────────────────────────────────────────────────────────────────

import os
import io
import sys
import csv
import re
import glob
from itertools import izip
from operator import itemgetter

from ij import IJ, Prefs
from ij.plugin import Duplicator, ZProjector, ImagesToStack, RGBStackMerge, StackReverser, ChannelArranger, ImageCalculator
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.gui import PointRoi
from ij.measure import ResultsTable

# Bioformats imports
from loci.plugins import BF
from loci.plugins import LociExporter
from loci.plugins.in import ImporterOptions
from loci.plugins.out import Exporter
from loci.formats.in import MetadataOptions
from loci.formats import ImageReader
from loci.formats import MetadataTools

from ch.epfl.biop.ij2command import Labels2Rois
from ch.epfl.biop.wrappers.cellpose.ij2commands import Cellpose_SegmentImgPlusAdvanced, CellposePrefsSet


# 3DSuite imports
from mcib3d.geom import Objects3DPopulation
from mcib3d.image3d import ImageInt, ImageHandler

# ─── FUNCTIONS ──────────────────────────────────────────────────────────────────

def sorted_alphanumeric(data):
    """Sort a list alphanumerically

    Parameters
    ----------
    data : list
        List containing all the files to sort

    Returns
    -------
    list
        List with filenames sorted
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(data, key=alphanum_key)

def list_all_filenames(source, filetype):
    """Get a sorted list of all files of specified filetype in a given directory

    Parameters
    ----------
    source : str
        Path to source dir
    filetype : str
        File extension to specify filetype

    Returns
    -------
    List
        List of all files of the given type in the source dir
    """


    # os.chdir(str(source))
    return sorted_alphanumeric(glob.glob(os.path.join(source, "*"+filetype)))

def open_single_series_with_BF(path_to_file, series_number):
    """Open a single serie for a file using Bio-Formats

    Parameters
    ----------
    path_to_file : str
        Path to the file
    series_number : int
        Number of the serie to open

    Returns
    -------
    ImagePlus
        ImagePlus of the serie
    """
    options = ImporterOptions()
    options.setColorMode(ImporterOptions.COLOR_MODE_COMPOSITE)
    options.setSeriesOn(series_number, True) # python starts at 0
    # options.setSpecifyRanges(True)
    # options.setCBegin(series_number-1, channel_number-1) # python starts at 0
    # options.setCEnd(series_number-1, channel_number-1)
    # options.setCStep(series_number-1, 1)
    options.setId(path_to_file)
    imps = BF.openImagePlus(options) # is an array of imp with one entry

    return imps[0]



def BFExport(implus, savepath):
    """Export using BioFormats

    Parameters
    ----------
    implus : ImagePlus
        ImagePlus of the file to save
    savepath : str
        Path where to save the image

    """
    paramstring = "outfile=[" + savepath + "] windowless=true compression=Uncompressed saveROI=false"


    print('Savepath: ', savepath)
    plugin     = LociExporter()
    plugin.arg = paramstring
    exporter   = Exporter(plugin, implus)
    exporter.run()

def get_series_count_from_ome_metadata(path_to_file):
    """Get the number of series from a file

    Parameters
    ----------
    path_to_file : str
        Path to the file

    Returns
    -------
    int
        Number of series for the file
    """
    reader = ImageReader()
    omeMeta = MetadataTools.createOMEXMLMetadata()
    reader.setMetadataStore(omeMeta)
    reader.setId(path_to_file)
    series_count = reader.getSeriesCount()
    reader.close()

    return series_count


def progressbar(progress, total, line_number, prefix=''):
    """Progress bar for the IJ log window

    Parameters
    ----------
    progress : int
        Current step of the loop
    total : int
        Total number of steps for the loop
    line_number : int
        Number of the line to be updated
    prefix : str, optional
        Text to use before the progress bar, by default ''
    """

    size = 30
    x    = int(size*progress/total)
    IJ.log("\\Update%i:%s[%s%s] %i/%i\r" % (line_number, prefix, "#"*x, "."*(size-x), progress, total))

def get_file_info(file):
    """Get information about a file

    Parameters
    ----------
    file : str
        Path to a file

    Returns
    -------
    str
        Folder and basename without extension of the file
    """

    folder   = os.path.dirname(file)
    basename = os.path.basename(file)
    basename = os.path.splitext(basename)[0]

    return folder, basename

# ─── VARIABLES ──────────────────────────────────────────────────────────────────

mcherry_chnl = 2
gfp_chnl = 3

min_volume = 0.1

# ─── MAIN CODE ──────────────────────────────────────────────────────────────────

IJ.log("\\Clear")
IJ.log("Script starting...")

# Retrieve list of files
src_dir = str(src_dir)
files   = list_all_filenames(src_dir, filename_filter)

if files:

    file = files[0]

    folder, basename = get_file_info(file)

    csv_output = os.path.join(folder, "Results.csv")

            # Set up the DL environment
    command.run(CellposePrefsSet, False,
            "envType", "conda",
            "cellposeEnvDirectory", omnipose_env,
            "version", "1.0")

    filenames_list   = []
    series_list      = []
    obj_names_list   = []
    gfp_mean_int_list     = []
    mcherry_mean_int_list = []
    gfp_median_int_list     = []
    mcherry_median_int_list = []

    for index, file in enumerate(files):

        progressbar(index + 1, len(files), 1, "Opening : ")
        folder, basename = get_file_info(file)

        series_count = get_series_count_from_ome_metadata(file)

        for series in range(series_count):
            progressbar(series + 1, series_count, 2, "Opening series : ")

            imp = open_single_series_with_BF(file, series)

            cal          = imp.getCalibration()
            nbr_slices   = imp.getNSlices()
            nbr_channels = imp.getNChannels()
            nbr_frames   = imp.getNFrames()

            series_name = imp.getTitle()

            out_zip = os.path.join(folder, basename + "_series" + str(series + 1) + ".zip")

            bact_diameter = int(round(bact_diameter_calibrated / cal.pixelWidth))

            imp_mcherry = Duplicator().run(imp,
                                    mcherry_chnl, mcherry_chnl,
                                    1, imp.getNSlices(),
                                    1, imp.getNFrames())

            imp_bact = Duplicator().run(imp,
                                    mcherry_chnl, mcherry_chnl,
                                    1, imp.getNSlices(),
                                    1, imp.getNFrames())


            imp_gfp = Duplicator().run(imp,
                                    gfp_chnl, gfp_chnl,
                                    1, imp.getNSlices(),
                                    1, imp.getNFrames())

            IJ.run(imp_bact, "Median...", "radius=3 stack")

            progressbar(series + 1, len(files), 2, "Running Cellpose on series : ")

            res = command.run(Cellpose_SegmentImgPlusAdvanced, False,
                                "imp", imp_bact, "diameter", bact_diameter,
                                "model", "bact_fluor_omni", "nuclei_channel", 0,
                                "cyto_channel", 1,
                                "flow_threshold", 0,
                                "cellproba_threshold", 0,
                                "diam_threshold", 0,
                                "dimensionMode", "3D").get()
            imp_label_bact = res.getOutput("cellpose_imp")

            img_bact = ImageInt.wrap(imp_label_bact)
            pop_bact = Objects3DPopulation(img_bact)
            unit     = imp.getCalibration().getUnits()
            nb_bact  = pop_bact.getNbObjects()

            obj_to_remove = []
            IH_gfp = ImageHandler.wrap(imp_gfp)
            IH_mcherry = ImageHandler.wrap(imp_mcherry)

            obj_index = 1

            progressbar(series + 1, series_count, 2, "Measuring intensities on series : ")

            for i in range(nb_bact):
                obj = pop_bact.getObject(i)
                if(obj.getVolumeUnit() < min_volume):
                    obj_to_remove.append(obj)
                    continue
                if obj.touchBorders(img_bact, False):
                    obj_to_remove.append(obj)
                    continue
                obj.setName("Bact_" + str(obj_index))
                obj_index = obj_index + 1

                obj_names_list.append(obj.getName())
                gfp_mean_int_list.append(obj.getPixMeanValue(IH_gfp))
                mcherry_mean_int_list.append(obj.getPixMeanValue(IH_mcherry))
                gfp_median_int_list.append(obj.getPixMedianValue(IH_gfp))
                mcherry_median_int_list.append(obj.getPixMedianValue(IH_mcherry))

            for obj in obj_to_remove:
                pop_bact.removeObject(obj)

            pop_bact.saveObjects(out_zip)
            nb_bact = pop_bact.getNbObjects()

            filenames_list.extend([basename] * nb_bact)
            series_list.extend([series + 1] * nb_bact)

            imp.close()
            imp_gfp.close()
            imp_bact.close()
            imp_mcherry.close()
            imp_label_bact.close()

    with open(csv_output, 'wb') as f:
            writer = csv.writer(f, delimiter = ';')
            writer.writerow(
                ["Image name; Series number; Object Name; GFP mean intensity; GFP median intensity; mCherry mean intensity; mCherry median intensity"])
            writer.writerows(izip(filenames_list, series_list,
                                  obj_names_list, gfp_mean_int_list,
                                  gfp_median_int_list,
                                  mcherry_mean_int_list,
                                  mcherry_median_int_list))

IJ.log("SCRIPT FINISHED")