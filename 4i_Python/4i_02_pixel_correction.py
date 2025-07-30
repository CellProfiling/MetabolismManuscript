import os.path

import hpacellseg.cellsegmentator as cellsegmentator # v 0.1.8, available at https://github.com/CellProfiling/HPA-Cell-Segmentation
import pandas as pd # v 2.1.4
import re
import scipy.stats # v 1.11.4
from hpacellseg.utils import label_cell, label_nuclei
from imageio import imwrite, imread # v 2.33.1
import numpy as np # v 1.26.3
import skimage as ski # v 0.22.0



pd.set_option('display.max_columns', None)


def segment(tubulin, dapi, full_path):
    # This function takes two paths to two channels of the same image (tubulin and DAPI), performs segmentation of those, and returns the nuclei mask and the cell mask. full_path is used to write the png files of the masks.

    NUC_MODEL = "./nuclei-model.pth"
    CELL_MODEL = "./cell-model.pth"
    segmentator = cellsegmentator.CellSegmentator(
        NUC_MODEL,
        CELL_MODEL,
        scale_factor=0.5,
        device="cuda",
        # NOTE: setting padding=True seems to solve most issues that have been encountered
        #       during our single cell Kaggle challenge.
        padding=True,
        multi_channel_model=False,
    )
    images = [
        [tubulin],
        None,
        [dapi]
    ]

    # For nuclei: taking in nuclei channels as inputs
    nuc_segmentations = segmentator.pred_nuclei(images[2])
    # npns = np.array(nuc_segmentations)

    # For full cells: taking in 2 channels as inputs
    cell_segmentations = segmentator.pred_cells(images)

    # post-processing nuclei mask
    nuclei_mask = label_nuclei(nuc_segmentations[0])

    # post-processing nuclei and cell mask
    for i, (nuc_segmentation, cell_segmentation) in enumerate(zip(nuc_segmentations, cell_segmentations)):
        nuclei_mask, cell_mask = label_cell(nuc_segmentation, cell_segmentation)
        imwrite(full_path + "ch01t01_nucleimask.png", nuclei_mask)
        imwrite(full_path + "ch01t01_cellmask.png", cell_mask)
    return nuclei_mask, cell_mask


def create_output_template(number_of_cells):
    # This function creates an empty dataframe, based on the number of cells detected in the images. This dataframe can then be populated by other functions.
    cells = pd.DataFrame(index=range(number_of_cells), columns=['Cell ID',
                                                                'Cell size',
                                                                'Nucleus size',
                                                                'NCR',
                                                                'Is complete'
                                                                ]
                         )
    return cells


def cell_labels(output_df, img_label):
    # Generates and assigns the labels for individual cells. Each label starts with the image ID, passed as a string, and ends in 3 digits showing the number of the cell [mask] in the image.
    for i in range(output_df.shape[0]):
        output_df.at[i, 'Cell ID'] = img_label + "_" + "{:03d}".format(i)

def find_cell(cell_mask, nuclei_mask, i):
    # Finds a specific cell: returns 2 NumPy arrays; first contains the rows (y-coordinates), second - the columns (x-coordinates)
    cell_i = np.where(cell_mask == i + 1)
    nucleus_i = np.where(nuclei_mask == i + 1)
    zipped_cell = list(zip(cell_i[0], cell_i[1]))
    zipped_nucleus = list(zip(nucleus_i[0], nucleus_i[1]))
    zipped_cytoplasm = list(set(zipped_cell) - set(zipped_nucleus))
    try:
        cytoplasm_row, cytoplasm_col = zip(*zipped_cytoplasm)
        cytoplasm_i = [np.array(cytoplasm_row), np.array(cytoplasm_col)]
    except ValueError:
        print("ValueError: cell size and nucleus size are identical. Whole cell is assigned as cytoplasm.")
        cytoplasm_i = cell_i
    return cell_i, nucleus_i, cytoplasm_i


def cell_morphology(output_df, cell_mask, nuclei_mask):
    # Assesses cell morphology: cell size (in pixels), nucleus size (in pixels), nucleus-to-cytoplasm ratio, whether the whole cell is in the image (False for the cell masks touching the borders of the image, True for all other cell masks).
    for i in range(output_df.shape[0]):
        cell_i, nucleus_i, cytoplasm_i = find_cell(cell_mask, nuclei_mask, i)
        output_df.at[i, 'Cell size'] = len(cell_i[0])
        output_df.at[i, 'Nucleus size'] = len(nucleus_i[0])
        output_df.at[i, 'NCR'] = len(nucleus_i[0]) / len(cytoplasm_i[0])
        if np.isin(0, cell_i) or np.isin((len(cell_mask) - 1), cell_i):
            output_df.at[i, 'Is complete'] = False
        else:
            output_df.at[i, 'Is complete'] = True

def image_analysis(path, plate, well, good_FOVs):
    output = create_output_template(0)
    for fov in good_FOVs["FOV"]:
        full_path = path + well + "\\" + fov + "-"
        dapi = full_path + "ch02t01_maxproj_corrected_aligned.tiff"
        tubulin = full_path + "ch03t01_maxproj_corrected_aligned.tiff"

        if (os.path.exists(path + well + "\\" + fov + "-ch01t01_" + "nucleimask.png")):
            nuclei_mask = imread(path + well + "\\" + fov + "-ch01t01_" + "nucleimask.png")
            cell_mask = imread(path + well + "\\" + fov + "-ch01t01_" + "cellmask.png")

        else:
            nuclei_mask, cell_mask = segment(tubulin, dapi, full_path)

        # The maximum value present in the cell mask determines the total number of detected cells
        cells = create_output_template(cell_mask.max())

        cell_labels(cells, (plate + "_" + fov))
        cell_morphology(cells, cell_mask, nuclei_mask)

        output = pd.concat([output, cells]).reset_index(drop=True)

    return output


def find_zcells(cells):
    zcells = cells.loc[(cells["Is complete"] == True) & (cells["Nucleus size"] > 1500)].sample(100)
    return zcells


def find_zfactor(cells, ch):
    well = cells[("Cell ID")].str.extract(r'(r\d+c\d+)').drop_duplicates().reset_index(drop=True).values
    fovs = cells[("Cell ID")].str.extract(r'(f\d+)').drop_duplicates().reset_index(drop=True).values
    zpixels = []
    for fov in fovs:
        channel = ski.io.imread(
            path + str(well[0][0]) + "\\" + str(well[0][0]) + str(fov[0]) + "-" + ch + "t01_maxproj_corrected_aligned.tiff")
        nuclei_mask = imread(path + str(well[0][0]) + "\\" + str(well[0][0]) + str(fov[0]) + "-ch01t01_" + "nucleimask.png")
        cell_mask = imread(path + str(well[0][0]) + "\\" + str(well[0][0]) + str(fov[0]) + "-ch01t01_" + "cellmask.png")
        for cell in cells["Cell ID"]:
            if fov[0] in cell:
                cell_id = int(cell[-3:])
                cell_i, nucleus_i, cytoplasm_i = find_cell(cell_mask, nuclei_mask, cell_id)
                array = channel[cell_i[0], cell_i[1]]
                zpixels = np.append(zpixels, array)
    zmean = np.average(zpixels)
    zstd = np.std(zpixels)
    return zmean, zstd


def zscore_img(impath, zmean, zstd, ch):
    channel = ski.util.img_as_float(ski.io.imread(impath + "-" + ch + "t01_maxproj_corrected_aligned.tiff"))
    cell_mask = ski.io.imread(impath + "-ch01t01_" + "cellmask.png")
    channel[cell_mask > 0] = scipy.stats.zscore(channel[cell_mask > 0], axis=None)
    channel[cell_mask > 0] = channel[cell_mask > 0] * zstd
    channel[cell_mask > 0] = channel[cell_mask > 0] + zmean
    # plt.imshow(channel)
    # plt.show()
    imwrite(impath + "-" + ch + "t01_z-scored.tiff", channel)

def run_fast_scandir(dir, ext):    # dir: str, ext: list
    subfolders, files = [], []

    for f in os.scandir(dir):
        if f.is_dir():
            subfolders.append(f.path)
        if f.is_file():
            if os.path.splitext(f.name)[1].lower() in ext:
                files.append(f.path)


    for dir in list(subfolders):
        sf, f = run_fast_scandir(dir, ext)
        subfolders.extend(sf)
        files.extend(f)
    return subfolders, files




path = "X:\CELL_PROFILING\\research\Alina\\4i\\4i_U2OS_set01\\U2OS_4i_round 01_00_stained_250312\hs\\2d366600-ce92-458c-b148-066badbd0aca\images\\"
plate = "U2OS_4i_round 01_00_stained_250312"
_, files = run_fast_scandir(path, [".png"])
wells = set([re.search(r'r\d\dc\d\d', f).group() for f in files])
channels = ["ch01", "ch02", "ch03", "ch04"]

for well in wells:
    good_FOVs = pd.DataFrame(['01', '02', '03', '04', '05', '06', '07', '08', '09'], columns=['FOV'])
    good_FOVs['FOV'] = well + 'f' + good_FOVs['FOV'].astype(str)

    analysis_output = image_analysis(path, plate, well, good_FOVs)

    for channel in channels:
        zcells = find_zcells(analysis_output) # randomly select a subset of cells to use for z-scoring
        zmean, zstd = find_zfactor(zcells, channel) # use the subset to find the z-factor for the whole set of fields of view
        print(zmean)
        print(zstd)
        for fov in good_FOVs['FOV']:
            impath = path + well + "\\" + fov
            if os.path.isfile(impath+ "-" + channel + "t01_z-scored.tiff")==False:
                zscore_img(impath, zmean, zstd, channel) # use the calculated z-factor to normalize the intensities across fields of view
            else:
                print(well+"_"+fov+"_"+channel+": z-scored image exists")
