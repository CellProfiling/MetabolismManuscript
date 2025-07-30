import os.path
import re

import pandas as pd # v 2.1.4
import scipy.stats # v 1.11.4
import numpy as np # v 1.26.3
import skimage as ski # v 0.22.0
from imageio import imread # v 2.33.1


pd.set_option('display.max_columns', None)


def create_output_template(number_of_cells):
    # This function creates an empty dataframe, based on the number of cells detected in the images. This dataframe can then be populated by other functions.
    cells = pd.DataFrame(index=range(number_of_cells), columns=['Cell ID'
                                                                ]
                         )
    return cells


def cell_labels(output_df,img_label):
    # Generates and assigns the labels for individual cells. Each label starts with the image ID, passed as a string, and ends in 3 digits showing the number of the cell [mask] in the image.
    for i in range(output_df.shape[0]):
        output_df.at[i, 'Cell ID'] = img_label + "_" + "{:03d}".format(i)


def average_intensity_75(channel,mask):
    # Calculates the average intensity of the upper quartile of the channel within the mask. The mask is a 2-column array, with column 0 showing X-values and column 1 showing corresponding Y-values of the pixels.
    array = channel[mask[0],mask[1]]
    try:
        upper_quart = np.quantile(array,0.75)
        ave_int_75 = array[array>upper_quart].mean()
    except IndexError:
        ave_int_75 = np.nan
    return ave_int_75


def find_cell(cell_mask,nuclei_mask, i):
    # Finds a specific cell: returns 2 NumPy arrays; first contains the rows (y-coordinates), second - the columns (x-coordinates)
    cell_i = np.where(cell_mask==i+1)
    nucleus_i = np.where(nuclei_mask==i+1)
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


def cell_morphology(output_df,cell_mask,nuclei_mask):
    # Assesses cell morphology: cell size (in pixels), nucleus size (in pixels), nucleus-to-cytoplasm ratio, whether the whole cell is in the image (False for the cell masks touching the borders of the image, True for all other cell masks).
    for i in range(output_df.shape[0]):
        cell_i, nucleus_i, cytoplasm_i = find_cell(cell_mask,nuclei_mask,i)
        output_df.at[i,'Cell size'] = len(cell_i[0])
        output_df.at[i,'Nucleus size'] = len(nucleus_i[0])
        output_df.at[i, 'NCR'] = len(nucleus_i[0]) / len(cytoplasm_i[0])
        if np.isin(0,cell_i) or np.isin((len(cell_mask)-1),cell_i):
            output_df.at[i, 'Is complete'] = False
        else:
            output_df.at[i, 'Is complete'] = True


def calculate_average_signal_top25(round, dapi, dapi_zmean, dapi_zstd, tubulin, tubulin_zmean, tubulin_zstd, marker, marker_zmean, marker_zstd, hpa, hpa_zmean, hpa_zstd, cell_mask, nuclei_mask, output_df):
    # Calculates the average fluorescence intensity for the selected channels in the whole cell, the nucleus and the cytoplasm, as well as the nucleus-to-cytoplasm ratio of the signal.
    position = re.search(r'r\d\dc\d\df\d\d', dapi[0]).group()

    dapi_img = ski.io.imread(dapi[0])
    dapi_img = zscore_img(dapi_img,cell_mask,dapi_zmean,dapi_zstd)
    tubulin_img = ski.io.imread(tubulin[0])
    tubulin_img = zscore_img(tubulin_img,cell_mask,tubulin_zmean,tubulin_zstd)
    marker_img = ski.io.imread(marker[0])
    marker_img = zscore_img(marker_img,cell_mask,marker_zmean,marker_zstd)
    hpa_img = ski.io.imread(hpa[0])
    hpa_img = zscore_img(hpa_img,cell_mask,hpa_zmean,hpa_zstd)

    fov_cells = output_df[output_df["Cell ID"].str.contains(position)].reset_index(drop=True)


    for i in range(fov_cells.shape[0]):
        cell_i, nucleus_i, cytoplasm_i = find_cell(cell_mask, nuclei_mask, i)
        fov_cells.at[i, round+'_'+'Average intensity of DAPI'] = average_intensity_75(dapi_img, cell_i)
        fov_cells.at[i, round+'_'+'Nucleus_DAPI_avg'] = average_intensity_75(dapi_img, nucleus_i)
        fov_cells.at[i, round+'_'+'Cytoplasm_DAPI_avg'] = average_intensity_75(dapi_img, cytoplasm_i)
        fov_cells.at[i, round+'_'+'Nucleus-to-cytoplasm DAPI intensity'] = fov_cells.loc[i, round+'_'+'Nucleus_DAPI_avg'] / fov_cells.loc[i, round+'_'+'Cytoplasm_DAPI_avg']

        fov_cells.at[i, round+'_'+'Average intensity of tubulin'] = average_intensity_75(tubulin_img, cell_i)
        fov_cells.at[i, round+'_'+'Nucleus_tubulin_avg'] = average_intensity_75(tubulin_img,nucleus_i)
        fov_cells.at[i, round+'_'+'Cytoplasm_tubulin_avg'] = average_intensity_75(tubulin_img,cytoplasm_i)
        fov_cells.at[i, round+'_'+'Nucleus-to-cytoplasm tubulin intensity'] = fov_cells.loc[i,round+'_'+'Nucleus_tubulin_avg'] / fov_cells.loc[i, round+'_'+'Cytoplasm_tubulin_avg']

        fov_cells.at[i, round+'_'+'Average intensity of marker'] = average_intensity_75(marker_img, cell_i)
        fov_cells.at[i, round+'_'+'Nucleus_marker_avg'] = average_intensity_75(marker_img,nucleus_i)
        fov_cells.at[i, round+'_'+'Cytoplasm_marker_avg'] = average_intensity_75(marker_img,cytoplasm_i)
        fov_cells.at[i, round+'_'+'Nucleus-to-cytoplasm marker intensity'] = fov_cells.loc[i,round+'_'+'Nucleus_marker_avg'] / fov_cells.loc[i, round+'_'+'Cytoplasm_marker_avg']

        fov_cells.at[i, round+'_'+'Average intensity of HPA'] = average_intensity_75(hpa_img, cell_i)
        fov_cells.at[i, round+'_'+'Nucleus_HPA_avg'] = average_intensity_75(hpa_img,nucleus_i)
        fov_cells.at[i, round+'_'+'Cytoplasm_HPA_avg'] = average_intensity_75(hpa_img,cytoplasm_i)
        fov_cells.at[i, round+'_'+'Nucleus-to-cytoplasm HPA intensity'] = fov_cells.loc[i,round+'_'+'Nucleus_HPA_avg'] / fov_cells.loc[i, round+'_'+'Cytoplasm_HPA_avg']

    return fov_cells


def image_analysis(path, wells, fovs, files, rounds, dapi_channel, tubulin_channel, marker_channel, hpa_channel):
    _, masks = run_fast_scandir(path, ".png")

    for well in wells:
        output = create_output_template(0)
        for fov in fovs:
            cell_mask = imread([m for m in masks if (well in m) & (fov in m) & ("cellmask" in m) & ("round 01" in m)][0])
            nuclei_mask = imread(
                [m for m in masks if (well in m) & (fov in m) & ("nucleimask" in m) & ("round 01" in m)][0])
            # The maximum value present in the cell mask determines the total number of detected cells
            cells = create_output_template(cell_mask.max())
            cell_labels(cells, (well + fov[:-1]))
            cell_morphology(cells, cell_mask, nuclei_mask)
            output = pd.concat([output, cells]).reset_index(drop=True)

        for r in rounds:
                round_output = create_output_template(0)
                zcells = find_zcells(output)
                dapi_zmean, dapi_zstd = find_zfactor(path,fovs,zcells,"ch01",r)
                tubulin_zmean, tubulin_zstd = find_zfactor(path,fovs,zcells,"ch02",r)
                marker_zmean, marker_zstd = find_zfactor(path,fovs,zcells,"ch03",r)
                hpa_zmean, hpa_zstd = find_zfactor(path,fovs,zcells,"ch04",r)
                for fov in fovs:
                    images = [f for f in files if (well in f) & (fov in f)]
                    cell_mask = imread([m for m in masks if (well in m) & (fov in m) & ("cellmask" in m) & ("round 01" in m)][0])
                    nuclei_mask = imread([m for m in masks if (well in m) & (fov in m) & ("nucleimask" in m) & ("round 01" in m)][0])
                    dapi = [i for i in images if (r in i) & (dapi_channel in i) & (fov in i)]
                    tubulin = [i for i in images if (r in i) & (tubulin_channel in i) & (fov in i)]
                    marker = [i for i in images if (r in i) & (marker_channel in i) & (fov in i)]
                    hpa = [i for i in images if (r in i) & (hpa_channel in i) & (fov in i)]
                    print(dapi)
                    results=calculate_average_signal_top25(r.split(" ")[1], dapi, dapi_zmean, dapi_zstd, tubulin, tubulin_zmean, tubulin_zstd, marker, marker_zmean, marker_zstd, hpa, hpa_zmean, hpa_zstd, cell_mask, nuclei_mask, output)
                    round_output = pd.concat([round_output,results]).reset_index(drop=True)
                output=pd.merge(output,round_output,on=output.columns.values.tolist())

        output.to_csv(path+well+"_all rounds_cell by cell.csv")


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


def find_zcells(cells):
    zcells = cells.loc[(cells["Is complete"] == True) & (cells["Nucleus size"] > 1500)].sample(50)
    return zcells


def find_zfactor(path, fovs, cells, ch, r):
    _, masks = run_fast_scandir(path, ".png")
    _, channels = run_fast_scandir(path, ".tiff")
    zpixels = []
    for fov in fovs:
        c = [c for c in channels if (fov in c) & (ch in c) & (r in c) & ("aligned" in c)]
        if len(c)>0:
            channel = ski.io.imread(c[0])
            nuclei_mask = imread([m for m in masks if (fov in m) & ("nucleimask" in m) & ("round 01" in m)][0])
            cell_mask = imread([m for m in masks if (fov in m) & ("cellmask" in m) & ("round 01" in m)][0])
            for cell in cells["Cell ID"]:
                if fov.split("-")[0] in cell:
                    cell_id = int(cell[-3:])
                    cell_i, nucleus_i, cytoplasm_i = find_cell(cell_mask, nuclei_mask, cell_id)
                    array = channel[cell_i[0], cell_i[1]]
                    zpixels = np.append(zpixels, array)
            zmean = np.average(zpixels)
            zstd = np.std(zpixels)
        else:
            zmean = 0
            zstd = 0
    return zmean, zstd


def zscore_img(image, cell_mask, zmean, zstd):
    channel = ski.util.img_as_float(image)
    channel[cell_mask > 0] = scipy.stats.zscore(channel[cell_mask > 0], axis=None)
    channel[cell_mask > 0] = channel[cell_mask > 0] * zstd
    channel[cell_mask > 0] = channel[cell_mask > 0] + zmean
    return channel


path = "X:\CELL_PROFILING\\research\Alina\\4i\\4i_U2OS_set01\\"
full_paths, files = run_fast_scandir(path, ".tiff")
files = [f for f in files if ("aligned" in f) & ("elution" not in f)]
fovs = set([re.search(r'f\d\d-', f).group() for f in files])
wells = set([re.search(r'r\d\dc\d\d', f).group() for f in files])
rounds = set([re.search(r'round \d\d', f).group() for f in full_paths])

dapi_channel = "ch01"
tubulin_channel = "ch02"
marker_channel = "ch03"
hpa_channel = "ch04"


image_analysis(path, wells, fovs, files, rounds, dapi_channel, tubulin_channel, marker_channel, hpa_channel)
