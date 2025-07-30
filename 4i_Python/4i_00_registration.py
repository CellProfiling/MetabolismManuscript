import os
import numpy as np # v 1.26.3
import cv2 as cv # opencv-python v 4.9.0.80
import re
import hpacellseg.cellsegmentator as cellsegmentator # v 0.1.8, available at https://github.com/CellProfiling/HPA-Cell-Segmentation
from hpacellseg.utils import label_cell, label_nuclei
from imageio import imwrite # v 2.33.1


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


def register(reference,input):
    ref = reference
    img = input

    #  Resize the image by a factor of 8 on each side. If your images are
    # very high-resolution, you can try to resize even more, but if they are
    # already small you should set this to something less agressive.
    resize_factor = 1.0/8.0

    img_rs = cv.resize(img, (0,0), fx=resize_factor, fy=resize_factor)
    ref_rs = cv.resize(ref, (0,0), fx=resize_factor, fy=resize_factor)
    img_rs_8 = cv.convertScaleAbs(img_rs)
    ref_rs_8 = cv.convertScaleAbs(ref_rs)

    # Initiate SIFT detector
    sift_detector = cv.SIFT_create()

    # Find the keypoints and descriptors with SIFT on the lower resolution images
    kp1, des1 = sift_detector.detectAndCompute(img_rs_8, None)
    kp2, des2 = sift_detector.detectAndCompute(ref_rs_8, None)

    # BFMatcher with default params
    bf = cv.BFMatcher()
    matches = bf.knnMatch(des1, des2, k=2)

    # Filter out poor matches
    good_matches = []
    for m,n in matches:
        if m.distance < 0.75*n.distance:
            good_matches.append(m)

    matches = good_matches

    points1 = np.zeros((len(matches), 2), dtype=np.float32)
    points2 = np.zeros((len(matches), 2), dtype=np.float32)

    for i, match in enumerate(matches):
        points1[i, :] = kp1[match.queryIdx].pt
        points2[i, :] = kp2[match.trainIdx].pt

    # Find homography
    try:
        H, mask = cv.findHomography(points1, points2, cv.RANSAC)
        # Get low-res and high-res sizes
        low_height, low_width = img_rs.shape
        height, width = img.shape
        low_size = np.float32([[0, 0], [0, low_height], [low_width, low_height], [low_width, 0]])
        high_size = np.float32([[0, 0], [0, height], [width, height], [width, 0]])

        # Compute scaling transformations
        scale_up = cv.getPerspectiveTransform(low_size, high_size)
        scale_down = cv.getPerspectiveTransform(high_size, low_size)

        #  Combine the transformations. Remember that the order of the transformation
        # is reversed when doing matrix multiplication
        # so this is actualy scale_down -> H -> scale_up
        h_and_scale_up = np.matmul(scale_up, H)
        scale_down_h_scale_up = np.matmul(h_and_scale_up, scale_down)

        return scale_down_h_scale_up
    except cv.error:
        print("Could not register the images")
        return input
    except ValueError:
        print("Could not register the images")
        return input




def warp(reference, image, function):

    # Warp image 1 to align with image 2
    img1Reg = cv.warpPerspective(
                image,
                function,
                (reference.shape[1], reference.shape[0])
              )
    return img1Reg


def segment(tubulin, dapi):
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

    image_path = os.path.dirname(dapi)
    image_name = os.path.basename(dapi).split("_")[0]

    if not os.path.exists(image_path+"\\"+image_name+"_nucleimask.png"):

        tub_img = (cv.normalize(cv.imread(tubulin,-1),None, alpha=0, beta=65535, norm_type=cv.NORM_MINMAX, dtype=cv.CV_16U)).astype(np.uint16)
        cv.imwrite(tubulin,tub_img)
        dapi_img = (cv.normalize(cv.imread(dapi,-1), None, alpha=0, beta=65535, norm_type=cv.NORM_MINMAX, dtype=cv.CV_16U)).astype(np.uint16)
        cv.imwrite(dapi,dapi_img)


        # For nuclei: taking in nuclei channels as inputs
        nuc_segmentations = segmentator.pred_nuclei(images[2])


        # For full cells: taking in 2 channels as inputs
        cell_segmentations = segmentator.pred_cells(images)

        # post-processing nuclei mask
        nuclei_mask = label_nuclei(nuc_segmentations[0])

        # post-processing nuclei and cell mask
        for i, (nuc_segmentation, cell_segmentation) in enumerate(zip(nuc_segmentations, cell_segmentations)):
            nuclei_mask, cell_mask = label_cell(nuc_segmentation, cell_segmentation)
            imwrite(image_path+"\\"+image_name+"_nucleimask.png", nuclei_mask)
            imwrite(image_path+"\\"+image_name+"_cellmask.png", cell_mask)



folder = "X:\CELL_PROFILING\\research\Alina\\4i\\4i_U2OS_set01\\"
reference_channel = "ch02" # channel for tubulin or DAPI that will be used for registration

# Find all relevant files in the directory (including subdirectories)
subfolders, files = run_fast_scandir(folder, [".tiff"])
_, masks = run_fast_scandir(folder, [".png"])
files = [f for f in files if ("maxproj" in f) & ("elution" not in f) & ("corrected" in f)]
masks = [m for m in masks if ("mask" in m)]

# Find the relevant fields of view
fovs = set([re.search(r'f\d\d-', f).group() for f in files])

# Find the relevant wells
wells = set([re.search(r'r\d\dc\d\d', f).group() for f in files])


for well in wells:
    for fov in fovs:
        images = [f for f in files if (well in f) & (fov in f)]
        channels = set([re.search(r'ch\d\d', f).group() for f in images])
        rounds = set([re.search(r'round \d\d', f).group() for f in images])
        reference = cv.imread(str([i for i in images if (reference_channel in i) & ("round 01" in i)][0]),-1)
        for r in rounds:
            tubulin = cv.imread(str([i for i in images if (reference_channel in i) & (r in i)][0]),-1)
            warp_function = register(reference, tubulin)
            for c in channels:
                mask = [m for m in masks if (well in m) & (fov in m) & (r in m) & (c in m)]
                aligned = [i for i in images if (c in i) & (r in i) & ("_aligned" in i)]
                if (len(mask)==0) & (len(aligned)==0):
                    try:
                        image = cv.imread(str([i for i in images if ("aligned" not in i) & (c in i) & (r in i)][0]),-1)
                        image_reg = warp(reference,image,warp_function)
                        cv.imwrite(str([i for i in images if ("aligned" not in i) & (c in i) & (r in i)][0][0:-5]+'_aligned.tiff'),
                                       image_reg)
                    except IndexError:
                        print("No channel")
                    except cv.error:
                        print("Could not register the images")



    full_paths_well, files_well = run_fast_scandir(folder, ".tiff")
    files_well = [f for f in files_well if ("aligned" in f) & ("elution" not in f) & (well in f)]
    dapi_channel = "ch01"
    tubulin_channel = "ch02"

    for fov in fovs:
        for f in files_well:
            if (fov in f) & ("round 01" in f) & (dapi_channel in f):
                dapi = f
            elif (fov in f) & ("round 01" in f) & (tubulin_channel in f):
                tubulin = f
        segment(tubulin,dapi)