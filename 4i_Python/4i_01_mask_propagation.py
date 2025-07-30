import os
import re
import shutil

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

folder = "X:\CELL_PROFILING\\research\Alina\\4i\\4i_U2OS_set01\\"
reference_channel = "ch02" # channel for tubulin or DAPI that will be used for registration

# Find all relevant files in the directory (including subdirectories)
subfolders, files = run_fast_scandir(folder, [".png"])
files = [f for f in files if ("mask" in f)]


rounds = set([re.search(r'round \d\d', f).group() for f in files])


# Find the relevant fields of view
fovs = set([re.search(r'f\d\d-', f).group() for f in files])

# Find the relevant wells
wells = set([re.search(r'r\d\dc\d\d', f).group() for f in files])


for well in wells:
    well_folders = [s for s in subfolders if ("round" in s)& (well in s)&("round 01" not in s)]
    for fov in fovs:
        images = [f for f in files if (well in f) & (fov in f)]
        for i in images:
            i_name = os.path.basename(i)
            for wf in well_folders:
                shutil.copy(i,wf) # copies the existing cell and nuclei masks to the directories of the consecutive 4i rounds