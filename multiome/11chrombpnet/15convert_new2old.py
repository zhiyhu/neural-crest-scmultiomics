# reticulate::use_condaenv("tfmodiscolit")
# Zhiyuan
# 4 jan 2023
import glob
import modiscolite.io

# Define the directory and file pattern
directory = "/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/10tfmodisco"
pattern = "mNC_head_mesenchymal"

# Construct the wildcard pattern for new files
new_file_pattern = f"{directory}/*/modisco_results_profile.h5"

# Iterate over files that match the pattern
for new_h5 in glob.glob(new_file_pattern):
    # Construct the corresponding old file path
    old_h5 = new_h5.replace("10tfmodisco", "15cwm_scan").replace("/modisco_results_profile.h5", ".modisco.hdf5")

    # Perform the conversion
    print(f"Converting {new_h5} to {old_h5}")
    modiscolite.io.convert_new_to_old(new_h5, old_h5)



