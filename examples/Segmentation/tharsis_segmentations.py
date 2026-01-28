# Import the run_workflow function from full_workflow
import os
from mineye.BayesianSegmentation.full_workflow import run_workflow

# %% Output directory
output_dir = "examples/Data/Segmentation_Output_Data"

# %%  20m resolution
bands = {
    "B4": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B04_20m.jp2",
    "B6": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B06_20m.jp2",
    "B7": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B07_20m.jp2",
    "B8A": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B8A_20m.jp2",
    "B11": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B11_20m.jp2",
    "B12": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B12_20m.jp2",
    "TCI": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_TCI_20m.jp2",
    "SCL": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_SCL_20m.jp2",
}
# %%  60m resolution
bands = {
    "B4": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/60m/merged_B04_60m.jp2",
    "B6": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/60m/merged_B06_60m.jp2",
    "B7": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/60m/merged_B07_60m.jp2",
    "B8A": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/60m/merged_B8A_60m.jp2",
    "B11": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/60m/merged_B11_60m.jp2",
    "B12": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/60m/merged_B12_60m.jp2",
    "TCI": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/60m/merged_TCI_60m.jp2",
    "SCL": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/60m/merged_SCL_60m.jp2",
}

# %% Tharsis ROI 1 configuration
bounds = (204498 , 4170995, 227828, 4187076)
n_classes = 6
beta_init = 30.0
n_iterations = 400
beta_jump_length = 0.1

run_workflow(
    bands=bands,
    bounds=bounds,
    n_classes=n_classes,
    beta_init=beta_init,
    n_iterations=n_iterations,
    beta_jump_length=beta_jump_length,
    use_soil_mask=True,
    save_npy=False,
    plot_tci=True,
    ref_band="B4",
    output_prefix=os.path.join(output_dir, "ROI1"),
)

# %% Tharsis ROI 2 configuration
bounds = (259547, 4211257 , 289534, 4235224)
n_classes = 6
beta_init = 30.0
n_iterations = 400
beta_jump_length = 0.1

run_workflow(
    bands=bands,
    bounds=bounds,
    n_classes=n_classes,
    beta_init=beta_init,
    n_iterations=n_iterations,
    beta_jump_length=beta_jump_length,
    use_soil_mask=True,
    save_npy=False,
    plot_tci=True,
    ref_band="B4",
    output_prefix=os.path.join(output_dir, "ROI2"),
)

# %% Tharsis ROI 3 configuration
bounds = (299348, 4185186 , 324839 ,4201160)
n_classes = 6
beta_init = 30.0
n_iterations = 400
beta_jump_length = 0.1

run_workflow(
    bands=bands,
    bounds=bounds,
    n_classes=n_classes,
    beta_init=beta_init,
    n_iterations=n_iterations,
    beta_jump_length=beta_jump_length,
    use_soil_mask=True,
    save_npy=False,
    plot_tci=True,
    ref_band="B4",
    output_prefix=os.path.join(output_dir, "ROI3"),
)

# %% Tharsis ROI 4 configuration
bounds = (309100, 4211465 , 394200, 4285027) #xmax was reduced because of area coverage. need to fix null value handling!
n_classes = 6
beta_init = 30.0
n_iterations = 400
beta_jump_length = 0.1

run_workflow(
    bands=bands,
    bounds=bounds,
    n_classes=n_classes,
    beta_init=beta_init,
    n_iterations=n_iterations,
    beta_jump_length=beta_jump_length,
    use_soil_mask=True,
    save_npy=False,
    plot_tci=True,
    ref_band="B4",
    output_prefix=os.path.join(output_dir, "ROI4"),
)

# %% Tharsis ROI 5 configuration


bounds = (247900, 4170000 , 265675, 4201751)
n_classes = 6
beta_init = 30.0
n_iterations = 400
beta_jump_length = 0.1

run_workflow(
    bands=bands,
    bounds=bounds,
    n_classes=n_classes,
    beta_init=beta_init,
    n_iterations=n_iterations,
    beta_jump_length=beta_jump_length,
    use_soil_mask=True,
    save_npy=False,
    plot_tci=True,
    ref_band="B4",
    output_prefix=os.path.join(output_dir, "ROI5"),
)