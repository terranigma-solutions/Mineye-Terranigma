# Import the run_workflow function from full_workflow
from Segmentation.full_workflow import run_workflow


# %% Tharsis ROI 1 configuration
bands = {
    "B4": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R20m/T29SQB_20230829T110621_B04_20m.jp2",
    "B6": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R20m/T29SQB_20230829T110621_B06_20m.jp2",
    "B7": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R20m/T29SQB_20230829T110621_B07_20m.jp2",
    "B8A": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R20m/T29SQB_20230829T110621_B8A_20m.jp2",
    "B11": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R20m/T29SQB_20230829T110621_B11_20m.jp2",
    "B12": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R20m/T29SQB_20230829T110621_B12_20m.jp2",
    "TCI": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R20m/T29SQB_20230829T110621_TCI_20m.jp2",
    "SCL": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R20m/T29SQB_20230829T110621_SCL_20m.jp2",
}

bounds = (733891.6, 4168988.9, 756038.65, 4186614.52)
n_classes = 4
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
    output_prefix="ROI1",
)
