import os
import numpy as np
import rasterio
from rasterio.io import MemoryFile
from scipy.signal import savgol_filter
from sklearn.decomposition import PCA
from sklearn.covariance import EmpiricalCovariance
import xml.etree.ElementTree as ET


# ============================================================
# 1. Parse wavelengths from ENMAP metadata XML
# ============================================================

def parse_enmap_wavelengths(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    wavelengths = []

    for band in root.findall(".//BAND"):
        wl_el = band.find("WAVELENGTH")
        if wl_el is not None:
            wavelengths.append(float(wl_el.text))

    return np.array(wavelengths)  # nm


# ============================================================
# 2. Load ENMAP cube, metadata, QA masks
# ============================================================

def load_enmap_dataset(path):
    files = os.listdir(path)

    cube_path = next(f for f in files if "SPECTRAL_IMAGE" in f and f.endswith(".TIF"))
    cube_path = os.path.join(path, cube_path)

    with rasterio.open(cube_path) as src:
        cube = src.read()      # (bands, y, x)
        meta = src.meta.copy()

    xml_path = next(f for f in files if "METADATA.XML" in f.upper())
    xml_path = os.path.join(path, xml_path)

    qa_masks = {}
    for f in files:
        if ("QL_PIXELMASK" in f) or ("QL_QUALITY" in f):
            fp = os.path.join(path, f)
            with rasterio.open(fp) as src:
                qa_masks[f] = src.read(1)

    return cube, meta, qa_masks, xml_path


# ============================================================
# 3. QA mask
# ============================================================

def build_mask(qa_masks):
    mask = None
    for arr in qa_masks.values():
        bad = (arr != 0)
        mask = bad if mask is None else (mask | bad)
    return mask


# ============================================================
# 4. Drop noisy ENMAP bands
# ============================================================

def drop_bad_bands(cube, wavelengths):
    keep = []
    for i, wl in enumerate(wavelengths):
        if wl < 430: continue
        if 1340 < wl < 1440: continue
        if 1800 < wl < 2000: continue
        if wl > 2450: continue
        keep.append(i)

    return cube[keep], wavelengths[keep], keep


# ============================================================
# 5. Savitzky–Golay smoothing
# ============================================================

def smooth_cube(cube, window=11, poly=2):
    return savgol_filter(
        cube, axis=0, mode="mirror",
        window_length=window, polyorder=poly
    )


# ============================================================
# 6. MNF transform
# ============================================================

def estimate_noise(cube):
    diffs = cube[:, 1:, :] - cube[:, :-1, :]
    diffs = diffs.reshape(cube.shape[0], -1).T
    cov = EmpiricalCovariance().fit(diffs)
    return cov.covariance_

def run_mnf(cube, n_components=12):
    bands, h, w = cube.shape
    X = cube.reshape(bands, -1).T

    noise_cov = estimate_noise(cube)
    noise_inv = np.linalg.inv(noise_cov)

    Xw = X @ noise_inv

    pca = PCA(n_components=n_components)
    mnf = pca.fit_transform(Xw)

    return mnf.T.reshape(n_components, h, w)


# ============================================================
# 7. Absorption depths
# ============================================================

def continuum_removed_depth(cube, wavelengths, center_range):
    start, end = center_range
    idx = np.where((wavelengths >= start) & (wavelengths <= end))[0]
    if len(idx) < 3:
        return None

    sub = cube[idx]
    left = sub[0]
    right = sub[-1]
    continuum = np.linspace(left, right, sub.shape[0])
    cr = 1 - (sub / continuum)
    return np.nanmax(cr, axis=0)


# ============================================================
# 8. Derivative PCA
# ============================================================

def derivative_pca(cube, wavelengths, n_components=3):
    deriv = np.gradient(cube, axis=0)
    bands, h, w = deriv.shape

    X = deriv.reshape(bands, -1).T
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(X)

    return pcs.T.reshape(n_components, h, w)


# ============================================================
# 9. Wrap 2D array into an in-memory rasterio dataset
# ============================================================

def create_rasterio_layer(array, meta):
    """Return a rasterio in-memory DatasetReader for a numpy array."""
    array = array.astype("float32")
    new_meta = meta.copy()
    new_meta.update({
        "count": 1,
        "dtype": "float32",
    })

    memfile = MemoryFile()
    with memfile.open(**new_meta) as dst:
        dst.write(array, 1)

    return memfile.open()   # <-- rasterio.io.DatasetReader


# ============================================================
# 10. Main pipeline
# ============================================================

def enmap_to_feature_stack(folder_path, n_mnf=12, n_deriv_pcs=3):
    cube, meta, qa_masks, xml = load_enmap_dataset(folder_path)

    wavelengths = parse_enmap_wavelengths(xml)

    mask = build_mask(qa_masks)

    cube, wavelengths, _ = drop_bad_bands(cube, wavelengths)

    cube = smooth_cube(cube)
    cube[:, mask] = np.nan

    mnf = run_mnf(cube, n_components=n_mnf)

    depth_2200 = continuum_removed_depth(cube, wavelengths, (2200, 2230))
    depth_2300 = continuum_removed_depth(cube, wavelengths, (2300, 2330))
    depth_2340 = continuum_removed_depth(cube, wavelengths, (2320, 2350))
    depth_1000 = continuum_removed_depth(cube, wavelengths, (900, 1030))
    depth_1750 = continuum_removed_depth(cube, wavelengths, (1700, 1800))

    deriv = derivative_pca(cube, wavelengths, n_components=n_deriv_pcs)

    features = {}

    # MNF components
    for i in range(mnf.shape[0]):
        features[f"MNF_{i+1:02d}"] = create_rasterio_layer(mnf[i], meta)

    # Absorption depths
    if depth_2200 is not None: features["Depth_2200"] = create_rasterio_layer(depth_2200, meta)
    if depth_2300 is not None: features["Depth_2300"] = create_rasterio_layer(depth_2300, meta)
    if depth_2340 is not None: features["Depth_2340"] = create_rasterio_layer(depth_2340, meta)
    if depth_1000 is not None: features["Depth_1000"] = create_rasterio_layer(depth_1000, meta)
    if depth_1750 is not None: features["Depth_1750"] = create_rasterio_layer(depth_1750, meta)

    # Derivative PCA
    for i in range(n_deriv_pcs):
        features[f"Deriv_PC{i+1}"] = create_rasterio_layer(deriv[i], meta)

    return features, meta