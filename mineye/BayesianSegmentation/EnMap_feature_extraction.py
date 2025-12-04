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
    """Parse center wavelengths from an EnMAP metadata XML file.

    Preference order:
    1) bandCharacterisation/bandID(@number)/wavelengthCenterOfBand (sorted by @number)
    2) BAND-like elements with wavelength child tags
    3) Vector list tags like CENTER_WAVELENGTHS
    Returns wavelengths in nanometers (nm).
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    def strip(tag):
        return tag.split('}')[-1] if '}' in tag else tag

    def to_float_safe(text):
        if text is None:
            return None
        s = text.strip()
        if not s:
            return None
        s_clean = (
            s.replace("nm", " ")
             .replace("nanometer", " ")
             .replace("nanometre", " ")
             .replace("µm", " ")
             .replace("um", " ")
             .replace("micrometer", " ")
             .replace("micrometre", " ")
        )
        try:
            return float(s_clean.split()[0])
        except Exception:
            return None

    units_detected = None  # 'nm' or 'um'

    # Strategy 1: Explicit EnMAP structure: bandCharacterisation > bandID number > wavelengthCenterOfBand
    band_map = {}
    for elem in root.iter():
        if strip(elem.tag).lower() == "bandid":
            num_txt = elem.attrib.get("number")
            try:
                num = int(num_txt) if num_txt is not None else None
            except Exception:
                num = None
            wl = None
            for child in list(elem):
                if strip(child.tag).lower() == "wavelengthcenterofband":
                    wl = to_float_safe(child.text)
                    # Unit attr if present
                    unit_attr = (child.attrib.get("unit") or child.attrib.get("units") or "").lower()
                    if unit_attr:
                        if "nm" in unit_attr:
                            units_detected = units_detected or "nm"
                        elif "um" in unit_attr or "µm" in unit_attr or "microm" in unit_attr:
                            units_detected = units_detected or "um"
                    break
            if num is not None and wl is not None:
                band_map[num] = wl
    if len(band_map) > 0:
        # Sort by band number to ensure exactly described bands are used
        numbers = sorted(band_map.keys())
        wavelengths = [band_map[n] for n in numbers]
    else:
        wavelengths = []

    # Strategy 2: Generic BAND-like scan (only if Strategy 1 yielded nothing)
    if len(wavelengths) == 0:
        band_like_names = {"BAND", "BAND_INFORMATION", "SPECTRALBAND", "CHANNEL", "BANDID"}
        wavelength_like_subtags = (
            "WAVELENGTH", "CENTER_WAVELENGTH", "CENTRAL_WAVELENGTH", "WAVELENGTH_CENTER",
            "BAND_CENTER", "CENTER", "CWAVELENGTH", "WAVELENGTHCENTEROFBAND"
        )
        for elem in root.iter():
            tag_u = strip(elem.tag).upper()
            if tag_u in band_like_names:
                wl_val = None
                for child in list(elem):
                    ctag = strip(child.tag).upper()
                    if any(k in ctag for k in wavelength_like_subtags):
                        wl_val = to_float_safe(child.text)
                        unit_attr = (child.attrib.get("unit") or child.attrib.get("units") or "").lower()
                        if unit_attr:
                            if "nm" in unit_attr:
                                units_detected = units_detected or "nm"
                            elif "um" in unit_attr or "µm" in unit_attr or "microm" in unit_attr:
                                units_detected = units_detected or "um"
                        if wl_val is not None:
                            break
                if wl_val is None:
                    for k, v in elem.attrib.items():
                        ku = k.upper()
                        if any(x in ku for x in wavelength_like_subtags):
                            wl_val = to_float_safe(v)
                            unit_attr = (elem.attrib.get("unit") or elem.attrib.get("units") or "").lower()
                            if unit_attr:
                                if "nm" in unit_attr:
                                    units_detected = units_detected or "nm"
                                elif "um" in unit_attr or "µm" in unit_attr or "microm" in unit_attr:
                                    units_detected = units_detected or "um"
                            break
                if wl_val is not None:
                    wavelengths.append(wl_val)

    # Strategy 3: Vector-like lists (only if still empty)
    if len(wavelengths) == 0:
        vector_like_tags = (
            "WAVELENGTHS", "CENTER_WAVELENGTHS", "CENTRAL_WAVELENGTHS",
            "BAND_CENTERS", "WAVELENGTH_LIST"
        )
        for elem in root.iter():
            tag_u = strip(elem.tag).upper()
            if any(vtag == tag_u or vtag in tag_u for vtag in vector_like_tags):
                text = (elem.text or "").strip()
                if text:
                    parts = text.replace(",", " ").replace(";", " ").split()
                    vals = []
                    for p in parts:
                        f = to_float_safe(p)
                        if f is not None:
                            vals.append(f)
                    if len(vals) > 0:
                        wavelengths = vals
                        unit_attr = (elem.attrib.get("unit") or elem.attrib.get("units") or "").lower()
                        if unit_attr:
                            if "nm" in unit_attr:
                                units_detected = units_detected or "nm"
                            elif "um" in unit_attr or "µm" in unit_attr or "microm" in unit_attr:
                                units_detected = units_detected or "um"
                        break

    wls = np.array(wavelengths, dtype=float) if len(wavelengths) > 0 else np.array([], dtype=float)

    if wls.size == 0:
        print(f"[EnMap] parse_enmap_wavelengths: No wavelengths parsed from {os.path.basename(xml_path)}")
        return wls

    converted = False
    if units_detected == "um" or (units_detected is None and np.nanmax(wls) < 10.0):
        wls = wls * 1000.0
        converted = True

    print(f"[EnMap] Wavelengths: count={wls.size}, min={np.nanmin(wls):.2f} nm, max={np.nanmax(wls):.2f} nm, converted_from_um={converted}")
    return wls  # in nm


# ============================================================
# 2. Load ENMAP cube, metadata, QA masks
# ============================================================

def load_enmap_dataset(path):
    files = os.listdir(path)

    cube_path = next(f for f in files if "SPECTRAL_IMAGE" in f and f.upper().endswith((".TIF", ".TIFF")))
    cube_path = os.path.join(path, cube_path)

    with rasterio.open(cube_path) as src:
        cube = src.read()      # (bands, y, x)
        meta = src.meta.copy()
    print(f"[EnMap] Loaded spectral cube: path={os.path.basename(cube_path)}, shape={cube.shape}, dtype={cube.dtype}")
    print(f"[EnMap] Geo: crs={meta.get('crs')}, transform={meta.get('transform')}, res={(abs(meta.get('transform').a), abs(meta.get('transform').e)) if meta.get('transform') is not None else None}")

    # Convert reflectance to float and apply scale if likely scaled integers
    if np.issubdtype(cube.dtype, np.integer):
        mx = float(np.nanmax(cube))
        mn = float(np.nanmin(cube))
        if 2.0 < mx <= 10000.0 and mn >= 0.0:
            cube = cube.astype(np.float32) / 10000.0
            print("[EnMap] Applied scale factor 1/10000 to convert integer reflectance to [0,1] float range.")
        else:
            cube = cube.astype(np.float32)
    else:
        cube = cube.astype(np.float32)

    xml_path = next(f for f in files if "METADATA.XML" in f.upper())
    xml_path = os.path.join(path, xml_path)
    print(f"[EnMap] Metadata XML: {os.path.basename(xml_path)}")

    qa_masks = {}
    for f in files:
        fu = f.upper()
        # Only consider actual raster files, ignore sidecar files like .aux.xml
        if ("QL_PIXELMASK" in fu or "QL_QUALITY" in fu) and fu.endswith((".TIF", ".TIFF")):
            fp = os.path.join(path, f)
            try:
                with rasterio.open(fp) as src:
                    qa_masks[f] = src.read(1)
                print(f"[EnMap] Loaded QA mask: {f} shape={qa_masks[f].shape} dtype={qa_masks[f].dtype}")
            except Exception as e:
                print(f"[EnMap] Skipped QA mask '{f}' due to error: {e}")
                continue

    if len(qa_masks) == 0:
        print("[EnMap] No QA mask rasters found in folder. Proceeding without QA mask.")

    return cube, meta, qa_masks, xml_path


# ============================================================
# 3. QA mask
# ============================================================

def build_mask(qa_masks):
    """Combine QA layers into a single boolean mask of bad pixels.

    Rules:
    - For generic QA layers (e.g., QL_PIXELMASK, QL_QUALITY_* except CLASSES): bad where value != 0
    - For QL_QUALITY_CLASSES: keep classes 0 (None) and 1 (Land); mark others (2=Water, 3=Background) as bad
    """
    mask = None
    for name, arr in qa_masks.items():
        name_u = name.upper()
        if "QUALITY_CLASSES" in name_u:
            # Keep 0 (None) and 1 (Land); everything else considered bad
            bad = ~np.isin(arr, [0, 1])
            # Log class distribution for diagnostics
            unique, counts = np.unique(arr, return_counts=True)
            stats = ", ".join([f"{int(u)}:{int(c)}" for u, c in zip(unique, counts)])
            print(f"[EnMap] QA '{name}' class histogram -> {stats}")
        else:
            # Default rule: 0=good, non-zero=bad
            bad = (arr != 0)
        count_bad = int(np.sum(bad))
        total = bad.size
        pct = (count_bad / total * 100.0) if total > 0 else 0.0
        print(f"[EnMap] QA '{name}': bad_pixels={count_bad} ({pct:.2f}%)")
        mask = bad if mask is None else (mask | bad)
    if mask is not None:
        overall_bad = int(np.sum(mask))
        total = mask.size
        pct = (overall_bad / total * 100.0) if total > 0 else 0.0
        print(f"[EnMap] Combined QA bad mask: bad_pixels={overall_bad} ({pct:.2f}%), shape={mask.shape}")
    return mask


# ============================================================
# 4. Drop noisy ENMAP bands
# ============================================================

def drop_bad_bands(cube, wavelengths):
    keep = []
    for i, wl in enumerate(wavelengths):
        if wl < 430:
            continue
        if 1340 < wl < 1440:
            continue
        if 1800 < wl < 2000:
            continue
        if wl > 2450:
            continue
        keep.append(i)

    if len(keep) == 0:
        print("[EnMap] drop_bad_bands: No bands remained after filtering.")
        if len(wavelengths) > 0:
            print(f"[EnMap] Wavelengths stats -> count={len(wavelengths)}, min={np.min(wavelengths):.2f}, max={np.max(wavelengths):.2f}")
        else:
            print("[EnMap] Wavelengths array is empty.")
        raise ValueError("No valid bands after filtering. Check wavelength parsing/units and filter ranges.")

    return cube[keep], wavelengths[keep], keep


# ============================================================
# 5. Savitzky–Golay smoothing
# ============================================================

def smooth_cube(cube, window=11, poly=2):
    return savgol_filter(
        cube, axis=0, mode="mirror",
        window_length=window, polyorder=poly
    )


def destripe_columns(cube, frac=1.0, polyfit_order=None):
    """Lightweight destriping along columns for each band.

    - Computes per-band column medians ignoring NaNs and subtracts them.
    - Optionally fits a low-order polynomial across column medians to
      remove very low-frequency trends (set polyfit_order to 1 or 2).
    - Adds back the global median per band to preserve overall level.

    Args:
        cube: array (bands, rows, cols)
        frac: scale factor for subtraction (1.0 = full, 0.5 = half)
        polyfit_order: None or int for polynomial fit across columns.
    Returns:
        cube_out with reduced column striping.
    """
    if frac is None or frac <= 0:
        return cube
    bands, rows, cols = cube.shape
    out = cube.astype(np.float32, copy=True)
    x = np.arange(cols, dtype=np.float32)
    for b in range(bands):
        img = out[b]
        # Column medians ignoring NaNs
        col_med = np.nanmedian(img, axis=0)
        # Replace non-finite with overall median
        base_med = np.nanmedian(col_med) if np.isfinite(col_med).any() else 0.0
        col_med = np.where(np.isfinite(col_med), col_med, base_med)
        trend = col_med
        if polyfit_order is not None and polyfit_order >= 1:
            try:
                coeffs = np.polyfit(x, col_med, deg=int(polyfit_order))
                trend = np.polyval(coeffs, x)
            except Exception:
                trend = col_med
        # Scale subtraction and preserve level
        img_med = np.nanmedian(img) if np.isfinite(img).any() else 0.0
        correction = (frac * trend)[None, :]
        out[b] = img - correction + img_med
    return out


# ============================================================
# 6. MNF transform
# ============================================================

def _impute_nan_columns(X):
    """Impute NaNs and infs column-wise with column means; if a column is all NaN/inf, fill with 0.
    Args:
        X: 2D array (n_samples, n_features)
    Returns:
        X_out: imputed float32 array
    """
    X = X.astype(np.float64, copy=True)
    # Treat inf as NaN
    non_finite = ~np.isfinite(X)
    if np.any(non_finite):
        X[non_finite] = np.nan
    col_means = np.nanmean(X, axis=0)
    # Columns that are all NaN -> set mean to 0
    col_means = np.where(np.isfinite(col_means), col_means, 0.0)
    # Find NaNs to replace
    inds = np.where(np.isnan(X))
    if inds[0].size > 0:
        X[inds] = np.take(col_means, inds[1])
    return X.astype(np.float32, copy=False)


def estimate_noise(cube):
    diffs = cube[:, 1:, :] - cube[:, :-1, :]
    diffs = diffs.reshape(cube.shape[0], -1).T  # (n_samples, n_features)
    # Remove any samples (rows) containing NaNs/Infs
    mask = np.all(np.isfinite(diffs), axis=1)
    valid = diffs[mask]
    if valid.shape[0] == 0:
        raise ValueError("[EnMap] estimate_noise: No valid (finite) samples available for noise covariance estimation.")
    cov = EmpiricalCovariance().fit(valid)
    return cov.covariance_

def _whitening_matrix_from_cov(noise_cov, eps=0.0, rtol=1e-8):
    """Compute symmetric whitening matrix W from noise covariance.
    Given noise covariance N (bands x bands), compute eigendecomposition N = E Λ E^T
    and return W = E Λ^{-1/2} E^T suitable for right-multiplication X @ W
    where X is (n_samples x bands).
    Small/negative eigenvalues are clamped using max(λ, eps, rtol*λ_max).
    """
    # Ensure symmetric
    N = 0.5 * (noise_cov + noise_cov.T)
    # Eigen decomposition for symmetric matrices
    vals, vecs = np.linalg.eigh(N)
    lam_max = float(np.max(vals)) if vals.size > 0 else 0.0
    floor = max(eps, rtol * lam_max)
    vals_clamped = np.clip(vals, floor, None)
    inv_sqrt = 1.0 / np.sqrt(vals_clamped)
    # Recompose whitening matrix: E * Λ^{-1/2} * E^T
    W = (vecs * inv_sqrt) @ vecs.T  # vecs @ diag(inv_sqrt) @ vecs.T using broadcasting
    return W.astype(np.float64, copy=False)


def run_mnf(cube, n_components=12):
    bands, h, w = cube.shape
    X = cube.reshape(bands, -1).T  # (n_pixels, bands)

    # Estimate noise covariance robustly (handles NaNs internally)
    noise_cov = estimate_noise(cube)

    # Build proper eigen-based whitening matrix
    W = _whitening_matrix_from_cov(noise_cov, eps=0.0, rtol=1e-8)

    # Noise-whiten features using X_white = X @ W
    Xw = X @ W

    # Impute NaNs/infs before PCA
    Xw = _impute_nan_columns(Xw)

    # Ensure n_components <= bands
    k = min(n_components, bands)
    pca = PCA(n_components=k)
    mnf = pca.fit_transform(Xw)

    return mnf.T.reshape(k, h, w)


# ============================================================
# 7. Absorption depths
# ============================================================

def continuum_removed_depth(cube, wavelengths, center_range):
    start, end = center_range
    idx = np.where((wavelengths >= start) & (wavelengths <= end))[0]
    if len(idx) < 3:
        return None

    sub = cube[idx].astype(np.float32, copy=False)
    left = sub[0]
    right = sub[-1]
    continuum = np.linspace(left, right, sub.shape[0]).astype(np.float32, copy=False)
    # Avoid division by near-zero continuum values
    eps = 1e-6
    continuum = np.where(np.isfinite(continuum) & (np.abs(continuum) > eps), continuum, np.nan)
    cr = 1.0 - (sub / continuum)
    # Clamp to finite range
    cr = np.where(np.isfinite(cr), cr, np.nan)
    return np.nanmax(cr, axis=0)


# ============================================================
# 8. Derivative PCA
# ============================================================

def derivative_pca(cube, wavelengths, n_components=3):
    deriv = np.gradient(cube, axis=0)
    bands, h, w = deriv.shape

    X = deriv.reshape(bands, -1).T  # (n_pixels, bands)
    # Impute NaNs/infs before PCA
    X = _impute_nan_columns(X)

    k = min(n_components, bands)
    pca = PCA(n_components=k)
    pcs = pca.fit_transform(X)

    return pcs.T.reshape(k, h, w)


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



def enmap_to_feature_stack(folder_path, n_mnf=12, n_deriv_pcs=3, destripe=True, destripe_frac=1.0, destripe_poly=1):
    cube, meta, qa_masks, xml = load_enmap_dataset(folder_path)

    wavelengths = parse_enmap_wavelengths(xml)
    # Enforce strict consistency with metadata (no assumptions/harmonization)
    if wavelengths.size != cube.shape[0]:
        raise ValueError(
            f"[EnMap] Wavelength-band mismatch: wavelengths={wavelengths.size}, "
            f"cube_bands={cube.shape[0]}. The parser only returns bands described in the metadata XML."
        )

    mask = build_mask(qa_masks)

    cube, wavelengths, _ = drop_bad_bands(cube, wavelengths)

    cube = smooth_cube(cube)
    if destripe:
        cube = destripe_columns(cube, frac=destripe_frac, polyfit_order=destripe_poly)
        print(f"[EnMap] Applied column destriping: frac={destripe_frac}, poly_order={destripe_poly}")
    if mask is not None:
        cube[:, mask] = np.nan

    mnf = run_mnf(cube, n_components=n_mnf)

    depth_2200 = continuum_removed_depth(cube, wavelengths, (2200, 2230))
    depth_2300 = continuum_removed_depth(cube, wavelengths, (2300, 2330))
    depth_2340 = continuum_removed_depth(cube, wavelengths, (2320, 2350))
    depth_1000 = continuum_removed_depth(cube, wavelengths, (900, 1030))
    depth_1750 = continuum_removed_depth(cube, wavelengths, (1700, 1800))

    deriv = derivative_pca(cube, wavelengths, n_components=n_deriv_pcs)

    features = {}

    # MNF components (return raw arrays; caller wraps into rasters)
    for i in range(mnf.shape[0]):
        features[f"MNF_{i+1:02d}"] = mnf[i]

    # Absorption depths
    if depth_2200 is not None: features["Depth_2200"] = depth_2200
    if depth_2300 is not None: features["Depth_2300"] = depth_2300
    if depth_2340 is not None: features["Depth_2340"] = depth_2340
    if depth_1000 is not None: features["Depth_1000"] = depth_1000
    if depth_1750 is not None: features["Depth_1750"] = depth_1750

    # Derivative PCA
    for i in range(n_deriv_pcs):
        features[f"Deriv_PC{i+1}"] = deriv[i]

    return features, meta