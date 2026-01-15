import os
import numpy as np
import time
import warnings
try:
    import rasterio
    from rasterio.io import MemoryFile
except Exception:  # pragma: no cover
    rasterio = None
    MemoryFile = None
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter, median_filter, zoom
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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
    if rasterio is None:
        raise ImportError("rasterio is required for load_enmap_dataset(). Please install rasterio.")
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


def _masked_gaussian_lowpass(img2d: np.ndarray, valid2d: np.ndarray, sigma: float) -> np.ndarray:
    """Mask-aware Gaussian low-pass.

    Computes:
        B = Gσ(I * V) / Gσ(V)
    where V is 1 for valid pixels and 0 otherwise.
    """
    img = np.asarray(img2d, dtype=np.float32)
    v = np.asarray(valid2d, dtype=np.float32)
    if img.shape != v.shape:
        raise ValueError(f"valid2d must match img2d shape {img.shape}, got {v.shape}")
    if sigma is None or float(sigma) <= 0:
        return img.astype(np.float32, copy=False)

    num = gaussian_filter(img * v, sigma=float(sigma), mode="nearest")
    den = gaussian_filter(v, sigma=float(sigma), mode="nearest")

    # Avoid division by ~0 in fully-masked regions.
    eps = 1e-6
    b = num / np.maximum(den, eps)
    return b.astype(np.float32, copy=False)


def _masked_gaussian_lowpass_downsampled(
        img2d: np.ndarray,
        valid2d: np.ndarray,
        sigma: float,
        downsample: int,
) -> np.ndarray:
    """Downsample→Gaussian→upsample mask-aware low-pass for speed.

    This is an approximation of `_masked_gaussian_lowpass`.
    """
    ds = int(downsample)
    if ds <= 1:
        return _masked_gaussian_lowpass(img2d, valid2d, sigma)

    img = np.asarray(img2d, dtype=np.float32)
    v = np.asarray(valid2d, dtype=np.float32)
    h, w = img.shape

    # Block-mean downsample (pad to multiples of ds)
    pad_h = (-h) % ds
    pad_w = (-w) % ds
    if pad_h or pad_w:
        img_pad = np.pad(img, ((0, pad_h), (0, pad_w)), mode="edge")
        v_pad = np.pad(v, ((0, pad_h), (0, pad_w)), mode="edge")
    else:
        img_pad = img
        v_pad = v

    hp, wp = img_pad.shape
    h2, w2 = hp // ds, wp // ds
    img_small = img_pad.reshape(h2, ds, w2, ds).mean(axis=(1, 3)).astype(np.float32, copy=False)
    v_small = v_pad.reshape(h2, ds, w2, ds).mean(axis=(1, 3)).astype(np.float32, copy=False)

    # Scale sigma into low-res space.
    sigma_small = float(sigma) / float(ds)
    sigma_small = max(sigma_small, 0.5)  # avoid too-tiny kernels

    b_small = _masked_gaussian_lowpass(img_small, v_small > 0.0, sigma_small)
    b_pad = zoom(b_small, zoom=(ds, ds), order=1, mode="nearest", prefilter=False)
    b = b_pad[:h, :w]
    return b.astype(np.float32, copy=False)


def remove_background_field(
        cube: np.ndarray,
        mask_bad: np.ndarray | None,
        sigma: float = 75.0,
        downsample: int | None = None,
        progress_every: int | None = 1,
) -> np.ndarray:
    """Remove a smooth 2D background field per band (detrending).

    Implements additive correction:
        I_corr = I - B
    where B is a very smooth, mask-aware low-frequency field.

    Args:
        cube: array (bands, rows, cols)
        mask_bad: optional bool (rows, cols) where True pixels are invalid.
        sigma: Gaussian sigma in pixels (e.g. 50–100 for EnMAP 30 m).
        downsample: optional int for downsampled low-pass computation (e.g. 2 or 4).
        progress_every: optional int. If set, prints a progress line every N bands.
    Returns:
        cube_out float32 with background removed.
    """
    cube = np.asarray(cube)
    if cube.ndim != 3:
        raise ValueError(f"cube must be 3D (bands, rows, cols), got shape {cube.shape}")
    bands, rows, cols = cube.shape
    out = cube.astype(np.float32, copy=True)

    mb = None
    if mask_bad is not None:
        if mask_bad.shape != (rows, cols):
            raise ValueError(
                f"[EnMap] remove_background_field: mask_bad must have shape (rows, cols)={(rows, cols)}, got {mask_bad.shape}"
            )
        mb = mask_bad.astype(bool, copy=False)

    ds = int(downsample) if downsample is not None else None
    if ds is not None and ds <= 1:
        ds = None

    pe = None
    if progress_every is not None:
        pe = int(progress_every)
        if pe <= 0:
            pe = None

    t0 = time.perf_counter()
    print(
        f"[EnMap] Detrending/background removal: bands={bands}, rows={rows}, cols={cols}, "
        f"sigma={sigma}, downsample={ds}, mask_bad={'yes' if mb is not None else 'no'}"
    )

    for b in range(bands):
        tb = time.perf_counter()
        img = out[b]

        finite = np.isfinite(img)
        if mb is not None:
            finite = finite & (~mb)

        valid_frac = float(np.mean(finite)) if img.size > 0 else 0.0

        # Fill non-finite with 0 for the numerator; validity mask handles normalization.
        filled = np.where(np.isfinite(img), img, 0.0).astype(np.float32, copy=False)

        if ds is not None:
            bg = _masked_gaussian_lowpass_downsampled(filled, finite, sigma=float(sigma), downsample=ds)
        else:
            bg = _masked_gaussian_lowpass(filled, finite, sigma=float(sigma))

        corr = img - bg
        out[b] = corr.astype(np.float32, copy=False)

        # Lightweight diagnostics
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered")
            bg_std = float(np.nanstd(bg))
            corr_std = float(np.nanstd(corr))

        if pe is not None and ((b + 1) % pe == 0 or (b + 1) == bands):
            dt_band = time.perf_counter() - tb
            dt_total = time.perf_counter() - t0
            print(
                f"[EnMap] Detrending progress: band {b + 1}/{bands} "
                f"(dt={dt_band:.2f}s, total={dt_total:.1f}s, valid={valid_frac * 100.0:.1f}%, "
                f"background_std={bg_std:.4g}, corrected_std={corr_std:.4g})"
            )

    return out


def _lowfreq_reference_median_filter(img, kernel):
    """Compute a low-frequency reference image via 2D median filter.

    Notes
    -----
    - NaNs/Infs are filled with the global median before filtering.
    - This function is intentionally lightweight and self-contained.

    Args:
        img: 2D array (rows, cols)
        kernel: int, median filter window size (odd recommended)
    Returns:
        ref: 2D float32 array
    """
    if kernel is None or int(kernel) <= 1:
        return img.astype(np.float32, copy=False)

    k = int(kernel)
    if k % 2 == 0:
        k += 1

    img = img.astype(np.float32, copy=False)
    finite = np.isfinite(img)
    fill = np.nanmedian(img[finite]) if np.any(finite) else 0.0
    filled = np.where(finite, img, fill)
    ref = median_filter(filled, size=(k, k), mode="nearest")
    return ref.astype(np.float32, copy=False)


def _lowfreq_reference_median_filter_downsampled(img, kernel, downsample):
    """Compute a low-frequency reference via downsample→median-filter→upsample.

    This is a speed optimization for residual-based destriping.

    Args:
        img: 2D array (rows, cols)
        kernel: int, median filter window size (odd recommended)
        downsample: int, factor (>=2) used for block-mean downsampling
    Returns:
        ref: 2D float32 array with the same shape as img
    """
    if kernel is None or int(kernel) <= 1:
        return img.astype(np.float32, copy=False)

    ds = int(downsample)
    if ds <= 1:
        return _lowfreq_reference_median_filter(img, kernel)

    k = int(kernel)
    if k % 2 == 0:
        k += 1

    img = img.astype(np.float32, copy=False)
    h, w = img.shape

    finite = np.isfinite(img)
    fill = np.nanmedian(img[finite]) if np.any(finite) else 0.0
    filled = np.where(finite, img, fill)

    # Pad to multiples of ds for block aggregation
    pad_h = (-h) % ds
    pad_w = (-w) % ds
    if pad_h or pad_w:
        filled_pad = np.pad(filled, ((0, pad_h), (0, pad_w)), mode="edge")
    else:
        filled_pad = filled

    hp, wp = filled_pad.shape
    h2, w2 = hp // ds, wp // ds
    small = filled_pad.reshape(h2, ds, w2, ds).mean(axis=(1, 3)).astype(np.float32, copy=False)

    # Apply median filter in low-res space
    ref_small = median_filter(small, size=(k, k), mode="nearest").astype(np.float32, copy=False)

    # Upsample back; use bilinear interpolation for a smooth low-frequency reference
    ref_pad = zoom(ref_small, zoom=(ds, ds), order=1, mode="nearest", prefilter=False)
    ref = ref_pad[:h, :w]
    return ref.astype(np.float32, copy=False)


def destripe_columns(
        cube,
        frac=1.0,
        polyfit_order=None,
        mask_bad=None,
        reference_kernel=None,
        reference_downsample=None,
        smooth_cols=None,
        progress_every: int | None = 1,
):
    """Lightweight destriping along columns for each band.

    - Computes per-band column medians ignoring NaNs and subtracts them.
    - Optionally fits a low-order polynomial across column medians to
      remove very low-frequency trends (set polyfit_order to 1 or 2).
    - Adds back the global median per band to preserve overall level.

    Optional extensions (backwards-compatible):
    - mask_bad: 2D boolean mask (rows, cols) where True pixels are excluded
      from stripe estimation (internally treated as NaN for estimating column
      medians). The correction is still applied to the full image.
    - reference_kernel: if set (e.g., 31 or 51), estimate the column bias on a
      high-pass residual E = I - R, where R is a low-frequency median-filter
      reference. This reduces the risk of removing real broad gradients.
    - reference_downsample: if set (e.g., 2 or 4), computes the low-frequency
      reference on a downsampled image and upsamples it back. This can drastically
      reduce runtime for large images/kernels at the cost of a small approximation.
    - smooth_cols: if set (e.g., 11 or 21), apply a 1D median filter to the
      estimated column bias vector across columns.

    Args:
        cube: array (bands, rows, cols)
        frac: scale factor for subtraction (1.0 = full, 0.5 = half)
        polyfit_order: None or int for polynomial fit across columns.
        mask_bad: optional bool array (rows, cols), True = exclude from estimation
        reference_kernel: optional int, median-filter window for low-frequency reference
        reference_downsample: optional int, factor for downsample→filter→upsample reference computation
        smooth_cols: optional int, median-filter window across columns
        progress_every: optional int. If set, prints a progress line every N bands.
            Default is 1 (print after each band) to give interactive feedback on long runs.
    Returns:
        cube_out with reduced column striping.
    """
    if frac is None or frac <= 0:
        return cube
    bands, rows, cols = cube.shape
    out = cube.astype(np.float32, copy=True)
    x = np.arange(cols, dtype=np.float32)

    # NOTE: `reference_kernel` triggers a 2D median filter per band, which can be very
    # expensive for large images / large kernels.
    t0 = time.perf_counter()
    rk = int(reference_kernel) if reference_kernel is not None else None
    rds = int(reference_downsample) if reference_downsample is not None else None
    if rds is not None and rds <= 1:
        rds = None
    if rk is not None and rk > 1:
        k_eff = rk + (1 if rk % 2 == 0 else 0)
        px = rows * cols
        print(
            f"[EnMap] Destriping {bands} bands (rows={rows}, cols={cols}), "
            f"mode=residual(reference_kernel={k_eff}, downsample={rds}), smooth_cols={smooth_cols}, poly_order={polyfit_order}, "
            f"frac={frac}, mask_bad={'yes' if mask_bad is not None else 'no'}. "
            f"This can take a long time for large kernels (median filter per band over {px:,} pixels)."
        )
    else:
        print(
            f"[EnMap] Destriping {bands} bands (rows={rows}, cols={cols}), "
            f"mode=direct, smooth_cols={smooth_cols}, poly_order={polyfit_order}, "
            f"frac={frac}, mask_bad={'yes' if mask_bad is not None else 'no'}."
        )

    if mask_bad is not None:
        if mask_bad.shape != (rows, cols):
            raise ValueError(
                f"[EnMap] destripe_columns: mask_bad must have shape (rows, cols)={(rows, cols)}, got {mask_bad.shape}"
            )
        mask_bad = mask_bad.astype(bool, copy=False)

    smooth_k = None
    if smooth_cols is not None:
        smooth_k = int(smooth_cols)
        if smooth_k <= 1:
            smooth_k = None
        elif smooth_k % 2 == 0:
            smooth_k += 1

    pe = None
    if progress_every is not None:
        pe = int(progress_every)
        if pe <= 0:
            pe = None

    for b in range(bands):
        tb = time.perf_counter()
        img = out[b]
        est_img = img
        if mask_bad is not None:
            est_img = img.copy()
            est_img[mask_bad] = np.nan

        if reference_kernel is not None and int(reference_kernel) > 1:
            if rds is not None:
                ref = _lowfreq_reference_median_filter_downsampled(est_img, reference_kernel, rds)
            else:
                ref = _lowfreq_reference_median_filter(est_img, reference_kernel)
            est = est_img - ref
        else:
            est = est_img

        # Lightweight per-band stats (computed on the estimation image) for logging.
        finite_est = np.isfinite(est)
        valid_frac = float(np.mean(finite_est)) if est.size > 0 else 0.0
        valid_cols = int(np.sum(np.any(finite_est, axis=0))) if cols > 0 else 0

        # Column medians ignoring NaNs.
        # Some columns may be fully masked -> expected "All-NaN slice" warnings.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered")
            col_med = np.nanmedian(est, axis=0)
        # Replace non-finite with overall median
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered")
            base_med = np.nanmedian(col_med) if np.isfinite(col_med).any() else 0.0
        col_med = np.where(np.isfinite(col_med), col_med, base_med)

        trend = col_med
        if smooth_k is not None:
            trend = median_filter(trend.astype(np.float32, copy=False), size=smooth_k, mode="nearest")

        if polyfit_order is not None and polyfit_order >= 1:
            try:
                coeffs = np.polyfit(x, trend, deg=int(polyfit_order))
                trend = np.polyval(coeffs, x)
            except Exception:
                # Fall back to the non-polyfit trend
                pass

        # Informative per-band destriping metrics.
        # - correction_std: magnitude of applied per-column correction (offset units)
        # - stripe_score: high-frequency column variation remaining in the estimated column medians
        #   after removing the (smoothed/polyfit) trend. Lower is better.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered")
            correction_std = float(np.nanstd(trend))
            stripe_score = float(np.nanstd(col_med - trend))

        correction = (frac * trend)[None, :]
        if reference_kernel is not None and int(reference_kernel) > 1:
            # Residual-based correction: apply directly (do not add band median).
            out[b] = img - correction
        else:
            # Original behavior: scale subtraction and preserve level.
            # If mask_bad is provided, preserve the level relative to the valid pixels
            # used for estimation (i.e., ignore masked-out outliers when computing the band median).
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="All-NaN slice encountered")
                img_med = np.nanmedian(est_img) if np.isfinite(est_img).any() else 0.0
            out[b] = img - correction + img_med

        if pe is not None and ((b + 1) % pe == 0 or (b + 1) == bands):
            dt_band = time.perf_counter() - tb
            dt_total = time.perf_counter() - t0
            print(
                f"[EnMap] Destriping progress: band {b + 1}/{bands} "
                f"(dt={dt_band:.2f}s, total={dt_total:.1f}s, "
                f"valid={valid_frac * 100.0:.1f}%, valid_cols={valid_cols}/{cols}, "
                f"correction_std={correction_std:.4g}, stripe_score={stripe_score:.4g})"
            )
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
    # Lazy import to avoid hard dependency during module import in minimal environments
    try:
        from sklearn.covariance import EmpiricalCovariance
    except Exception as e:
        raise ImportError("scikit-learn is required for estimate_noise(). Please install scikit-learn.") from e
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

    # Mask-aware PCA fitting: fit only on fully finite rows
    finite_rows = np.all(np.isfinite(Xw), axis=1)
    Xw_valid = Xw[finite_rows]

    # Handle edge cases: if no valid samples, raise; if too few, adjust components
    if Xw_valid.shape[0] == 0:
        raise ValueError("[EnMap] run_mnf: No valid samples to fit PCA (after whitening).")

    k = int(min(n_components, bands, Xw_valid.shape[1], Xw_valid.shape[0]))
    if k < 1:
        raise ValueError("[EnMap] run_mnf: n_components invalid after masking.")

    # Lazy import to avoid hard dependency at module import time
    try:
        from sklearn.decomposition import PCA
    except Exception as e:
        raise ImportError("scikit-learn is required for run_mnf(). Please install scikit-learn.") from e
    pca = PCA(n_components=k)
    pca.fit(Xw_valid)

    # For transforming all pixels, impute only non-finite entries using means from valid set
    col_means = np.nanmean(Xw_valid, axis=0)
    col_means = np.where(np.isfinite(col_means), col_means, 0.0)
    Xw_for_transform = Xw.copy()
    non_finite = ~np.isfinite(Xw_for_transform)
    if np.any(non_finite):
        # replace each non-finite entry with corresponding column mean
        rows, cols = np.where(non_finite)
        if rows.size > 0:
            Xw_for_transform[rows, cols] = col_means[cols]

    mnf = pca.transform(Xw_for_transform)

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
    # Safe column-wise max without RuntimeWarning for all-NaN slices
    cr_finite = np.where(np.isfinite(cr), cr, -np.inf)
    m = np.max(cr_finite, axis=0)
    m[~np.isfinite(m) | (m == -np.inf)] = np.nan
    return m


# ============================================================
# 8. Derivative PCA
# ============================================================

def _select_swir_window(wavelengths, preferred=(2000.0, 2400.0), fallback=(1900.0, 2450.0), min_bands=3):
    """Return indices for an automatically selected SWIR window based on wavelengths.

    Strategy:
    - Prefer [2000, 2400] nm if it contains >= min_bands.
    - Else try a broader window [1900, 2450] nm if it contains >= min_bands.
    - Else return None to indicate 'use all bands'.
    Returns: (idx, chosen_range or None)
    """
    wl = np.asarray(wavelengths, dtype=float)
    if wl.ndim != 1 or wl.size == 0:
        return None, None
    if not np.isfinite(wl).any():
        return None, None
    # Clean inf/nan to avoid comparisons failing
    finite_mask = np.isfinite(wl)
    wl_clean = wl.copy()
    wl_clean[~finite_mask] = np.nan

    def pick(window):
        lo, hi = float(window[0]), float(window[1])
        idx = np.where((wl_clean >= lo) & (wl_clean <= hi))[0]
        return idx

    for window in (preferred, fallback):
        idx = pick(window)
        if idx.size >= int(min_bands):
            return idx, window
    return None, None


def derivative_pca(cube, wavelengths, n_components=3, swir_range="auto"):
    """Compute derivative PCA restricted to a SWIR absorption region to reduce noise.

    Args:
        cube: array (bands, h, w)
        wavelengths: array (bands,) in nm
        n_components: number of PCs to return
        swir_range: one of
            - "auto" (default): automatically select a SWIR window from wavelengths.
            - None: use all bands.
            - (min_nm, max_nm) tuple: explicit window.
        If fewer than 3 bands fall in the selected range, falls back to using all bands.
    """
    wl = np.asarray(wavelengths).astype(float)

    use_all = False
    subcube = cube

    if swir_range == "auto":
        idx, chosen = _select_swir_window(wl)
        if idx is not None and idx.size >= 3:
            subcube = cube[idx]
            print(f"[EnMap] derivative_pca(auto): using SWIR window {chosen} with {idx.size} bands.")
        else:
            use_all = True
            print("[EnMap] derivative_pca(auto): no suitable SWIR window found (>=3 bands). Using all bands.")
    elif swir_range is None:
        use_all = True
    else:
        # Explicit tuple provided
        lo, hi = swir_range
        idx = np.where((wl >= float(lo)) & (wl <= float(hi)))[0]
        if idx.size >= 3:
            subcube = cube[idx]
        else:
            use_all = True
            print(f"[EnMap] derivative_pca: SWIR range {swir_range} yields {idx.size} bands (<3). Falling back to all bands.")

    if use_all:
        subcube = cube

    deriv = np.gradient(subcube, axis=0)
    bands, h, w = deriv.shape

    X = deriv.reshape(bands, -1).T  # (n_pixels, bands)

    # Fit PCA on valid (fully finite) rows only
    finite_rows = np.all(np.isfinite(X), axis=1)
    X_valid = X[finite_rows]
    if X_valid.shape[0] == 0:
        raise ValueError("[EnMap] derivative_pca: No valid samples to fit PCA (after derivative).")

    k = int(min(n_components, bands, X_valid.shape[1], X_valid.shape[0]))
    if k < 1:
        raise ValueError("[EnMap] derivative_pca: n_components invalid after masking.")

    # Lazy import to avoid hard dependency at module import time
    try:
        from sklearn.decomposition import PCA
    except Exception as e:
        raise ImportError("scikit-learn is required for derivative_pca(). Please install scikit-learn.") from e
    pca = PCA(n_components=k)
    pca.fit(X_valid)

    # Prepare full transform: impute non-finite with column means from valid set
    col_means = np.nanmean(X_valid, axis=0)
    col_means = np.where(np.isfinite(col_means), col_means, 0.0)
    X_for_transform = X.copy()
    non_finite = ~np.isfinite(X_for_transform)
    if np.any(non_finite):
        rows, cols = np.where(non_finite)
        if rows.size > 0:
            X_for_transform[rows, cols] = col_means[cols]

    pcs = pca.transform(X_for_transform)

    return pcs.T.reshape(k, h, w)


# ============================================================
# 9. Wrap 2D array into an in-memory rasterio dataset
# ============================================================

def validate_features_dict(features, meta):
    """Validate that features is a dict[str, rasterio DatasetReader] matching meta.

    Requirements:
    - features is a dict with string keys
    - each value is a rasterio DatasetReader with:
      - count == 1
      - dtype == float32
      - CRS and transform present
      - height/width match meta
    Returns True if valid; otherwise raises ValueError with an informative message.
    """
    if not isinstance(features, dict):
        raise ValueError("features must be a dict of name -> rasterio dataset")
    req_h = meta.get("height")
    req_w = meta.get("width")
    if req_h is None or req_w is None:
        raise ValueError("meta must include 'height' and 'width'")
    for k, v in features.items():
        if not isinstance(k, str):
            raise ValueError(f"Feature key must be str, got {type(k)} for key {k!r}")
        if not (hasattr(v, "read") and hasattr(v, "profile")):
            raise ValueError(f"Feature '{k}' must be a rasterio dataset, got {type(v)}")
        try:
            h, w = int(v.height), int(v.width)
            count = int(v.count)
            dtypes = list(v.dtypes)
            crs = v.crs
            transform = v.transform
        except Exception as e:
            raise ValueError(f"Feature '{k}' dataset is not readable: {e}")
        if h != req_h or w != req_w:
            raise ValueError(
                f"Feature '{k}' has size ({h},{w}) but expected ({req_h},{req_w}) from meta"
            )
        if count != 1:
            raise ValueError(f"Feature '{k}' must be single-band (count=1), got count={count}")
        if not dtypes or dtypes[0] != "float32":
            got = dtypes[0] if dtypes else "unknown"
            raise ValueError(f"Feature '{k}' must have dtype float32, got {got}")
        if crs is None:
            raise ValueError(f"Feature '{k}' is missing CRS")
        if transform is None:
            raise ValueError(f"Feature '{k}' is missing transform")
    return True


def create_rasterio_layer(array, meta):
    """Return a rasterio in-memory DatasetReader for a numpy array.

    Important: Keeps the underlying MemoryFile alive by attaching it to
    the returned dataset object (dataset._memfile), so callers can use the
    dataset safely without immediately managing the MemoryFile lifetime.

    Note on masking values:
    - Internally, model fitting uses NaNs to exclude invalid pixels.
    - For outputs consumed by downstream Bayesian segmentation, we replace
      NaNs/Infs with a sentinel (-9999.0) because BaySeg cannot handle NaNs.
    """
    if MemoryFile is None:
        raise ImportError("rasterio is required for create_rasterio_layer(). Please install rasterio.")
    # Ensure float32 and replace non-finite values with sentinel
    arr = np.asarray(array).astype("float32", copy=False)
    # Replace NaN/Inf with sentinel to make downstream consumers robust
    sentinel = -9999.0
    arr = np.nan_to_num(arr, nan=sentinel, posinf=sentinel, neginf=sentinel)

    new_meta = meta.copy()
    new_meta.update({
        "count": 1,
        "dtype": "float32",
        "height": arr.shape[-2],
        "width": arr.shape[-1],
        "nodata": sentinel,
    })

    memfile = MemoryFile()
    # Write then reopen read-only; attach memfile to keep it alive
    with memfile.open(**new_meta) as dst:
        dst.write(arr, 1)
    ds = memfile.open()  # rasterio.io.DatasetReader
    # Attach the memfile so it doesn't get garbage-collected
    setattr(ds, "_memfile", memfile)
    return ds


# ============================================================
# 10. Main pipeline
# ============================================================



def enmap_to_feature_stack(
        folder_path,
        n_mnf=12,
        n_deriv_pcs=3,
        detrend=False,
        detrend_sigma=75.0,
        detrend_downsample=None,
        destripe=True,
        destripe_frac=1.0,
        destripe_poly=1,
        destripe_reference_kernel=None,
        destripe_reference_downsample=None,
        destripe_smooth_cols=None,
):
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

    if detrend:
        print(
            f"[EnMap] Detrending config: enabled={detrend}, sigma={detrend_sigma}, downsample={detrend_downsample}, "
            f"mask_bad={'yes' if mask is not None else 'no'}"
        )
        cube = remove_background_field(
            cube,
            mask_bad=mask,
            sigma=float(detrend_sigma),
            downsample=detrend_downsample,
            progress_every=1,
        )

    # Savitzky–Golay smoothing (spectral). Must run before injecting NaNs.
    cube = smooth_cube(cube)

    # Apply QA mask for downstream steps that are NaN-aware (MNF, depths, derivative PCA).
    if mask is not None:
        cube[:, mask] = np.nan

    if destripe:
        print(
            f"[EnMap] Destriping config: enabled={destripe}, frac={destripe_frac}, poly_order={destripe_poly}, "
            f"reference_kernel={destripe_reference_kernel}, reference_downsample={destripe_reference_downsample}, "
            f"smooth_cols={destripe_smooth_cols}, mask_bad={'yes' if mask is not None else 'no'}"
        )
        cube = destripe_columns(
            cube,
            frac=destripe_frac,
            polyfit_order=destripe_poly,
            mask_bad=mask,
            reference_kernel=destripe_reference_kernel,
            reference_downsample=destripe_reference_downsample,
            smooth_cols=destripe_smooth_cols,
        )
        print(
            f"[EnMap] Applied column destriping: frac={destripe_frac}, poly_order={destripe_poly}, "
            f"reference_kernel={destripe_reference_kernel}, reference_downsample={destripe_reference_downsample}, "
            f"smooth_cols={destripe_smooth_cols}"
        )

    # Keep as a safeguard in case future pipeline edits insert operations that reintroduce non-NaNs.
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

    # MNF components (wrapped into rasterio MemoryFile datasets)
    for i in range(mnf.shape[0]):
        features[f"MNF_{i+1:02d}"] = create_rasterio_layer(mnf[i], meta)

    # Absorption depths (wrapped if available)
    if depth_2200 is not None:
        features["Depth_2200"] = create_rasterio_layer(depth_2200, meta)
    if depth_2300 is not None:
        features["Depth_2300"] = create_rasterio_layer(depth_2300, meta)
    if depth_2340 is not None:
        features["Depth_2340"] = create_rasterio_layer(depth_2340, meta)
    if depth_1000 is not None:
        features["Depth_1000"] = create_rasterio_layer(depth_1000, meta)
    if depth_1750 is not None:
        features["Depth_1750"] = create_rasterio_layer(depth_1750, meta)

    # Derivative PCA
    # Use the number of computed components to avoid indexing beyond available
    n_deriv_out = deriv.shape[0]
    for i in range(n_deriv_out):
        features[f"Deriv_PC{i+1}"] = create_rasterio_layer(deriv[i], meta)

    # Validate output conforms to expected dictionary schema (name -> rasterio dataset)
    try:
        validate_features_dict(features, meta)
        print(f"[EnMap] Feature dict validated: {len(features)} layers, size=({meta.get('height')},{meta.get('width')})")
    except Exception as e:
        raise ValueError(f"[EnMap] Feature dict validation failed: {e}")

    return features, meta


# ============================================================
# 11. Visualization utility: save per-band plots
# ============================================================

def save_enmap_band_plots(folder_path, apply_mask=True, pmin=2.0, pmax=98.0, vmin=None, vmax=None, dpi=150, figsize=(6, 5)):
    """Save PNG plots for all bands in an EnMAP dataset using viridis colorscale.

    Parameters
    ----------
    folder_path : str
        Path to the EnMAP product folder containing SPECTRAL_IMAGE*.tif and METADATA*.xml.
    apply_mask : bool
        If True, combine available QA layers and mask bad pixels (set to NaN) before plotting.
    pmin, pmax : float
        Percentiles (0-100) for per-band dynamic range if vmin/vmax are not provided.
    vmin, vmax : float or None
        Global display range for all bands. If both provided, overrides per-band percentiles.
    dpi : int
        Output image DPI.
    figsize : tuple
        Matplotlib figure size in inches.

    Output
    ------
    Creates a directory '<folder_path>/intermediate_Data/plots_<dataset_name>' and saves one PNG per band:
      band_###_wlXXXXnm.png (wavelength shown if available).
    """
    if rasterio is None:
        raise ImportError("rasterio is required for save_enmap_band_plots(). Please install rasterio.")
    # Load cube and metadata
    cube, meta, qa_masks, xml = load_enmap_dataset(folder_path)

    # Parse wavelengths (may be empty)
    try:
        wavelengths = parse_enmap_wavelengths(xml)
    except Exception as e:
        print(f"[EnMap] Wavelength parse failed: {e}")
        wavelengths = np.array([], dtype=float)

    # Build mask
    mask = build_mask(qa_masks) if apply_mask else None

    # Prepare output directory
    dataset_name = os.path.basename(os.path.abspath(folder_path)) or "dataset"
    out_dir = os.path.join(folder_path, "intermediate_Data", f"plots_{dataset_name}")
    os.makedirs(out_dir, exist_ok=True)
    print(f"[EnMap] Saving per-band plots to: {out_dir}")

    H = meta.get('height')
    W = meta.get('width')
    nbands = cube.shape[0]

    # Flatten mask for fast nan assignment
    if mask is not None:
        if mask.shape != (H, W):
            print(f"[EnMap] Warning: QA mask shape {mask.shape} != raster shape {(H, W)}; ignoring mask.")
            mask = None

    # Determine global limits if requested
    if vmin is not None and vmax is not None:
        global_vmin, global_vmax = float(vmin), float(vmax)
    else:
        global_vmin = global_vmax = None

    for i in range(nbands):
        img = cube[i].astype(np.float32)
        if mask is not None:
            img = img.copy()
            img[mask] = np.nan

        # Compute display limits
        if global_vmin is not None:
            lo, hi = global_vmin, global_vmax
        else:
            finite = np.isfinite(img)
            if not np.any(finite):
                print(f"[EnMap] Band {i+1}: all values NaN/inf, skipping percentile stretch; using defaults 0..1")
                lo, hi = 0.0, 1.0
            else:
                lo = float(np.nanpercentile(img[finite], pmin))
                hi = float(np.nanpercentile(img[finite], pmax))
                if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
                    lo, hi = 0.0, 1.0

        # Title and filename
        if wavelengths.size == nbands:
            wl = wavelengths[i]
            title = f"Band {i+1} — {wl:.1f} nm"
            fname = f"band_{i+1:03d}_wl{int(round(wl)):04d}nm.png"
        else:
            title = f"Band {i+1}"
            fname = f"band_{i+1:03d}.png"

        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
        im = ax.imshow(img, cmap='viridis', vmin=lo, vmax=hi)
        ax.set_title(title)
        ax.set_axis_off()
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Reflectance')
        out_path = os.path.join(out_dir, fname)
        fig.tight_layout()
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)

    print(f"[EnMap] Saved {nbands} band plot(s) to {out_dir}")