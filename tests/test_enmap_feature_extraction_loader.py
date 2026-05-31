import os

import pytest

from mineye.BayesianSegmentation.EnMap_feature_extraction import (
    _find_required_enmap_file,
    _resolve_enmap_product_path,
)


def test_resolve_enmap_product_path_uses_nested_product_directory(tmp_path):
    product_name = "ENMAP01-____L2A-DT0000026661_20230712T114038Z_001_V010402_20240818T134118Z"
    product_dir = tmp_path / product_name
    product_dir.mkdir()
    (tmp_path / "__MACOSX").mkdir()
    (product_dir / f"{product_name}-SPECTRAL_IMAGE.TIF").touch()
    (product_dir / f"{product_name}-METADATA.XML").touch()

    assert _resolve_enmap_product_path(tmp_path) == os.fspath(product_dir)


def test_find_required_enmap_file_reports_available_files(tmp_path):
    files = ["README.txt", "QL_VNIR.TIF"]
    for name in files:
        (tmp_path / name).touch()

    with pytest.raises(FileNotFoundError, match="Could not find required EnMAP file") as exc_info:
        _find_required_enmap_file(tmp_path, files, "SPECTRAL_IMAGE", (".TIF", ".TIFF"))

    message = str(exc_info.value)
    assert "README.txt" in message
    assert "QL_VNIR.TIF" in message