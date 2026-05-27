"""Local conftest for model2 inversion tests.

Provides the ``data_dir`` fixture used by the magnetic point
extraction tests without importing gempy (which requires the full
GemPy stack).
"""

import os

import dotenv
import pyro
import torch

dotenv.load_dotenv()
os.environ["VALIDATE_SERIALIZATION"] = ""

seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)

import pytest


@pytest.fixture(scope="session")
def base_dir():
    """Project root directory."""
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..', '..', '..'),
    )


@pytest.fixture(scope="session")
def data_dir(base_dir):
    """Directory containing general input data."""
    return os.path.join(
        base_dir, 'examples', 'Data', 'General_Input_Data',
    )