import os
from os import path
from setuptools import setup, find_packages


def read_requirements(file_name):
    """Read requirements from a file."""
    requirements = []
    if not os.path.exists(file_name):
        return requirements

    with open(file_name, "r", encoding="utf-8") as f:
        for line in f:
            # Strip whitespace and ignore comments
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            requirements.append(line)

    return requirements


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='mineye-terranigma',
    packages=find_packages(exclude=('test', 'tests', 'docs', 'examples')),
    install_requires=read_requirements("requirements.txt"),
    url='https://github.com/leguark/Mineye-Terranigma',
    license='MIT',
    author='Terranigma Solutions GmbH',
    author_email='miguel@terranigma-solutions.com',
    description='Bayesian segmentation of satellite imagery and probabilistic geological modeling with GemPy',
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords=['geology', 'probabilistic modeling', 'bayesian inference', 'satellite imagery', 'gempy'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    setup_requires=['setuptools_scm'],
    use_scm_version={
        "root": ".",
        "relative_to": __file__,
        "write_to": path.join("mineye", "_version.py"),
        "fallback_version": "0.1.0"
    },
    python_requires='>=3.10',
)
