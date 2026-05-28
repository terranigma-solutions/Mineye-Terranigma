#!/bin/bash

# Profile test_run_predictive for Model2 magnetic inversion
set -e  # Exit on any error

echo "=== Starting profiling script ==="

# Activate virtual environment
echo "Activating virtual environment..."
VENV_PATH="$HOME/.venv/2025/bin/activate"
if [ -f "$VENV_PATH" ]; then
    source "$VENV_PATH"
else
    echo "Virtual environment not found at $VENV_PATH"
    exit 1
fi

# Navigate to test directory
echo "Navigating to test directory..."
TEST_DIR="$HOME/Projects/gempy_project/Mineye-Terranigma/tests/tests_inversions/model2"
cd "$TEST_DIR"

# Check if test file exists
if [ ! -f "test_magnetics_invesion.py" ]; then
    echo "Error: test_magnetics_invesion.py not found in current directory"
    exit 1
fi

# Set PyKeOps environment variables for debugging
export PYKEOPS_BUILD_TYPE="Debug"
export CUDA_LAUNCH_BLOCKING=1  # Synchronous CUDA calls for better debugging

# Generate timestamp for unique filenames
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_FILE="profile_magnetics_predictive_${TIMESTAMP}"

echo "Starting py-spy profiling with GPU debugging..."
echo "Output file: ${OUTPUT_FILE}"
echo "CUDA_LAUNCH_BLOCKING enabled for debugging"
echo "Press Ctrl+C to stop profiling early"
echo ""

# Run py-spy profiling test_run_predictive
timeout 60s py-spy record -o "${OUTPUT_FILE}" -f speedscope -d 60 -r 10 --idle -- python -m pytest test_magnetics_invesion.py::TestMagneticInversion::test_run_predictive -v -s

# Check if profile was created
if [ -f "${OUTPUT_FILE}" ]; then
    echo ""
    echo "=== Profiling completed successfully ==="
    echo "Profile saved to: ${OUTPUT_FILE}"
    echo "Visit https://www.speedscope.app/ to view the profile"
    echo "File location: $(pwd)/${OUTPUT_FILE}"
else
    echo ""
    echo "=== Warning: Profile file was not created ==="
    echo "Check the output above for errors"
fi

echo "=== Script finished ==="