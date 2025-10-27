#!/bin/bash

# Profile test script for PyCharm with enhanced GPU debugging
set -e  # Exit on any error

# Configuration variables
PROFILE_TIMEOUT=60  # Timeout in seconds
SAMPLING_RATE=1009  # py-spy sampling rate in Hz
TEST_FILE="test_benchmark_I.py"
OUTPUT_PREFIX="profile_precission"
SPEEDSCOPE_URL="https://www.speedscope.app/"

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

# Navigate to project root directory
echo "Navigating to project directory..."
PROJECT_DIR="$HOME/PycharmProjects/Mineye-Terranigma"
cd "$PROJECT_DIR"

# Check if test file exists
TEST_PATH="tests/test_benchmarking/$TEST_FILE"
if [ ! -f "$TEST_PATH" ]; then
    echo "Error: $TEST_PATH not found"
    exit 1
fi

# Set PyKeOps environment variables for debugging
export PYKEOPS_BUILD_TYPE="Debug"
export CUDA_LAUNCH_BLOCKING=1  # Synchronous CUDA calls for better debugging

# Generate timestamp for unique filenames
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_FILE="${OUTPUT_PREFIX}_${TIMESTAMP}"

echo "Starting py-spy profiling with GPU debugging..."
echo "Output file: ${OUTPUT_FILE}"
echo "CUDA_LAUNCH_BLOCKING enabled for debugging"
echo "Press Ctrl+C to stop profiling early"
echo ""

# Run py-spy with shorter duration and more frequent sampling for GPU debugging
timeout ${PROFILE_TIMEOUT}s py-spy record -o "${OUTPUT_FILE}" -f speedscope -d ${PROFILE_TIMEOUT} -r ${SAMPLING_RATE} --idle -- python -m pytest ${TEST_PATH} -v -s

# Check if profile was created
if [ -f "${OUTPUT_FILE}" ]; then
    echo ""
    echo "=== Profiling completed successfully ==="
    echo "Profile saved to: ${OUTPUT_FILE}"
    echo "Visit ${SPEEDSCOPE_URL} to view the profile"
    echo "File location: $(pwd)/${OUTPUT_FILE}"
else
    echo ""
    echo "=== Warning: Profile file was not created ==="
    echo "Check the output above for errors"
fi

echo "=== Script finished ==="