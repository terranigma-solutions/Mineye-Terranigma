#!/bin/bash

# Script to run Sphinx documentation build with the correct virtual environment
# Recommended to be used as a 'Shell Script' run configuration in PyCharm

# Target virtual environment
VENV_PATH="$HOME/.venv/2025"
DOCS_DIR="$(dirname "$0")/docs"

# 1. Check if the virtual environment exists
if [ ! -d "$VENV_PATH" ]; then
    echo "❌ Error: Virtual environment not found at $VENV_PATH"
    exit 1
fi

# 2. Activate the virtual environment
echo "🔄 Activating virtual environment: $VENV_PATH"
source "$VENV_PATH/bin/activate"

# 3. Navigate to the docs directory
if [ -d "$DOCS_DIR" ]; then
    cd "$DOCS_DIR"
else
    echo "❌ Error: Docs directory not found at $DOCS_DIR"
    exit 1
fi

# 4. Run make with provided arguments or default to 'html'
TARGET=${1:-html}
echo "🚀 Running Sphinx build: make $TARGET"
make "$TARGET"

# Deactivate is not strictly necessary in a script but good practice
deactivate
echo "✅ Done."
