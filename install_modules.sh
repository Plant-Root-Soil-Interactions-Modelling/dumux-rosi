#!/usr/bin/env bash
set -euo pipefail




# 1. Get Python's architecture-specific site-packages directory
# Detect if we are in a virtual environment
IN_VENV=$(python3 - <<'EOF'
import sys
print(sys.prefix != sys.base_prefix)
EOF
)

if [[ "$IN_VENV" == "True" ]]; then
    # Install into the active virtual environment
    USER_SITE=$(python3 -c "import sysconfig; print(sysconfig.get_path('platlib'))")
    echo "Virtual environment detected."
else
    # Install into user's site-packages
    USER_SITE=$(python3 -c "import site; print(site.getusersitepackages())")
    echo "No virtual environment detected; installing into user site-packages."
fi


# 2. Define destination folder
DEST="$USER_SITE/rosi"

echo "Installing .so and .py modules into: $DEST"

# 3. Create destination directory
mkdir -p "$DEST"

# 4. Copy all .so and .py files
cp build-cmake/cpp/python_binding/*.so "$DEST/"
cp python/modules/*.py "$DEST/"
cp -r python/modules/fv "$DEST/fv"

# 5. Create __init__.py to make it an importable package
if [[ ! -f "$DEST/__init__.py" ]]; then
    touch "$DEST/__init__.py"
fi

echo "Done."
