#!/usr/bin/env bash
set -euo pipefail

# 1. Get Python's architecture-specific site-packages directory
USER_SITE=$(python3 -c "import site; print(site.getusersitepackages())")

# 2. Define destination folder
DEST="$USER_SITE/rosi"

echo "Installing .so and .py modules into: $DEST"

# 3. Create destination directory
mkdir -p "$DEST"

# 4. Copy all .so and .py files
cp build-cmake/cpp/python_binding/*.so "$DEST/"
cp python/modules/*.py "$DEST/"

# 5. Create __init__.py to make it an importable package
if [[ ! -f "$DEST/__init__.py" ]]; then
    touch "$DEST/__init__.py"
fi

echo "Done."
