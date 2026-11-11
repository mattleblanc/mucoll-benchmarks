#!/bin/bash
# Build the calibration-hit DDG4 plugin.
# Run INSIDE the mucoll sim container that will run ddsim: the plugin must be
# built against the same DD4hep/Geant4. Each software stack gets its own build,
# named after its DD4hep install (e.g. install/dd4hep-1.37-5u6u23u/), so builds
# for different images (v3.1, main, ...) coexist and setup.sh picks the right one.
#
#   ./build.sh            # -> ./install/<stack>/lib/libCalibHits.so (+ .components)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -z "${DD4hep_DIR:-}" ]]; then
    set +u
    # shellcheck disable=SC1091
    source /opt/setup_mucoll.sh
    set -u
fi
# shellcheck source=stack_tag.sh
source "$HERE/stack_tag.sh"
STACK="$(calibhits_stack_tag)"

BUILD_DIR="${CALIBHITS_BUILD_DIR:-$HERE/build/$STACK}"
INSTALL_DIR="${CALIBHITS_INSTALL_DIR:-$HERE/install/$STACK}"

cmake -S "$HERE" -B "$BUILD_DIR" \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DCMAKE_INSTALL_LIBDIR=lib
cmake --build "$BUILD_DIR" -j "${CALIBHITS_JOBS:-4}"
cmake --install "$BUILD_DIR"

echo "Installed to $INSTALL_DIR/lib (stack $STACK). Before ddsim:  source $HERE/setup.sh"
