# Source (inside the container, after /opt/setup_mucoll.sh) to make the
# CalibHits plugin visible to ddsim. DD4hep finds plugins through the
# .components files on LD_LIBRARY_PATH.
#
# Uses the build for the running software stack (install/<stack>/, see
# build.sh), or $CALIBHITS_INSTALL_DIR if set.
_CALIBHITS_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=stack_tag.sh
source "$_CALIBHITS_HERE/stack_tag.sh"
if [[ -n "${CALIBHITS_INSTALL_DIR:-}" ]]; then
    _CALIBHITS_LIB="$CALIBHITS_INSTALL_DIR/lib"
else
    _CALIBHITS_LIB="$_CALIBHITS_HERE/install/$(calibhits_stack_tag)/lib"
fi
if [[ ! -f "$_CALIBHITS_LIB/libCalibHits.components" ]]; then
    echo "ERROR: no CalibHits build for this software stack ($_CALIBHITS_LIB)." >&2
    echo "       Run $_CALIBHITS_HERE/build.sh inside this image first." >&2
    unset _CALIBHITS_HERE _CALIBHITS_LIB
    return 1
fi
export LD_LIBRARY_PATH="$_CALIBHITS_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
unset _CALIBHITS_HERE _CALIBHITS_LIB
