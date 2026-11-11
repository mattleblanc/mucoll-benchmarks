# Name of the software stack in the running container, used to keep plugin
# builds for different images apart: the DD4hep install directory name with a
# short hash, e.g. "dd4hep-1.37-5u6u23u". Needs DD4hep_DIR (/opt/setup_mucoll.sh).
calibhits_stack_tag() {
    if [[ -z "${DD4hep_DIR:-}" ]]; then
        echo "ERROR: DD4hep_DIR not set; source /opt/setup_mucoll.sh first" >&2
        return 1
    fi
    basename "$DD4hep_DIR" | sed -E 's/^(dd4hep-[^-]+-.{7}).*$/\1/'
}
