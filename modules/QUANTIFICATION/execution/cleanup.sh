# Remove default output dir and smk internal dir

cleanup () {
    rc=$?
    rm -rf .snakemake/
    rm -rf ../output/
    cd "$user_dir"
    echo "Exit status: $rc"
}
trap cleanup EXIT

set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands

user_dir=$PWD
pipeline_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd "$pipeline_dir"
