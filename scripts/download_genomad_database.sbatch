#!/usr/bin/env bash
#SBATCH --job-name=pks-genomad-db
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=download_genomad_database.%j.log

set -euo pipefail

# Destination for the downloaded databases. Export PKS_DB_ROOT before submitting:
#   PKS_DB_ROOT=/path/to/dbs sbatch scripts/download_genomad_database.sbatch
shared_root="${PKS_DB_ROOT:?PKS_DB_ROOT must be set to the directory that will hold the geNomad databases}"
download_root="${shared_root}/downloads"
database_root="${shared_root}/genomad_db_v1.9"
archive="${download_root}/genomad_db_v1.9.tar.gz"
url=https://zenodo.org/api/records/14886553/files/genomad_db_v1.9.tar.gz/content
checksum=67244b528bb8bed464d1ca147136d33e

mkdir -p "${download_root}" "${database_root}"
exec 9>"${download_root}/pks_genomad_database.lock"
flock -n 9 || {
    echo "Another geNomad database installation is already running." >&2
    exit 1
}

curl -L --fail --retry 10 --retry-delay 30 --continue-at - "${url}" -o "${archive}"
printf '%s  %s\n' "${checksum}" "${archive}" | md5sum -c -

if [[ ! -f "${database_root}/.installed.md5" ]]; then
    tar -xzf "${archive}" -C "${database_root}" --strip-components=1
    test -f "${database_root}/genomad_db"
    printf '%s\n' "${checksum}" > "${database_root}/.installed.md5"
fi

test -f "${database_root}/genomad_db"
printf 'geNomad\t1.12.0\tdatabase-v1.9\t%s\t%s\n'     "${database_root}" "${checksum}"     > "${shared_root}/pksProfiler_GENOMAD_DATABASE.tsv"

echo "geNomad database installation complete: ${database_root}"
