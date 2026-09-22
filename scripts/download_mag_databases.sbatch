#!/usr/bin/env bash
#SBATCH --job-name=pks-mag-dbs
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=download_mag_databases.%j.log

set -euo pipefail

# Destination for the downloaded databases. Export PKS_DB_ROOT before submitting:
#   PKS_DB_ROOT=/path/to/dbs sbatch scripts/download_mag_databases.sbatch
shared_root="${PKS_DB_ROOT:?PKS_DB_ROOT must be set to the directory that will hold the CheckM2 and GTDB-Tk databases}"
download_root="${shared_root}/downloads"
checkm2_root="${shared_root}/checkm2_v1.1.0"
gtdbtk_root="${shared_root}/gtdbtk_r232"

checkm2_archive="${download_root}/checkm2_database.v1.1.0.tar.gz"
gtdbtk_archive="${download_root}/gtdbtk_r232_data.tar.gz"

checkm2_url=https://zenodo.org/api/records/14897628/files/checkm2_database.tar.gz/content
gtdbtk_url=https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz

mkdir -p "${download_root}" "${checkm2_root}" "${gtdbtk_root}"
exec 9>"${download_root}/pks_mag_databases.lock"
flock -n 9 || {
    echo "Another MAG database installation is already running." >&2
    exit 1
}

echo "Downloading CheckM2 v1.1.0 database"
curl -L --fail --retry 10 --retry-delay 30 --continue-at - \
    "${checkm2_url}" -o "${checkm2_archive}"
printf '07c10655620843b517d0df0c160d911f  %s\n' "${checkm2_archive}" | md5sum -c -

if ! find "${checkm2_root}" -type f -name 'uniref100.KO.1.dmnd' -size +0c -print -quit | grep -q .; then
    tar -xzf "${checkm2_archive}" -C "${checkm2_root}"
fi
checkm2_db=$(find "${checkm2_root}" -type f -name 'uniref100.KO.1.dmnd' -size +0c -print -quit)
test -n "${checkm2_db}"
printf '%s\n' "${checkm2_db}" > "${checkm2_root}/DATABASE_PATH.txt"

echo "Downloading GTDB-Tk R232 reference package"
curl -L --fail --retry 10 --retry-delay 30 --continue-at - \
    "${gtdbtk_url}" -o "${gtdbtk_archive}"
printf '25a59e0352b1fd150c589f56559767d4  %s\n' "${gtdbtk_archive}" | md5sum -c -

if [[ ! -f "${gtdbtk_root}/.installed.md5" ]]; then
    top_count=$(tar -tzf "${gtdbtk_archive}" | awk -F/ 'NF {print $1}' | sort -u | wc -l)
    if [[ "${top_count}" -eq 1 ]]; then
        tar -xzf "${gtdbtk_archive}" -C "${gtdbtk_root}" --strip-components=1
    else
        tar -xzf "${gtdbtk_archive}" -C "${gtdbtk_root}"
    fi
    printf '25a59e0352b1fd150c589f56559767d4\n' > "${gtdbtk_root}/.installed.md5"
fi

test -d "${gtdbtk_root}/markers"
test -d "${gtdbtk_root}/taxonomy"
test -d "${gtdbtk_root}/pplacer"

{
    printf 'CheckM2\t1.1.0\t%s\t07c10655620843b517d0df0c160d911f\n' "${checkm2_db}"
    printf 'GTDB-Tk\tR232\t%s\t25a59e0352b1fd150c589f56559767d4\n' "${gtdbtk_root}"
} > "${shared_root}/pksProfiler_MAG_DATABASES.tsv"

echo "MAG database download and checksum validation complete."
echo "CheckM2: ${checkm2_db}"
echo "GTDB-Tk: ${gtdbtk_root}"

