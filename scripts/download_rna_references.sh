#!/usr/bin/env bash
# Download and verify the references for --sample_type rna (public GRCh38.d1.vd1 + GENCODE v36; ~12 GB):
#   GRCh38.d1.vd1.splice_sr.mmi   minimap2 2.31 index built with -x splice:sr   -> --rna_splice_db
#   gencode.v36.junctions.bed     paftools.js gff2bed of GENCODE v36 (BED12)   -> --rna_junctions_bed
#   gencode.v36.annotation.gtf.gz the source annotation
#
#   bash scripts/download_rna_references.sh /path/to/dbs/rna_grch38_splice_sr_gencode_v36
#
# Downloads the shared Google Drive folder with gdown (pip install gdown), then checks every file
# against the MD5s below, which are fixed here rather than read from the download. Already-verified
# files are not fetched again. If gdown is unavailable or Drive rate-limits the download, fetch the
# folder in a browser, put the three files in the destination, and rerun: the script then only verifies.
# To build the same files yourself instead, see docs/running/wes_rna.md.
set -euo pipefail

FOLDER_URL=https://drive.google.com/drive/folders/10np5NSeAPRpHz1a22drGybs4enTP-OdR
dest="${1:?usage: download_rna_references.sh <destination directory>}"

declare -A MD5=(
    [GRCh38.d1.vd1.splice_sr.mmi]=a3eef1343261ca39eed311bc0a8b35c5
    [gencode.v36.junctions.bed]=26df7e2772b9ac99319960a3b4c10af0
    [gencode.v36.annotation.gtf.gz]=c03931958d4572148650d62eb6dec41a
)

verified() { [[ -s "$dest/$1" ]] && [[ $(md5sum "$dest/$1" | cut -d' ' -f1) == "${MD5[$1]}" ]]; }

mkdir -p "$dest"
missing=()
for f in "${!MD5[@]}"; do verified "$f" || missing+=("$f"); done

if (( ${#missing[@]} )); then
    if ! command -v gdown >/dev/null 2>&1; then
        echo "gdown not found. Install it (pip install gdown), or download $FOLDER_URL by hand into $dest and rerun." >&2
        exit 1
    fi
    staging="$dest/.download"
    mkdir -p "$staging"
    gdown --folder "$FOLDER_URL" -O "$staging" --remaining-ok
    # The shared folder may hold the files directly or one level down; take each by name.
    for f in "${missing[@]}"; do
        src=$(find "$staging" -maxdepth 3 -type f -name "$f" | sed -n 1p)
        [[ -n $src ]] || { echo "ERROR: $f not found in the downloaded folder" >&2; exit 1; }
        mv "$src" "$dest/$f"
    done
fi

status=0
for f in "${!MD5[@]}"; do
    if verified "$f"; then echo "ok        $f"; else echo "MISMATCH  $f" >&2; status=1; fi
done
(( status == 0 )) || { echo "Checksum failure: delete the bad file(s) and rerun." >&2; exit 1; }
rm -rf "$dest/.download"
echo
echo "Use:  --rna_splice_db $dest/GRCh38.d1.vd1.splice_sr.mmi --rna_junctions_bed $dest/gencode.v36.junctions.bed"
