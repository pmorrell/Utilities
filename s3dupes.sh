#!/usr/bin/env bash
# usage: s3dupes.sh                         scan every bucket you can list
#        s3dupes.sh s3://bucket [s3://b2/prefix ...]   scan specific buckets/prefixes
# Finds objects with the same size and ETag (across all scanned locations).
set -euo pipefail

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

if [[ $# -eq 0 ]]; then
  # columns: DATE TIME s3://bucket
  targets=()
  while IFS= read -r b; do targets+=("$b"); done < <(s3cmd ls | awk '{print $NF}')
else
  targets=("$@")
fi
[[ ${#targets[@]} -gt 0 ]] || { echo "no buckets found" >&2; exit 1; }

# one line per object: etag <TAB> size <TAB> url (empty objects are skipped)
for t in "${targets[@]}"; do
  echo "scanning $t ..." >&2
  if ! s3cmd ls -r --list-md5 "$t" 2> "$tmp/err" | awk '
      $4 ~ /^[0-9a-f]{32}(-[0-9]+)?$/ && $3 > 0 {
        print $4 "\t" $3 "\t" substr($0, index($0, "s3://"))
      }' >> "$tmp/all"; then
    echo "warning: failed to list $t: $(head -n1 "$tmp/err")" >&2
  fi
done

[[ -s $tmp/all ]] || { echo "no objects found" >&2; exit 0; }

LC_ALL=C sort -t$'\t' -k1,1 -k2,2n -k3,3 "$tmp/all" |
awk -F'\t' '
  function flush(   i) {
    if (n > 1)
      for (i = 1; i <= n; i++) print (n - 1) * size "\t" etag "\t" size "\t" n "\t" url[i]
  }
  $1 != etag || $2 != size { flush(); etag = $1; size = $2; n = 0 }
  { url[++n] = $3 }
  END { flush() }
' | LC_ALL=C sort -t$'\t' -s -k1,1nr -k2,2 -k3,3n |
awk -F'\t' '
  function human(b,   i, u) {
    split("B KiB MiB GiB TiB PiB", u, " "); i = 1
    while (b >= 1024 && i < 6) { b /= 1024; i++ }
    return sprintf("%.2f %s", b, u[i])
  }
  $2 != etag || $3 != size {
    etag = $2; size = $3; groups++; wasted += $1
    printf "\n%s copies, %d bytes each, %s reclaimable (ETag %s)\n", $4, $3, human($1), $2
  }
  { files++; print "  " $5 }
  END {
    if (!groups) { print "No duplicates found."; exit }
    printf "\n%d duplicate groups, %d files, %d bytes (%s) reclaimable\n", groups, files, wasted, human(wasted)
  }'

