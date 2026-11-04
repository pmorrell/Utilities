#!/usr/bin/env bash
# usage: s3tree.sh s3://bucket/prefix/            file count and total size
#        s3tree.sh -t s3://bucket/prefix/         summary plus tree and ETags
#        s3tree.sh -c s3://a/prefix/ s3://b/x/    compare two directories
set -euo pipefail

usage() { sed -n '2,4p' "$0" | sed 's/^# //' >&2; exit 1; }

# prints: relative_path <TAB> size <TAB> etag
list() {
  local url="${1%/}/"
  s3cmd ls -r --list-md5 "$url" | awk -v n="${#url}" '
    $4 ~ /^[0-9a-f]{32}(-[0-9]+)?$/ {
      key = substr($0, index($0, "s3://"))
      print substr(key, n + 1) "\t" $3 "\t" $4 "\t" $1 " " $2
    }' | LC_ALL=C sort -t$'\t' -k1,1
}

summary() {
  awk -F'\t' -v u="$1" '
    { n++; s += $2
      if (min == "" || $4 < min) min = $4
      if ($4 > max) max = $4 }
    END {
      split("B KiB MiB GiB TiB PiB", unit, " "); h = s; i = 1
      while (h >= 1024 && i < 6) { h /= 1024; i++ }
      printf "%s\n  files:  %d\n  size:   %d bytes (%.2f %s)\n", u, n, s, h, unit[i]
      if (n) printf "  oldest: %s (last modified)\n  newest: %s (last modified)\n", min, max
    }' "$2"
}

mode=summary
case "${1:-}" in
  -t) mode=tree; shift ;;
  -c) mode=compare; shift ;;
  -h|--help|"") usage ;;
esac

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

if [[ $mode == compare ]]; then
  [[ $# -eq 2 ]] || usage
  list "$1" > "$tmp/a"; list "$2" > "$tmp/b"
  summary "$1" "$tmp/a"; summary "$2" "$tmp/b"; echo
  awk -F'\t' '
    NR == FNR { size[$1] = $2; etag[$1] = $3; next }
    { seen[$1] = 1
      if (!($1 in size)) print "only in B:     " $1
      else if (size[$1] != $2) print "size differs:  " $1
      else if (etag[$1] != $3) print "ETag differs:  " $1 " (same size; may still be identical if multipart)"
    }
    END { for (k in size) if (!(k in seen)) print "only in A:     " k }
  ' "$tmp/a" "$tmp/b" | LC_ALL=C sort > "$tmp/diff"
  if [[ -s $tmp/diff ]]; then
    cat "$tmp/diff"; exit 1
  fi
  echo "Identical (same paths, sizes and ETags)"
  exit 0
fi

[[ $# -eq 1 ]] || usage
list "$1" > "$tmp/a"
summary "$1" "$tmp/a"

if [[ $mode == tree ]]; then
  echo
  cut -f1 "$tmp/a" | tree --fromfile -a --noreport .
  echo; echo "ETags:"
  awk -F'\t' '{ print $3 "  " $1 }' "$tmp/a"
  echo; echo "Duplicate ETags:"
  cut -f3 "$tmp/a" | sort | uniq -d
fi

