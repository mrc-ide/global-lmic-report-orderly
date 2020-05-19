#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}
DEFAULT_SHORT="FALSE"
SHORT_RUN=${2:-$DEFAULT_SHORT}

echo "*** Date: $DATE"

echo "*** Short Run: $SHORT_RUN"

echo "*** Updating country list"
./update_run_sh.R $DATE

echo "*** Running country reports"

# Parallel
grep -E '^[A-Z]{3}\s*' countries | \
parallel --progress -j 32 ./orderly run lmic_reports_google iso3c={} date=$DATE short_run=$SHORT_RUN
