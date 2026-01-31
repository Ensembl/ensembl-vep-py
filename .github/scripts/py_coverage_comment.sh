#!/usr/bin/env bash

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -euo pipefail

LCOV_FILE=$1
COVERAGE_THRESHOLD=${COVERAGE_THRESHOLD:-90}
MARKER="<!-- python-lcov-coverage-report -->"

: "${GITHUB_TOKEN:?Missing GITHUB_TOKEN}"
: "${GITHUB_REPOSITORY:?Missing GITHUB_REPOSITORY}"
: "${GITHUB_EVENT_PATH:?Missing GITHUB_EVENT_PATH}"

# Get PR number
PR_NUMBER=$(jq -r '.pull_request.number' "$GITHUB_EVENT_PATH")

if [[ -z "$PR_NUMBER" || "$PR_NUMBER" == "null" ]]; then
  echo "Not running on a pull request. Exiting."
  exit 0
fi

# ---- Parse LCOV totals ----
TOTAL_LINES=$(grep -h '^LF:' "$LCOV_FILE" | awk -F: '{sum += $2} END {print sum}')
HIT_LINES=$(grep -h '^LH:' "$LCOV_FILE" | awk -F: '{sum += $2} END {print sum}')

if [[ -z "$TOTAL_LINES" || "$TOTAL_LINES" -eq 0 ]]; then
  echo "Failed to parse LCOV file"
  exit 1
fi

TOTAL_COVERAGE=$(awk -v hit="$HIT_LINES" -v total="$TOTAL_LINES" \
  'BEGIN { printf "%.2f", (hit / total) * 100 }')

COVERAGE_OK=$(awk -v cov="$TOTAL_COVERAGE" -v th="$COVERAGE_THRESHOLD" \
  'BEGIN { print (cov >= th) ? "yes" : "no" }')

if [[ "$COVERAGE_OK" == "yes" ]]; then
  STATUS_ICON="✅"
  STATUS_TEXT="Coverage meets threshold"
else
  STATUS_ICON="❌"
  STATUS_TEXT="Coverage below ${COVERAGE_THRESHOLD}% threshold"
fi

# ---- Per-file coverage table ----
TABLE=$(awk '
BEGIN {
  print "| File | Coverage |"
  print "|------|----------|"
}
$1 ~ /^SF:/ { file = substr($1,4); lf=0; lh=0 }
$1 ~ /^LF:/ { lf = substr($1,4) }
$1 ~ /^LH:/ { lh = substr($1,4) }
$1 ~ /^end_of_record/ {
  if (lf > 0) {
    printf "| %s | %.2f%% |\n", file, (lh / lf) * 100
  }
}
' "$LCOV_FILE")

# ---- Comment body ----
COMMENT_BODY=$(cat <<EOF
$MARKER
### 🐍 Python Coverage Report (coverage.py / lcov)

**Total coverage:** **${TOTAL_COVERAGE}%**  
**Lines:** $HIT_LINES / $TOTAL_LINES  
**Status:** $STATUS_ICON $STATUS_TEXT

<details>
<summary>Per-file coverage</summary>

$TABLE

</details>
EOF
)

# ---- GitHub API ----
API_URL="https://api.github.com"
COMMENTS_URL="$API_URL/repos/$GITHUB_REPOSITORY/issues/$PR_NUMBER/comments"
AUTH_HEADER="Authorization: Bearer $GITHUB_TOKEN"

EXISTING_COMMENT_ID=$(curl -s \
  -H "$AUTH_HEADER" \
  -H "Accept: application/vnd.github+json" \
  "$COMMENTS_URL" \
  | jq -r ".[] | select(.body | contains(\"$MARKER\")) | .id")

if [[ -n "$EXISTING_COMMENT_ID" ]]; then
  echo "Updating existing Python coverage comment"
  curl -s -X PATCH \
    -H "$AUTH_HEADER" \
    -H "Accept: application/vnd.github+json" \
    "$API_URL/repos/$GITHUB_REPOSITORY/issues/comments/$EXISTING_COMMENT_ID" \
    -d "$(jq -nc --arg body "$COMMENT_BODY" '{body: $body}')"
else
  echo "Creating Python coverage comment"
  curl -s -X POST \
    -H "$AUTH_HEADER" \
    -H "Accept: application/vnd.github+json" \
    "$COMMENTS_URL" \
    -d "$(jq -nc --arg body "$COMMENT_BODY" '{body: $body}')"
fi

# ---- Fail CI if coverage too low ----
if [[ "$COVERAGE_OK" != "yes" ]]; then
  echo "Coverage ${TOTAL_COVERAGE}% is below threshold ${COVERAGE_THRESHOLD}%"
  exit 1
fi
