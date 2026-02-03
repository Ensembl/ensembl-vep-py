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

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if coverage file exists
COVERAGE_FILE="${1:-coverage.info}"
OUTPUT_DIR="${2:-.github/badges}"

if [[ ! -f "$COVERAGE_FILE" ]]; then
    echo -e "${RED}Error: Coverage file '$COVERAGE_FILE' not found${NC}" >&2
    exit 1
fi

echo -e "${GREEN}Processing coverage file: $COVERAGE_FILE${NC}"

# Extract coverage percentage using awk
COVERAGE=$(awk -F: '
    /^LH:/ { hit += $2 }
    /^LF:/ { total += $2 }
    END { 
        if (total > 0) 
            printf "%.2f", (hit/total)*100 
        else 
            print "0.00"
    }
' "$COVERAGE_FILE")

# Extract line counts for reporting
LINES_HIT=$(grep "^LH:" "$COVERAGE_FILE" | awk -F: '{sum+=$2} END {print sum}')
LINES_TOTAL=$(grep "^LF:" "$COVERAGE_FILE" | awk -F: '{sum+=$2} END {print sum}')

echo -e "${GREEN}Coverage: $COVERAGE% ($LINES_HIT/$LINES_TOTAL lines)${NC}"

# Determine badge color based on coverage
if (( $(echo "$COVERAGE >= 90" | bc -l) )); then
    COLOR="brightgreen"
    COLOR_HEX="44cc11"
elif (( $(echo "$COVERAGE >= 70" | bc -l) )); then
    COLOR="yellow"
    COLOR_HEX="dfb317"
elif (( $(echo "$COVERAGE >= 50" | bc -l) )); then
    COLOR="orange"
    COLOR_HEX="fe7d37"
else
    COLOR="red"
    COLOR_HEX="e05d44"
fi

echo -e "${GREEN}Badge color: $COLOR${NC}"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Create shields.io endpoint JSON
cat > "$OUTPUT_DIR/coverage.json" <<EOF
{
  "schemaVersion": 1,
  "label": "coverage",
  "message": "${COVERAGE}%",
  "color": "${COLOR}"
}
EOF

echo -e "${GREEN}Created badge JSON: $OUTPUT_DIR/coverage.json${NC}"

# Create a simple HTML badge (alternative)
cat > "$OUTPUT_DIR/coverage.html" <<EOF
<!DOCTYPE html>
<html>
<head>
    <meta http-equiv="refresh" content="0; url=https://img.shields.io/badge/coverage-${COVERAGE}%25-${COLOR}">
</head>
<body>
    <p>Redirecting to badge...</p>
</body>
</html>
EOF

# Create a markdown file with coverage info
cat > "$OUTPUT_DIR/coverage.md" <<EOF
# Coverage Report

- **Coverage:** ${COVERAGE}%
- **Lines Hit:** ${LINES_HIT}
- **Lines Total:** ${LINES_TOTAL}
- **Color:** ${COLOR}

Last updated: $(date -u +"%Y-%m-%d %H:%M:%S UTC")
EOF

echo -e "${GREEN}Created markdown summary: $OUTPUT_DIR/coverage.md${NC}"

# Export variables for GitHub Actions
if [[ -n "$GITHUB_OUTPUT" ]]; then
    echo "coverage=$COVERAGE" >> "$GITHUB_OUTPUT"
    echo "color=$COLOR" >> "$GITHUB_OUTPUT"
    echo "lines_hit=$LINES_HIT" >> "$GITHUB_OUTPUT"
    echo "lines_total=$LINES_TOTAL" >> "$GITHUB_OUTPUT"
fi

echo -e "${GREEN}✓ Coverage badge generation complete${NC}"