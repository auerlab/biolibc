#!/bin/sh -e

cat << EOM > functions.md
# Biolibc function list

Each function below is documented by a man page.  To view the documentation,
install biolibc using your chosen package manager and run \`man function\`
(e.g. \`man bl_fastq_read\`).

| Function | Purpose |
|----------|---------|
EOM

auto-man2man Man/* | awk -F - '$1 !~ "BIOLIBC" { printf("| %s | %s |\n", $1, $2); }' \
    >> functions.md

grip --export functions.md
firefox ./functions.html
