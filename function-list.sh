#!/bin/sh -e

cat << EOM > functions.md
# Biolibc function list

| Function | Purpose |
|----------|---------|
EOM

auto-man2man Man/* | awk -F - '{ printf("| %s | %s |\n", $1, $2); }' \
    >> functions.md

# grip --export functions.md
# firefox ./functions.html
