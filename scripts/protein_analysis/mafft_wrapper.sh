#!/bin/bash
# MAFFT wrapper to handle WSL2 /dev/stderr issues
# This script fixes the stderr redirection problem in WSL2 environments

# Create temporary stderr file if needed
if [ ! -e /dev/stderr ]; then
    # Create a temporary file to act as stderr
    TEMP_STDERR=$(mktemp)
    exec 2>"$TEMP_STDERR"
fi

# Run MAFFT with all provided arguments
exec /usr/bin/mafft "$@"