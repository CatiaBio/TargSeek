#!/bin/bash
# MAFFT WSL2 wrapper to handle /dev/stderr issues
# This script creates a temporary stderr redirection and runs MAFFT

# Create a temporary directory for stderr handling
TEMP_DIR=$(mktemp -d)
TEMP_STDERR="$TEMP_DIR/stderr.tmp"

# Export the temp stderr for MAFFT to use
export MAFFT_STDERR="$TEMP_STDERR"

# Run MAFFT with stderr redirection
exec 2>"$TEMP_STDERR"
/usr/bin/mafft "$@"
MAFFT_EXIT=$?

# Clean up
rm -rf "$TEMP_DIR"

exit $MAFFT_EXIT