#!/usr/bin/env sh

guix environment -m ../scripts/gcc-latest.scm --container -- ./gcc-build.sh
