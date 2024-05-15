#!/usr/bin/env bash
podman run -it --rm -v $(pwd):/app -w /app coaler "$@"