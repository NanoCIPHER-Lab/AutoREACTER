import os

# Root directory for the LUNAR installation.
# Configure via the LUNAR_ROOT_DIR environment variable; defaults to None if unset.
LUNAR_ROOT_DIR = os.environ.get("LUNAR_ROOT_DIR") or None
