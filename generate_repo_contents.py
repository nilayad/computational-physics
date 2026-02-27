#!/usr/bin/env python3
"""
Generates repo_contents.txt containing the raw content of every file
in this repository (excluding the .git directory and the output file itself).
"""

import os

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(REPO_ROOT, "repo_contents.txt")
OUTPUT_FILENAME = os.path.basename(OUTPUT_FILE)

SEPARATOR = "=" * 80


def collect_files(root):
    """Yield all regular (non-symlink) file paths under root, skipping .git and the output file."""
    abs_output = os.path.abspath(OUTPUT_FILE)
    abs_root = os.path.abspath(root)
    for dirpath, dirnames, filenames in os.walk(root):
        # Skip .git directory in-place so os.walk won't descend into it
        dirnames[:] = [d for d in dirnames if d != ".git"]
        for filename in sorted(filenames):
            filepath = os.path.join(dirpath, filename)
            abs_filepath = os.path.abspath(filepath)
            # Skip symlinks to avoid following paths outside the repository
            if os.path.islink(filepath):
                continue
            # Ensure the file is actually inside the repository root
            if not abs_filepath.startswith(abs_root + os.sep) and abs_filepath != abs_root:
                continue
            # Skip the output file itself to avoid self-inclusion
            if abs_filepath == abs_output:
                continue
            yield filepath


def main():
    # Validate that the output file will be written inside the repository
    abs_output = os.path.abspath(OUTPUT_FILE)
    abs_root = os.path.abspath(REPO_ROOT)
    if not abs_output.startswith(abs_root + os.sep):
        raise ValueError(f"OUTPUT_FILE {abs_output!r} is outside REPO_ROOT {abs_root!r}")

    files = sorted(collect_files(REPO_ROOT))
    with open(OUTPUT_FILE, "w", encoding="utf-8", errors="replace") as out:
        out.write("REPOSITORY CONTENTS\n")
        out.write(f"Root: {REPO_ROOT}\n")
        out.write(f"Total files: {len(files)}\n")
        out.write(SEPARATOR + "\n\n")

        for filepath in files:
            rel_path = os.path.relpath(filepath, REPO_ROOT)
            out.write(f"FILE: {rel_path}\n")
            out.write(SEPARATOR + "\n")
            try:
                with open(filepath, "r", encoding="utf-8", errors="replace") as f:
                    out.write(f.read())
            except (OSError, UnicodeDecodeError) as exc:
                out.write(f"[Could not read file: {exc}]\n")
            out.write("\n" + SEPARATOR + "\n\n")

    print(f"Written: {OUTPUT_FILE}")
    print(f"Files included: {len(files)}")


if __name__ == "__main__":
    main()
