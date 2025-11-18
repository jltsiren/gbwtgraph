#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# set-pggname.py name graph.gfa
# set-pggname.py $(pggname graph.gfa)

import sys

if len(sys.argv) != 3:
    print("Usage: set-pggname.py name graph.gfa")
    sys.exit(1)

name = sys.argv[1]
gfa_file = sys.argv[2]

# Read the lines from the GFA file
lines = []
with open(gfa_file, 'r') as f:
    lines = f.readlines()

# If there is already a header line with a NM:Z tag, do nothing.
for i in range(len(lines)):
    if lines[i].startswith('H'):
        if 'NM:Z:' in lines[i]:
            sys.exit(0)

# Check if the first line is a file header
file_header = None
if lines and lines[0].startswith('H'):
    if 'VN:Z:' in lines[0]:
        file_header = lines[0]
        lines = lines[1:]

# Write the possible file header, the new header line, and the rest of the lines back to the file
new_header = f"H\tNM:Z:{name}\n"
with open(gfa_file, 'w') as f:
    if file_header:
        f.write(file_header)
    f.write(new_header)
    f.writelines(lines)
