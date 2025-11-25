"""
This script looks at an example tar file to prep grabbing instruments

"""
import pandas as pd
import tarfile
from datetime import datetime
import gzip

tardir = "data/5-MVP_download/phs002453.MVP_R4.1000G_AGR.GIA.PheCodes_CirculatorySystem_batch1.analysis-PI.MULTI.tar"

# inspect tar contents
with tarfile.open(tardir, 'r') as tar:
    members = tar.getmembers()
    print(f"Total members in tar: {len(members)}\n")
    print("Contents:")
    for member in members:
        print(f"  {member.name} (size: {member.size} bytes, isdir: {member.isdir()})")

# Extract and read the first .txt.gz file
print("\n" + "="*60)
print("Extracting first .txt.gz file...")
print("="*60 + "\n")

with tarfile.open(tardir, 'r') as tar:
    members = tar.getmembers()
    # Find first .txt.gz file
    txt_gz_file = next((m for m in members if m.name.endswith('.txt.gz')), None)
    
    if txt_gz_file:
        print(f"Found: {txt_gz_file.name}\n")
        
        # Extract the file from tar
        f = tar.extractfile(txt_gz_file)
        
        # Decompress and read
        with gzip.GzipFile(fileobj=f) as gz:
            content = gz.read().decode('utf-8')
            lines = content.split('\n')
            print(f"Total lines: {len(lines)}\n")
            print("First 50 lines:")
            print("\n".join(lines[:50]))
    else:
        print("No .txt.gz files found in tar")

pass