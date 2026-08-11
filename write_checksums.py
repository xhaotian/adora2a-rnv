#!/usr/bin/env python3
from pathlib import Path
import hashlib
import pandas as pd
root=Path(__file__).resolve().parent
rows=[]
for path in sorted(root.rglob('*')):
    if path.is_file() and '.git' not in path.parts and path.name != 'SHA256SUMS.tsv':
        rows.append({'path':str(path.relative_to(root)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size})
pd.DataFrame(rows).to_csv(root/'SHA256SUMS.tsv',sep='\t',index=False)
print(f"Wrote checksums for {len(rows)} files")
