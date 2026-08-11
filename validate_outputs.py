#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
from PIL import Image

root=Path(__file__).resolve().parent
r=root/'review_remediation'
checks=[]
def add(name,ok,detail=''): checks.append({'check':name,'status':'PASS' if ok else 'FAIL','detail':str(detail)})
effects=pd.read_csv(r/'02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_STUDY_EFFECTS.tsv',sep='\t')
meta=pd.read_csv(r/'02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_META_RESULTS.tsv',sep='\t').iloc[0]
comp=pd.read_csv(r/'02_mouse_compartment_audit/MOUSE_COMPARTMENT_META_RESULTS.tsv',sep='\t').set_index('analysis_set')
models=pd.read_csv(r/'04_human_rebuild/human_primary_model_full_results.tsv',sep='\t')
add('mouse_effect_count',len(effects)==16,len(effects)); add('mouse_positive_count',int((effects.hedges_g>0).sum())==14)
add('mouse_pooled_g',np.isclose(meta.pooled_hedges_g,1.12244963229513,atol=1e-10),meta.pooled_hedges_g)
add('compartment_direction',comp.loc['WHOLE_RETINA_OR_LYSATE','pooled_hedges_g']>0 and comp.loc['ENRICHED_OR_ISOLATED_CELL_COMPARTMENT','pooled_hedges_g']>0)
add('human_models',set(models.model)=={'M0','M1','M2'} and len(models)==18,len(models))
for n in range(1,5):
 p=r/'07_figures'/f'Fig{n}.tif'; im=Image.open(p)
 add(f'Fig{n}_600dpi',im.info.get('dpi')==(600.0,600.0),im.info.get('dpi'))
 add(f'Fig{n}_LZW',im.tag_v2.get(259)==5,im.tag_v2.get(259)); add(f'Fig{n}_size',p.stat().st_size<10_000_000,p.stat().st_size)
out=pd.DataFrame(checks); out.to_csv(root/'FRESH_RUN_CHECKS.tsv',sep='\t',index=False)
summary={'checks':len(out),'pass':int(out.status.eq('PASS').sum()),'fail':int(out.status.eq('FAIL').sum())}
(root/'FRESH_RUN_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
if summary['fail']: raise SystemExit(1)
