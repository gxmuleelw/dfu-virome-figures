import math
import numpy as np

def summarize(raw, raw_rpkm, breadth, clean, samples):
    counts=raw*(breadth>=.5);rpkm=raw_rpkm*(breadth>=.5)
    raw_totals=np.sum(raw,axis=0,dtype=np.int64);totals=np.sum(counts,axis=0,dtype=np.int64)
    rpkm_totals=[math.fsum(rpkm[:,j]) for j in range(len(samples))]
    assert np.all(raw_totals<=clean)
    composition=(totals>0)&(np.asarray(rpkm_totals)>0)
    assert np.array_equal(totals>0,np.asarray(rpkm_totals)>0)
    return raw_totals, totals, rpkm_totals, composition
