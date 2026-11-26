import math
from collections import defaultdict

def number15(x):return float(format(float(x),'.15g'))
def text_bool(x):return 'TRUE' if x else 'FALSE'

def assign(ids,records):
    idset=set(ids)
    edges=defaultdict(list)
    for h in records:
        assert h['Virus'] in idset and 'g__' in h['Host genus'] and math.isfinite(float(h['Confidence score']))
        if float(h['Confidence score'])>=90:edges[h['Virus']].append((h['Host genus'].rsplit('g__',1)[1],float(h['Confidence score'])))
    hosts={}
    for vid in ids:
        z=edges[vid];best=max((x[1] for x in z),default=None);winners={x[0] for x in z if x[1]==best};group='Unassigned' if not z else next(iter(winners)) if len(winners)==1 else 'Ambiguous'
        status='unassigned' if group=='Unassigned' else 'tied_best_genera' if group=='Ambiguous' else 'best_genus_at_score_ge_90'
        hosts[vid]=dict(predicted_host_genus=group,host_iPHoP_best_score='NA' if best is None else number15(best),host_n_qualifying_genera=len({x[0] for x in z}),host_prediction_status=status,has_host_genus_prediction=text_bool(status=='best_genus_at_score_ge_90'))
    return hosts
