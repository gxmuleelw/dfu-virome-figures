def neighbourhood(tree,target,labels,nearest,anchors):
    focal=next(t for t in tree.get_terminals() if t.name==target)
    distances={t.name:tree.distance(focal,t) for t in tree.get_terminals()}
    order=sorted(distances,key=lambda n:(distances[n],n))
    nearest_ids=[n for n in order if n!=target][:nearest]
    anchor_ids=[n for n in order if labels[n]['category']=='ICTV reference'][:anchors]
    chosen={target}|set(nearest_ids)|set(anchor_ids)
    rows=[]
    for n in order:
        reasons=[]
        if n==target:reasons.append('target')
        if n in nearest_ids:reasons.append('nearest_by_tree_distance')
        if n in anchor_ids:reasons.append('nearest_official_reference_context')
        rows.append({'sequence_id':n,'distance_to_target':distances[n],'shown_as_individual_tip':n in chosen,
                     'selection_reason':';'.join(reasons),'zero_or_tied_distance_note':'Ties ordered lexically only for display; no resolved biological ordering implied.'})
    return chosen,rows

def union_bases(intervals):
    merged=[]
    for a,b in sorted((min(a,b),max(a,b)) for a,b in intervals):
        if merged and a<=merged[-1][1]+1:merged[-1][1]=max(merged[-1][1],b)
        else:merged.append([a,b])
    return sum(b-a+1 for a,b in merged),merged
