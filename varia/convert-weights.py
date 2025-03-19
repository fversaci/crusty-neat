#!/usr/bin/env python
import numpy as np
import json
with open("neat_quality_score_model.json") as f:
    j = json.load(f)
s = np.array(j['seed_weights'])
w = np.array(j['weights_from_one'])
# normalize seeds
see = s/s.sum()
see
# set 0 as sink state
w[:,0,0] = 1
w[:,39:42,0] = 1
# normalize weights
wei = w.astype(np.float64)  # Convert to double floats
wei /= wei.sum(axis=2, keepdims=True)  # Normalize each row
# save as json
d = dict()
d['seeds'] = see.tolist()
d['weights'] = wei.tolist()
with open("qs_model.json", "w") as f:
    json.dump(d, f)
