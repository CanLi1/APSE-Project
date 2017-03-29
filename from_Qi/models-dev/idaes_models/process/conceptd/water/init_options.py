"""
This class determines possible initialization options for the model.
Feasible initialization choices are saved to a text file in JSON format.
"""
from __future__ import division
from itertools import chain, combinations
__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

# Determined to be feasible with BARON
feasible_pts = (
	(1, 5), (1, 4), (1, 4, 5), (1, 3), (1, 3, 5),
	(1, 3, 4), (1, 3, 4, 5), (1, 2, 5), (1, 2, 4),
	(1, 2, 4, 5), (1, 2, 3), (1, 2, 3, 5), (1, 2, 3, 4),
	(1, 2, 3, 4, 5)
)

# Determined to be feasible with CONOPT
# feasible_pts = (
# 	(1, 4), (1, 3), (1, 3, 5),
# 	(1, 3, 4), (1, 3, 4, 5), (1, 2, 5), (1, 2, 4),
# 	(1, 2, 4, 5), (1, 2, 3, 5), (1, 2, 3, 4),
# 	(1, 2, 3, 4, 5)
# )

# from Python documentation: powerset recipe. Modified for nonempty elements.
def powerset(base_set):
	"powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
	return chain.from_iterable(combinations(base_set, r) for r in range(1, 4 + 1))

def covering(s):
	"""
	Checks if the given set covers all five possible units
	and that no element of s is completely contained in another element of s
	"""
	return all(sum(1 for x in s if i in x) >= 1 for i in range(1, 5 + 1)) \
		and all(not set(x).issubset(set(y)) for y in s for x in s if x != y)

init_choices = [s for s in powerset(feasible_pts) if covering(s)]

import json
with open('init_opts.txt', 'w') as outfile:
	json.dump(init_choices, outfile)

with open('init_opts.txt', 'r') as infile:
	in_data = json.loads(infile.read())
