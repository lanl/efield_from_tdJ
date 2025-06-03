import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import re

# Get input/output
i_file = sys.argv[1]
o_file = sys.argv[2]

if not os.path.exists(i_file):
    raise ValueError("Invalid input file: {i_file}")

content = []
component_labels = []
components = []
times = np.array([])
time_label = ""

with open(i_file, 'r') as f:
    content = f.readlines()
    hdr = content[0]

    # Get indices of arrays
    time_idx = -1
    times = np.zeros((len(content)-1,))
    for i, tok in enumerate(re.split('\s{2,}', hdr)):
        if "Time" in tok:
            time_idx = 0
            time_label = tok
        elif len(tok) > 0:
            component_labels.append(tok)
            components.append(np.zeros((len(content)-1,)))

    # Verify a time was found
    if time_idx == -1:
        raise RuntimeError("Misformatted input. Could not find a time axis.")

    # Extract the time series
    for i, line in enumerate(content):
        if i == 0:
            continue

        toks = line.split()
        # Otherwise, add to arrays
        times[i-1] = toks[0]
        for j in range(len(toks) - 1):
            components[j][i-1] = toks[j + 1]

fig, ax = plt.subplots(1, 1, figsize=(12,8))

# Plot the results                                                                                                      
fields = components

colors = ['black', 'blue', 'red']
for i, fe in enumerate(fields):
    ax.plot(times, fe, color=colors[i], label=component_labels[i])

ax.legend()
ax.set_xlabel(time_label)

plt.savefig(o_file)
