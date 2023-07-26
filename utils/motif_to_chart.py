#  ['0', 528,1348,0,780,30,8701,0,0,0,33,0,18115,0,0,0,47587,0,96,4,22917,30,520934,0,0,16,134541],
#         ['Eyebrow pencil', 13012, 5067, 3987, 3932],
#         ['Rouge', 11624, 7004, 3574, 5221],
#         ['Pomade', 8814, 9054, 4376, 9256],
#         ['Eyeshadows', 12998, 12043, 4572, 3308],
#         ['Eyeliner', 12321, 15067, 3417, 5432],
#         ['Foundation', 10342, 10119, 5231, 13701],
#         ['Lip gloss', 22998, 12043, 4572, 4008],
#         ['Mascara', 11261, 10419, 6134, 18712]
# 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
import os, sys, re
import pandas as pd
import numpy as np

df = pd.read_csv("/home/isis/MoCHy/motif_data_political_all.csv", header=None)

# add a column in first position with the index of the row
df.insert(0, 'index', range(len(df)))

df['index'] = df['index'].astype(str)
#first_column = df.pop('Name')

df.to_csv("/home/isis/MoCHy/chart.csv", index=False, header=None)

#remove empty lines
with open('/home/isis/MoCHy/chart.csv') as infile, open('/home/isis/MoCHy/chart_final.csv', 'w') as outfile:
    for line in infile:
        if not line.strip(): continue  # skip the empty line
        # add [' ad the start of each line
        #remove the new line from line
        line = line.replace("\n", "")
        #replace first number in line with the same number but with apostrophe and a comma
        line = re.sub(r'^(\d+)', r"\1'", line)
        outfile.write("['" + line + "],\n")
