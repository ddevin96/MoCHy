import os, sys, re

with open('/home/isis/MoCHy/email.hgf') as infile, open('/home/isis/MoCHy/email_clean.txt', 'w') as outfile:
    for line in infile:
        if not line.strip(): continue  # skip the empty line
        # place all the numbers in a list
        numbers = re.findall(r'\d+', line)
        #convert numbers to int
        # create a string from numbers with comma separator
        line = ','.join(numbers) + '\n'
        outfile.write(str(line))
