import os, sys, re
from numpy import einsum
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

counterNode = 0
hc = {}
# initialize a vector of vector of int

mf = [[None] for i in range(0, 26)]

def index(id, set):
    global hc
    if id not in hc:
        hc[id] = set
    return hc[id]

def insert(id, vec):
    global mf
    if (id < 1):
        id -= 1
    elif (id < 4):
        id -= 2
    elif (id < 6):
        id -= 3
    else:
        id -= 4   
    
    if mf[id][0] == None:
        mf[id].append(vec)
        mf[id].pop(0)
    else:
        mf[id].append(vec)

def plot_heatmap(path2, mat2, mat, avg_mat2, avg_mat):
    marks = np.array(mat2)
    
    names=['1', '2', '3', '4', '5', '6', '7']
    subjects=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26']
    
    plt.xticks(ticks=np.arange(len(names)),labels=names,rotation=90)
    plt.yticks(ticks=np.arange(len(subjects)),labels=subjects)
    hm=plt.imshow(marks, cmap='Blues',interpolation="nearest")
    plt.colorbar(hm)
    plt.savefig(path2 + 'hm1.png')

    fig = plt.figure()
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

    marks2 = np.array(mat)
    
    names2=['ei', 'ej', 'ek']
    subjects2=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26']
    
    plt.xticks(ticks=np.arange(len(names2)),labels=names2,rotation=90)
    plt.yticks(ticks=np.arange(len(subjects2)),labels=subjects2)
    hm2=plt.imshow(marks2, cmap='Blues',interpolation="nearest")
    plt.colorbar(hm2)
    plt.savefig(path2 + 'hm2.png')

    fig = plt.figure()
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

    # avg matrix

    marks3 = np.array(avg_mat2)
    
    names3=['1', '2', '3', '4', '5', '6', '7']
    subjects3=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26']
    
    plt.xticks(ticks=np.arange(len(names3)),labels=names3,rotation=90)
    plt.yticks(ticks=np.arange(len(subjects3)),labels=subjects3)
    hm3=plt.imshow(marks3, cmap='Blues',interpolation="nearest")
    plt.colorbar(hm3)
    plt.savefig(path2 + 'hm3.png')

    fig = plt.figure()
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

    marks4 = np.array(avg_mat)
    
    names4=['ei', 'ej', 'ek']
    subjects4=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26']
    
    plt.xticks(ticks=np.arange(len(names4)),labels=names4,rotation=90)
    plt.yticks(ticks=np.arange(len(subjects4)),labels=subjects4)
    hm4=plt.imshow(marks4, cmap='Blues',interpolation="nearest")
    plt.colorbar(hm4)
    plt.savefig(path2 + 'hm4.png')
    # ax = plt.subplots()
    # bar of frequency of each column
    # ax.bar(range(len(subjects)), np.sum(marks, axis=0), align='center')

    fig = plt.figure()
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

d = {
    'id': [0],
    'hyperedges': [
        [1, 2, 3],
    ],
}

# path first argument is the name of the folder to create
path1 = '/home/isis/MoCHy/analysis_result/' + sys.argv[1] + '/'
    
try: 
    os.mkdir(path1)
    print("Created folder " + path1)
except OSError as error: 
    print(error) 

print("creating hashmap...")
with open('/home/isis/MoCHy/intermediate/all_tweets_clean.txt') as infile :
#with open('/home/isis/MoCHy/testing/test1.txt') as infile :
    counter = 0
    for line in infile:
        if not line.strip(): continue  # skip the empty line
        
        numbers = re.findall(r'\d+', line)
        #convert numbers to int
        numbers = [int(i) for i in numbers]
        my_set = set(numbers)
        index(counter, my_set)
        counter += 1
print("done")

print("inserting all sets...")
with open('/home/isis/MoCHy/new_results_he/all_tweets_he') as infile :
    counter = 0
    for line in infile:
        if not line.strip(): continue  # skip the empty line
        # place all the numbers in a list
        numbers = re.findall(r'\d+', line)
        # pick the first number
        id = int(numbers.pop(0))
        numbers = ','.join(numbers)   
        # convert numbers to int
        numbers = [int(i) for i in numbers.split(',')]
        insert(id, numbers)
print("done")

# ei, ej, ek
mat = [[0 for x in range(3)] for y in range(26)]
# intersections of sets
mat2 = [[0 for x in range(7)] for y in range(26)]
# avg ei, ej, ek
avg_mat = [[0 for x in range(3)] for y in range(26)]
# avg intersections
avg_mat2 = [[0 for x in range(7)] for y in range(26)]

total_size_motifs = 0
# total size of each vector in mf
for i in range(0, 26):
    if mf[i] != [None] :
        total_size_motifs += len(mf[i])

print("doing the math...")
for i in range(0, 26):
    avg_size_ei = 0
    avg_size_ej = 0
    avg_size_ek = 0
    counter = 0
    temp = [0 for x in range(7)]
    if mf[i] != [None] :
        for j in range(0, len(mf[i])):
            #id hyperdges
            ei = mf[i][j][0]
            ej = mf[i][j][1]
            ek = mf[i][j][2]
            
            #avg_size += (len(hc[ei]) + len(hc[ej]) + len(hc[ek]))
            #counter += 3
            avg_size_ei += len(hc[ei])
            avg_size_ej += len(hc[ej])
            avg_size_ek += len(hc[ek])
            temp[0] += len(hc[ei] - hc[ej] - hc[ek])
            temp[1] += len(hc[ej] - hc[ek] - hc[ei])
            temp[2] += len(hc[ek] - hc[ei] - hc[ej])
            temp[3] += len(hc[ei] & hc[ej] - hc[ek])
            temp[4] += len(hc[ej] & hc[ek] - hc[ei])
            temp[5] += len(hc[ek] & hc[ei] - hc[ej])
            temp[6] += len(hc[ei] & hc[ej] & hc[ek])

        avg_mat2[i][0] = temp[0] / len(mf[i])
        avg_mat2[i][1] = temp[1] / len(mf[i])
        avg_mat2[i][2] = temp[2] / len(mf[i])
        avg_mat2[i][3] = temp[3] / len(mf[i])
        avg_mat2[i][4] = temp[4] / len(mf[i])
        avg_mat2[i][5] = temp[5] / len(mf[i])
        avg_mat2[i][6] = temp[6] / len(mf[i])    

        mat2[i][0] = ((temp[0] / len(mf[i])) * len(mf[i])) / total_size_motifs
        mat2[i][1] = ((temp[1] / len(mf[i])) * len(mf[i])) / total_size_motifs
        mat2[i][2] = ((temp[2] / len(mf[i])) * len(mf[i])) / total_size_motifs
        mat2[i][3] = ((temp[3] / len(mf[i])) * len(mf[i])) / total_size_motifs
        mat2[i][4] = ((temp[4] / len(mf[i])) * len(mf[i])) / total_size_motifs
        mat2[i][5] = ((temp[5] / len(mf[i])) * len(mf[i])) / total_size_motifs
        mat2[i][6] = ((temp[6] / len(mf[i])) * len(mf[i])) / total_size_motifs

        avg_size_ei /= len(mf[i])
        avg_size_ej /= len(mf[i])
        avg_size_ek /= len(mf[i])

        avg_mat[i][0] = avg_size_ei
        avg_mat[i][1] = avg_size_ej
        avg_mat[i][2] = avg_size_ek
        # for cnt in range(0, 7):
        #     mat[i][cnt] = avg_size * len(mf[i]) / total_size_motifs
        mat[i][0] = (avg_size_ei * len(mf[i])) / total_size_motifs
        mat[i][1] = (avg_size_ej * len(mf[i])) / total_size_motifs
        mat[i][2] = (avg_size_ek * len(mf[i])) / total_size_motifs
    # else we don't have any of this motif
    else:
        # for at in range(0, 7):
        #     mat[i][at] = 0
        mat[i][0] = 0
        mat[i][1] = 0
        mat[i][2] = 0
        avg_mat[i][0] = 0
        avg_mat[i][1] = 0
        avg_mat[i][2] = 0
        mat2[i][0] = 0
        mat2[i][1] = 0
        mat2[i][2] = 0
        mat2[i][3] = 0
        mat2[i][4] = 0
        mat2[i][5] = 0
        mat2[i][6] = 0
        avg_mat2[i][0] = 0
        avg_mat2[i][1] = 0
        avg_mat2[i][2] = 0
        avg_mat2[i][3] = 0
        avg_mat2[i][4] = 0
        avg_mat2[i][5] = 0
        avg_mat2[i][6] = 0

print("done")
print("writing results on file...")
with open(path1 + 'temp.txt', 'w') as outfile:
    outfile.write('ei, ej, ek\n')
    for i in range(0, 26):
        outfile.write(str(mat[i]) + '\n')
    outfile.write('\n\n')
    outfile.write('intersections\n')
    outfile.write('1, 2, 3, 4, 5, 6, 7\n')
    for i in range(0, 26):
        outfile.write(str(mat2[i]) + '\n')
    outfile.write(str(mat2[i]) + '\n\n')

with open(path1 + 'matrix.txt', 'w') as outfile:
    outfile.write('avg mat 1\n')
    for i in range(0, 26):
        outfile.write(str(avg_mat[i]) + '\n')
    outfile.write('\n\n')
    outfile.write('avg mat 2\n')
    for i in range(0, 26):
        outfile.write(str(avg_mat2[i]) + '\n')

print("done")

path2 = path1 + '/heatmaps/'
    
try: 
    os.mkdir(path2)
    print("Created subfolder " + path2)
except OSError as error: 
    print(error) 

print("creating heatmaps...")
plot_heatmap(path2, mat2, mat, avg_mat2, avg_mat)
print("done")

print("cleaning up...")
#remove [] from the txt file
with open(path1 + 'temp.txt', 'r') as infile, open(path1 + 'all_tweets_analysis.txt', 'w') as outfile:
    for line in infile:
        outfile.write(line.replace('[', '').replace(']', ''))

if os.path.exists(path1 + 'temp.txt'):
    os.remove(path1 + 'temp.txt') 

print("Terminated!")

# clean txt from empty lines
# with open('realpoli_clean_to_mochy.txt') as infile, open('output.txt', 'w') as outfile:
#     for line in infile:
#         if not line.strip(): continue  # skip the empty line
#         outfile.write(line)  # non-empty line. Write it to output

