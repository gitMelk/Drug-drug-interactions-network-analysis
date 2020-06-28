import numpy as np
import csv
import pandas as pd

#%% Load the raw data
data = pd.read_csv('ChCh-Miner_durgbank-chem-chem.tsv.csv', sep=',',dtype='str',encoding="utf-8")
info_col = 5
N_col = data.shape[1]
N_rows = data.shape[0]
data = data.fillna(value=0)

print('Total interactions:',N_rows)
#%% Put every drug in a single list

all_drugs = [];

for i in range(N_rows):
    for j in range(N_col):
        if(data.iloc[i][j] != 0):
            all_drugs.append(data.iloc[i][j])

#%% Get a file with the drug ID and drug name 
drug_vec = []
id_vec = []
counter = 1
with open('drugs_nodes.txt', 'w', encoding="utf-8") as f:
    f.write('Id\tLabel\n')
    for drug_tmp in all_drugs:
        if drug_tmp not in drug_vec:
            
            f.write(str(counter) + '\t' + drug_tmp + '\n')      
            drug_vec.append(drug_tmp)
            id_vec.append(counter)
            counter = counter + 1
#%% Replace a drug name with a ID 

for i in range(N_rows):
    for j in range(N_col):
        if(data.iloc[i][j] != 0):
           data.iloc[i][j] = drug_vec.index(data.iloc[i][j]) + 1 
#%% Build the edge file
           
link_all_source = []
link_all_target = []
with open('drugs_edges.txt', 'w', encoding="utf-8") as f:
    f.write('Source\tTarget\n')
    
    for x in range(N_rows):
        counter = 0
        for y in range(N_col):
            if(counter == 0):
                first_tmp = data.iloc[x][0]
                counter = counter + 1
            elif (counter>0 and data.iloc[x][y] != 0):
                                   
                    if(first_tmp < data.iloc[x][y]):
                        f.write(str(first_tmp) + '\t' + str(data.iloc[x][y]) + '\n')    
                    else:
                        f.write(str(data.iloc[x][y]) + '\t' + str(first_tmp) + '\n')
                    # link_all_source.append(first_tmp)
                    # link_all_target.append(data.iloc[x][y])

