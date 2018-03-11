import os 
import sys
import math
import decimal
import numpy as np
from collections import Counter
from collections import OrderedDict

# GLOBAL VARIABLES
sum_of_length = 0
CHAIN = ['#']
amino_acid = ''
TITLE = ''

#LISTS USED
list_4  = []
list_3  = []
list_2  = []
ligands = []
# FUNCTION TO CALCULATE FREQUENCY
def frequency_amino_acid_acid(line):
    frequency = dict(Counter(line.split()))
    frequency = OrderedDict(sorted(frequency.items(), key=lambda a:[0]))
    return frequency

# FUNCTION TO CALCULATE THE DIHEDRAL ANGLE 
def calc_dihedral_angle(a):
    p1 = np.array(a[1])
    p3 = np.array(a[3])
    p0 = np.array(a[0])
    p2 = np.array(a[2])
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
        # normalize b1 so that it does not influence magnitude of vector
        # v = projection of b0 onto plane perpendicular to b1
        # w = projection of b2 onto plane perpendicular to b1
        # angle between v and w in a plane is the dihedral angle
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))
    
# OPEN OUTPUT FILE
output_file = open("2wsc_output.txt",'w')
#READING PDB FILE AND STORING IN LIST
for line in open('2wsc.pdb'):
    list_1 = list(line.split()) 
    list_3.append(list_1)


for i in list_3:
    if (i[0] == "SEQRES"):
        if (CHAIN[-1] != i[2]):
            CHAIN.append(i[2])
        for j in range (0,4,1):
            i.pop(0)
        sum_of_length += len(i)
        for k in i:
            amino_acid +=" " + k
    
    elif (i[0] == "FORMUL"):   
        ligands.append(i[2])
    
    elif (i[0] == "TITLE"):
            list_2 = []
            list_2.append(i)
            list_2[0].pop(0)

            if (TITLE == ''):
                TITLE += ' '.join(list_2[0])
            
            else:
                list_2[0].pop(0)
                TITLE += ' ' + ' '.join(list_2[0])
    
    else:  
        continue


ligands.sort()
ligandsandcons = ' '.join(ligands)

frequency1 = frequency_amino_acid_acid(amino_acid)

if 'UNK' in frequency1:
    a = int(frequency1['UNK'])
else:
    a = 0;

CHAIN.pop(0)
CHAIN.sort()
CHAINconst = ','.join(CHAIN)

output_file.write (TITLE)
output_file.write ('\n')
str1 = 'LENGTH:'+'\t'+str(sum_of_length)+'\n'
output_file.write (str1)

str1 = 'CHAINS:'+'\t'+str(len(CHAIN))+'\t'+str(CHAINconst)+'\n'
output_file.write (str1)

for key in frequency1:
        list_4.append(str(key)+'\t'+str(round((float(frequency1[key])/sum_of_length),10))+'\n')
        
list_4.sort()
for string in list_4:        
        output_file.write (string)

str1 = 'UNKNOWN: '+str(a)
output_file.write (str1)
output_file.write ('\n')
output_file.write ('LIGANDS:')
output_file.write (ligandsandcons)
output_file.write ('\n')
C = []
N = []
CA = []
ch = []

for i in list_3:
    if (i[0] == "ATOM"):
        just_another_list = []
        just_another_list.extend([float(i[6]),float(i[7]),float(i[8])])
        if (i[2] == "C"):
            C.append(just_another_list)
        elif (i[2] == "N"):
            ch.append(i[4])
            N.append(just_another_list)
        elif (i[2] == "CA"):
            CA.append(just_another_list)
        else:
            continue
    else:
        continue           
strl = "CHAIN-"+ch[1]
phi = [strl]
phi.append('NA')
psi = ['']
omega = ['']

for i in range(0,len(C)-1,1):
    phi_list = [];
    phi_list.extend([C[i],N[i+1],CA[i+1],C[i+1]])
    phi.append(calc_dihedral_angle(phi_list))
    psi_list = [];
    psi_list.extend([N[i],CA[i],C[i],N[i+1]])
    psi.append(calc_dihedral_angle(psi_list))
    omega_list = [];
    omega_list.extend([CA[i],C[i],N[i+1],CA[i+1]])
    omega.append(calc_dihedral_angle(omega_list))
    if ( ch[i] != ch [i+1]):
        strl = "CHAIN - " + ch[i+1]
        phi[-1] = (strl)
        phi.append('NA')
        omega.pop(-1)
        psi.pop(-1)
    else:
        continue

for i in (zip(phi,psi,omega)):
    str1 =str(i[0])+'\t\t'+str(i[1]) +'\t\t'+str(i[2])
    output_file.write (str1)
    output_file.write ('\n')
output_file.write ('\n')