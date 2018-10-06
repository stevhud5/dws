# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 12:40:15 2018

@author: sdh
"""

#replace wrong value in file

import fileinput

path = 'C:\\Users\\sdh\\Documents\\Projects\\2 DWS\\J_J\\temp'
fnameBase='\\023-25C-Tscan_'

minNr=1
maxNr=17

for i in range(minNr,maxNr+1):
    fname_in=path+fnameBase+str(i)+'.dat'
    
    #edit text
    with open(fname_in, 'r+') as f_in:
        #goto line
        for index in range (14):
            f_in.readline()
        line_to_edit = f_in.readline()
        
    line_to_edit = line_to_edit.rstrip('\r\n')
    new_line = line_to_edit.replace('465', '300', 1)
    change = {line_to_edit: new_line}
    print(line_to_edit)
    print(new_line)
    print()
    
    #record change in file
    for line in fileinput.input(fname_in, inplace=True):
        line = line.rstrip('\r\n')
        print(change.get(line,line))
        