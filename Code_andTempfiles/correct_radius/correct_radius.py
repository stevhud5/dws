# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 12:40:15 2018

@author: sdh
"""

#replace wrong value in file

import fileinput

path = 'C:\\Users\\sdh\\Documents\\Projects\\2 DWS\\J_J\\temp'
fnameBase='\\023-25C-Tscan_'
i=0
fname_in=path+fnameBase+str(i)+'.dat'
fname_out=path+fnameBase+str(i)+'c.dat'


with open(fname_in, 'r+') as f_in:
    #goto line
    for index in range (14):
        f_in.readline()
    line_to_edit = f_in.readline()
    print(line_to_edit)
    new_line = line_to_edit.replace('465', '300', 1)
    print(new_line)
    
    #f_in.writelines(line_to_edit)
    