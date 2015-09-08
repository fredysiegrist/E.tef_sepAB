__author__ = 'fsiegris'
#date: 07.09.2015

#from recipe-577444-1 import getColumns
import csv

#test importing the dag file and print out the first significant line
try:
    print('Starting test phase:')
    file = open('../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go', 'r+')
    print(file.readline())
    print(file.readline())
    file.close()
except NameError:
    print('File not readable')
except FileNotFoundError:
    print('File not located where it should')
except:
    print('I give up on you!')
finally:
    print('Starting proof of principle:')
#getColumns('../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go')

#print(cols)

#get out the different columns and calculate the difference of the scaffold and reference chromosome coordinates

lol = list(csv.reader(open('../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go', 'r+'), delimiter='\t'))
for n in range(1, 10):
    print(lol[n][0]+' '+str(int(lol[n][2])-int(lol[n][6]))+' '+lol[n][4])

print('\n\nTest phase over! Starting real life: Раз Два Три\n\n')

#calculate for same Sorghum chromosome gene the distances