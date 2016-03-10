from __future__ import print_function, division, absolute_import
import re

with open("filenames.txt",'r') as f:
    filenames = f.read()

filenames = filenames.split("f1")

#Regex used to flexibly identify wells, rows, and columns from filenames
well_regex = re.compile("(?<=_)\d+(?=-)") #Matches any number of digits preceded by "_" and followed by "-"
column_regex = re.compile("(?<=_)\d{3}(?=\.)") #Matches three digits preceded by "_" and followed by "."
row_regex = re.compile("(?<=_)\d{3}(?=_)") #Matches three digits preceded by "_" and followed by "_"

wellDict = {}

for f in filenames:
    #Extract positioning information from filename with regex
    wellID = int(well_regex.search(f).group())
    rowNumber = int(row_regex.search(f).group())
    columnNumber = int(column_regex.search(f).group())

    #If there is no key for wellID in wellDict -> create a dict of properties
    if wellDict.get(wellID) == None:
        wellDict[wellID] = {'nrows':0, 'ncols':0, 'positions':{}, 'files':[]}

    #Populate Properties
    wellDict[wellID]['nrows'] = max(rowNumber+1, wellDict[wellID]['nrows'])
    wellDict[wellID]['ncols'] = max(columnNumber+1, wellDict[wellID]['ncols'])

    #List of filenames for the well
    wellDict[wellID]['files'].append(f)

    #Dict with (row, column):filename
    wellDict[wellID]['positions'][(rowNumber, columnNumber)]=f

k = wellDict.keys()[0]
print(wellDict[k])
