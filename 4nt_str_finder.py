#individual identification tool based on Short Tandem Repeat detection (AGAT, AATG, TATC)
import sys
import csv
from sys import argv
from csv import reader

#Load genome text and database csv and check for additional/lack of args
argc = len(argv)
if(argc!=3):
    print("ERROR: Please enter STR library (csv) and suspect genome (txt) respectively.")
    sys.exit()

#(0): Iterate through csv, save every row as a list
genomedatabase = open(argv[1])
database = csv.reader(genomedatabase, delimiter=',')
matrix = []
for row in database:
    matrix.append(row)

#Genome text file
genome= open(argv[2], "r")
genome = genome.read()

#(1): Substring counter: count len(largest consecutive STR sequence)

#Substring count arrays
AGAT = []
AGAT_counter = 0

AATG = []
AATG_counter = 0

TATC = []
TATC_counter = 0

#STR Counter: indexes across genome, scans for 4-Nt substrings AGAT, AATG, TATC and stores consistent reads in an array for each STR
a = 0
b = 4

for a in range(0,len(genome)-4):
    if("AGAT" in genome[a:b]):
        AGAT_counter+=1

        AATG.append(AATG_counter)
        AATG_counter = 0
        TATC.append(TATC_counter)
        TATC_counter = 0

    elif("AATG" in genome[a:b]):
        AATG_counter+=1

        TATC.append(TATC_counter)
        TATC_counter = 0
        AGAT.append(AGAT_counter)
        AGAT_counter = 0

    elif("TATC" in genome[a:b]):
        TATC_counter+=1

        AATG.append(AATG_counter)
        AATG_counter = 0
        AGAT.append(AGAT_counter)
        AGAT_counter = 0
    a+=1
    b+=1

#Ensure no values get left out from the list
AGAT.append(AGAT_counter)
TATC.append(TATC_counter)
AATG.append(AATG_counter)

#Take highest vals in each array as longest repeat
AGAT_val = max(AGAT)
AATG_val = max(AATG)
TATC_val = max(TATC)

#(2): Compare genome substring counts to database
suspect_found = 0
for x in range(1,len(matrix)):
    if(int(matrix[x][1]) == AGAT_val and int(matrix[x][2]) == AATG_val and int(matrix[x][3]) == TATC_val):
        suspect_found=1
        print("Predicted suspect: {}" .format(matrix[x][0]))
    x+=1

if(suspect_found==0):
    print("Suspect genome does not match any database entry")