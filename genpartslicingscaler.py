import argparse
import os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input')
parser.add_argument('-o', '--output', help='Output')
args = parser.parse_args()
files = os.listdir(args.input)
'''
with open('datalist.txt') as f: #Opens file to be read
    lines = f.readlines() #Puts lines of file into list
    f.close()
'''
with open('notnewdatalist.txt', 'w') as j: #Opens file to be written in
    count = 0 
    #for i in lines: #Loops through list of lines, "i" is a specific line
    for i in files: 
        
        #ttdil = i.strip()
        string = 'python genpartslicing.py -i ' + args.input + i + " -o " + args.output + "postcuts" + str(count) + ".root" #String to be written
        j.write(string + "\n") #Writes string to file          ^ was ttdil
        count += 1
    
    j.close()