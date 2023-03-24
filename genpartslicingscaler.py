import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input')
parser.add_argument('-o', '--output', help='Output')
args = parser.parse_args()

with open('datalist.txt') as f: #Opens file to be read
    lines = f.readlines() #Puts lines of file into list
    f.close()

with open('notnewdatalist.txt', 'w') as j: #Opens file to be written in
    count = 0 
    for i in lines: #Loops through list of lines, "i" is a specific line
        
        ttdil = i.strip()
        string = 'python cutevents.py -i ' + args.input + ttdil + " -o " + args.output + "postcuts" + str(count) + ".root" #String to be written
        j.write(string + "\n") #Writes string to file
        count += 1
    
    j.close()