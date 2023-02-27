with open('datalist.txt') as f: #Opens file to be read
    lines = f.readlines() #Puts lines of file into list
    f.close()

with open('notnewdatalist.txt', 'w') as j: #Opens file to be written in
    count = 0 
    for i in lines: #Loops through list of lines, "i" is a specific line
        string = 'python cutevents.py -i ' + i + "-o postcuts" + str(count) + ".root" #String to be written
        j.write(string) #Writes string to file
        count += 1
    
    j.close()