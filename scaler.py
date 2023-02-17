with open('datalist.txt') as f:
    lines = f.readlines()
    f.close()

with open('notnewdatalist', 'w') as j:
    for i in lines:
        j.write('python cutevents.py -i ' + i)
    
    j.close()