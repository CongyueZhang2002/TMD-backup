import os


for i in range(201):
    if len(str(i)) == 1:
        numb = '00'+str(i)
    elif len(str(i)) == 2:
        numb = '0'+str(i)
    else:
        numb = str(i)
    src = './fit/tools/readdata.f'
    dest = "replica_"+numb
    f = open(src, "r")
    copy = open(dest+"/tools/readdata.f", "w")
    for line in f:
        #if 'expdata' in line:
        #    line = line.replace('expdata/', 'expdata/rep'+numb+'/')
        if 'NLINES ' in line:
            line = '      NLINES = 1\n'
        copy.write(line)
    f.close()
    copy.close()
