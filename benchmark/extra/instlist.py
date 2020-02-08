# Simple program to get all instances names
# and write ":" next to them

fin = open("bm.txt", "r")
fout = open("target.txt", "w+")

for line in fin.readlines():
    word_list = line.split(":")
    
    if len(word_list) == 1:
        no_endl = word_list[0].split("\n")

        fout.write(no_endl[0])
        fout.write(":")
        fout.write("\n")

fin.close()
fout.close()