fa = open("summary.txt", "a")
fr = open("summary.txt", "r")

good = 0
optimal = 0
bad = 0

for line in fr.readlines():
    word_list = line.split(" --- ")

    if len(word_list) == 1: # Last line contains only "-"
        break

    rating = word_list[1]

    if rating == "Good\n":
        good += 1
    elif rating == "Optimal!\n":
        optimal += 1
    elif rating == "Needs improvement\n":
        bad += 1
    
fa.write("\nTOTAL COUNT\n")
fa.write("Optimal: ")
fa.write(str(optimal))
fa.write("\n")
fa.write("Good: ")
fa.write(str(good))
fa.write("\n")
fa.write("Needs improvement: ")
fa.write(str(bad))
fa.write("\n")