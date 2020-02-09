fin = open("bm.txt", "r")               # Results from running mlp 10 times for each instance
fout = open("results.txt", "w+")        # Output file
ftarget = open("target.txt", "r")       # Target average results
fsummary = open("summary.txt", "w+")    # Results summary, looking at off percentage

isFirst = True
flag = 0 # 0 is cost 1 is time

for line in fin.readlines():

    #Treating read line
    no_endl = line.split("\n")
    my_list = no_endl[0].split()
    
    if len(my_list) == 1: # Line is an instance's name

        instance_name = my_list[0]
        
        if isFirst == False: # There's no cost or time data at the first name
            avg_cost = sum_cost/10
            avg_time = sum_time/10

            off_pct = ( (avg_cost - target) / target ) * 100

            # Verbose - Writes info to results.txt
            fout.write("\nAvg cost: ")
            fout.write(str(avg_cost))
            fout.write("\nAvg time: ")
            fout.write(str(avg_time))
            fout.write("\nOff percentage: ")
            fout.write(str(off_pct))
            fout.write("\n\n")

            # Summary - Writes info to summary.txt
            fsummary.write(" --- ")
            if off_pct == 0:
                fsummary.write("Optimal!")
            elif off_pct <= 0.5:
                fsummary.write("Good")
            else:
                fsummary.write("Needs improvement")
            fsummary.write("\n")

        fout.write(instance_name)
        fsummary.write(instance_name)

        # Gets target value from target.txt for given instance
        # so the error % can be calculated later.
        for tgt in ftarget.readlines():
            
            #Treating read line
            tgt_list = tgt.split(":")
            tgt_name = tgt_list[0]

            if instance_name == tgt_name:
                
                # Continue treating the line (after if to save resources)
                tgt_value_endl = tgt_list[1]
                tgt_value_list = tgt_value_endl.split("\n")
                tgt_value = tgt_value_list[0]
                
                target = float(tgt_value)

                ftarget.seek(0,0)

                break

        #Reset total cost and time
        sum_cost = 0
        sum_time = 0



    else: # Else line is cost or time
        num = my_list[1]
        
        if flag == 0:
            sum_cost += int(num)
            flag = 1
        elif flag == 1:
            sum_time += float(num)
            flag = 0
        
    isFirst = False


fin.close()
fout.close()
ftarget.close()
fsummary.close()
