#Oct. 12, 2016

#this program will use input as Database.txt and a file that contains free chromatin regions to determine
#what percentage of the guide fall into a nucleosome free region

Database = open("Filtered_GATA3_database_S20s_GLs_varying.txt", 'r')
FreeChrom = open("Ad2_13_S9_peaks.xls", 'r')

#remove header line
Database = Database.readlines()[1:]

#print Database
#extract chrom:start-end
Guidelocations = []
gRNA_seq = []
for l in Database:
    Guidelocations.append(l.split('\t')[0])
    gRNA_seq.append(l.split('\t')[2])

i = 0
#print gRNA_seq[1]
#removes '>'
#print Guidelocations
while i < len(Guidelocations):
    Guidelocations[i] = Guidelocations[i][1:len(Guidelocations[i])]
    i = i+1
#print Guidelocations
counter = 0
chrom_o = []
gstart_o = []
gend_o = []
#print len(Guidelocations)
while counter < len(Guidelocations):
    
    chrom_o.append(Guidelocations[counter].split(':')[0])
    startend = Guidelocations[counter].split(':')[1]
    gstart_o.append(startend.split('-')[0])
    gend_o.append(startend.split('-')[1])
    
    counter = counter + 1


#print chrom
#print start
#print end
chromfree = []

filetest = [n for n in FreeChrom.readlines() if not n.startswith('#')]
for k in filetest:
    
    chromfree.append([x for x in k.split('\t')])

#chrom free is a list
#within this list you have lists of info from .xls file
chromfree = chromfree[2:len(chromfree)]

#print len(chrom)
big_counter = 0
new_list = []
final_percent_list = []


gstart_o = map(int,gstart_o)
gend_o = map(int,gend_o)

while big_counter < len(chrom_o):
    current_gRNA_length = len(gRNA_seq[big_counter])
    print big_counter
    chrom = []
    gstart = []
    gend = []
    chrom.append(chrom_o[big_counter])
    gstart.append(gstart_o[big_counter])
    gend.append(gend_o[big_counter])
    #print chrom
    #print gstart
    #print gend
    m = 0
    percent_in_chromatin_free = []
    #l = 0
    while m < len(chromfree):
        l = 0
      
        while l < len(chrom):
            
            
            if chromfree[m][0] == chrom[l]:
            
                chromfree[m][1] = int(chromfree[m][1])
                chromfree[m][2] = int(chromfree[m][2])
                if chromfree[m][1] >= gend[l] :
                    percent_in_chromatin_free.append(0)
                    l = l + 1
                elif gstart[l] >= chromfree[m][2]:
                    percent_in_chromatin_free.append(0)
                    l = l+1
            

                elif gstart[l] >= chromfree[m][1]:
                    if chromfree[m][2] >= (gend[l]):
                        percent_in_chromatin_free.append(100)
                   
                        l = l + 1
                    elif (chromfree[m][2] < (gend[l])):
                        calculation = ((int(chromfree[m][2]) - int(gstart[l])+1)/float(current_gRNA_length))*100
                        if calculation > 0:
                            percent_in_chromatin_free.append(int(calculation))
                            l = l + 1
                     

                elif gend[l] <= chromfree[m][2]:
                    if ((int(gend[l]) - int(current_gRNA_length)) >= chromfree[m][1]):
                        percent_in_chromatin_free.append(100.00)
                        l = l + 1
                        

                    elif ((int(gend[l]) - int(current_gRNA_length)) < chromfree[m][1]):
                        calculation = (int(gend[l] - chromfree[m][1])/float(current_gRNA_length))*100
                    
                        if calculation > 0:
                            percent_in_chromatin_free.append(int(calculation))
                            l = l + 1
                         

                elif gend[l] > chromfree[m][2]:

                    if ((int(gend[l]) - int(current_gRNA_length)) >= chromfree[m][1]):
                        calculation =( (chrom[m][2] - gstart[l])/float(current_gRNA_length))
                        percent_in_chromatin_free.append(int(calculation))
                        l = l + 1
                        

                    elif ((int(gend[l]) - int(current_gRNA_length)) < chromfree[m][1]):
                        calculation =((chromfree[m][2] - chromfree[m][1]+1)/float(current_gRNA_length))*100
                        if calculation > 0:
                            percent_in_chromatin_free.append(int(calculation))
                            l = l + 1
                           


                elif gstart[l] < chrom[m][1]:
                    if chromfree[m][2] >= (gend[l]):
                        calculation = ((gend[l] - chrom[m][1])/float(current_gRNA_length)) * 100
                        percent_in_chromatin_free.append(int(calculation))
                        

                        l = l + 1
                    elif (chromfree[m][2] < (gend[l])):
                        calculation = ((int(chromfree[m][2]) - int(chromfree[m][1]))/float(current_gRNA_length))*100
                        if calculation > 0:
                            percent_in_chromatin_free.append(int(calculation))
                            l = l + 1
                            

                            
            elif chromfree[m][0] != chrom[l]:
                
                percent_in_chromatin_free.append(0)
                l = l +1
        f = 0
        if len(percent_in_chromatin_free) == len(chromfree):

            while f < len(percent_in_chromatin_free):
                if percent_in_chromatin_free[f] != 0:
                    
                    maxpercent = max(percent_in_chromatin_free)
                    
                    new_list.append(maxpercent)
           
                    f = f +1
                f = f+1

            if len(new_list) != (l):
                
                new_list.append(0)
            percent_in_chromatin_free = []
        

        

        m = m + 1
    
    big_counter = big_counter + 1
    
    final_percent_list.append(new_list)
    new_list = []

#print "%chromatin free for each guide is"
#print final_percent_list




output = ""

i = 0
f = open("GATA3_Database_ChromatinRegions_GLs_varying.txt","w")
f.write('Location' + '\t' + '\t' + '\t'  + 'gRNA' + '\t' + '\t' + 'Direction' + '\t' + '#other perfect aligns' + '\t' +  'Specificity Score' + '\t' + '%GC Content' + '\t' + '#MM1'+ '\t' + '#MM2'+ '\t' + '#MM3' + '\t' + 'Cut Start' + '\t' + 'Cut End'+  '%Chromatin Free'+ '\n')

while i < len(Database):
    output = ""
    var = final_percent_list[i][0]
    
    output =  ("%s\t%s" %(Database[i].strip(),var))
    f.write(output+ '\n')
    i = i +1
    
f.close()
