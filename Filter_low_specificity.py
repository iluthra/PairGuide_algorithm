#Filter 1 removes any guides with Specificity < 20
#creates a new file with cut start and cut end columns
#gives values for Sguide bar graph to input into excel
#created a file with cut start differences to use in excel to plot


import operator
Database = open("Database_GATA3.txt", 'r')
Database = Database.readlines()[1:]
chromfree = []
for line in Database:
    
    chromfree.append([x for x in line.split('\t')])
print "Total possible guides:"
print len(chromfree)
filtered_list = []
i = 0
#Set filter threshold here
filter_threshold = 20.00
while i < len(chromfree):
    
    if float(chromfree[i][8]) > filter_threshold:
        
        filtered_list.append(chromfree[i])
        i = i +1
    else:
        i = i + 1
        
print "After removing guides with Specificity Score less than 20 we are left with:"
print len(filtered_list)
filtered_list2 = []
counter = 0
while counter < len(filtered_list):
    info = []
    startend = filtered_list[counter][0].split(':')[1]
    start = (startend.split('-')[0])
    start = int(start)
    end = (startend.split('-')[1])
    end = int(end)

    filtered_list[counter][14] = filtered_list[counter][14].rstrip('\n')

    if filtered_list[counter][3] == 'fwd':
        filtered_list[counter].append(end-3)
        filtered_list[counter].append(end + 1 -3)
        counter = counter + 1
    elif filtered_list[counter][3] == 'rev':
        filtered_list[counter].append(start-1 + 3)
        filtered_list[counter].append(start +3)
        counter = counter + 1

#print filtered_list[0:5]
#sort by cut start
sorted_list = sorted(filtered_list, key=operator.itemgetter(15))

j = 0

S23 = 0
S34 = 0 
S45 = 0
S56 = 0
S67 = 0
S78 = 0
S89 = 0
S910 = 0
k = 0
#print sorted_list[0][8]
while k < len(sorted_list):
    if (float(sorted_list[k][8]) > 20.00) & (float(sorted_list[k][8]) < 30.00):
        S23 = S23 + 1
        k = k +1
    elif (float(sorted_list[k][8]) >= 30.00) & (float(sorted_list[k][8]) < 40.00):
        S34 = S34 + 1
        k = k +1
    elif (float(sorted_list[k][8]) >= 40.00) & (float(sorted_list[k][8]) < 50.00):
        S45 = S45 + 1
        k = k +1
    elif (float(sorted_list[k][8]) >= 50.00) & (float(sorted_list[k][8]) < 60.00):
        S56 = S56 + 1
        k = k +1
    elif (float(sorted_list[k][8]) >= 60.00) & (float(sorted_list[k][8]) < 70.00):
        S67 = S67 + 1
        k = k +1
    elif (float(sorted_list[k][8]) >= 70.00) & (float(sorted_list[k][8]) < 80.00):
        S78 = S78 + 1
        k = k +1
    elif (float(sorted_list[k][8]) >= 80.00) & (float(sorted_list[k][8]) < 90.00):
        S89 = S89 + 1
        k = k +1
    elif (float(sorted_list[k][8]) >= 90.00) & (float(sorted_list[k][8]) < 100.00):
        S910 = S910 + 1
        k = k +1


print "Sguide between 20- 30"
print S23
print "Sguide between 30- 40"
print S34
print "Sguide between 40- 50"
print S45
print "Sguide between 50- 60"
print S56
print "Sguide between 60- 70"
print S67
print "Sguide between 70- 80"
print S78
print "Sguide between 80- 90"
print S89
print "Sguide between 90- 100"
print S910




f2 = open('Distance_between_guides_S20_GLs_varying.txt', 'w')
f3 = open('Bad_Regions_Info_S20s_GLs_varying.txt', 'w')
f = open('Filtered_GATA3_database_S20s_GLs_varying.txt','w')
f.write('Location' + '\t' + '\t' + '\t' + '\t' + 'gRNA' + '\t' + '\t' + 'Direction' + '\t' + 'other.perfect.aligns' + '\t' +  'Specificity.Score' + '\t' + '%GC Content' + '\t' + '#MM1'+ '\t' + '#MM2'+ '\t' + '#MM3' + '\t' + 'Cut.Start' + '\t' + 'Cut.End'+  '\n')
f3.write('Location' + '\t' + '\t' + '\t' + '\t' + 'gRNA' + '\t' + '\t' + 'Direction' + '\t' + 'other.perfect.aligns' + '\t' +  'Specificity.Score' + '\t' + '%GC.Content' + '\t' + '#MM1'+ '\t' + '#MM2'+ '\t' + '#MM3' + '\t' + 'Cut.Start' + '\t' + 'Cut.End'+  '\n')

for sublist in sorted_list:
    for item in sublist:
       f.write(str(item) + '\t')
    f.write('\n')

cut_start_diff = []
m = 0
#test = sorted_list[1][15] -sorted_list[0][15]
#print test
#add list that has info about really BAD coverage regions cut_start_diff > 1500
bad_areas = []
count = 0
print len(sorted_list)
while (m < len(sorted_list)-1):
    a = sorted_list[m+1][15] - sorted_list[m][15]
    cut_start_diff.append(a)
    
    if a > 1500:
        #print "BAD"
        #print a
        bad_areas.append(sorted_list[m])
        bad_areas.append(sorted_list[m+1])

    m = m +1

print bad_areas


    
for item in cut_start_diff:
    f2.write(str(item) + '\n')

for element in bad_areas:
    for item2 in element:
       f3.write(str(item2) + '\t')
    f3.write('\n')
    
