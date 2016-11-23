#Oct 27 2016
#this is the fastest version yet
#this version will check MM1, MM2 and MM3 for our thresholds and set Specificty to 1 if less than cut-off
import os
import sys
import argparse
import tables
import string
import numpy as np
import cProfile
dna_comp = None
import fileinput

def comp(seq_str):
    """complements the provided DNA sequence and returns it"""
    global dna_comp

    if dna_comp is None:
        dna_comp = string.maketrans("ATCGMRWSYKNatcgmrwsykn",
                                    "TAGCKYWSRMNtagckywsrmn")
    return seq_str.translate(dna_comp)


def revcomp(seq_str):
    """returns reverse complement of provided DNA sequence"""
    return comp(seq_str)[::-1]
    

def str_from_nparray(vals):
    """converts a numpy array into a sequence string"""
    return "".join(chr(x) for x in vals)


def get_seq(chrom,start,end,strand,genome_file):
     
    #seq_h5 = tables.openFile("/iblm/netapp/data1/external/GRC37/GRC37.h5", "r")
    #strand = "+"


    if chrom not in genome_file.root:
        known_chrom = [node.name for node in genome_file.root]
        raise ValueError("unknown chromosome %s, possible chromosomes are: "
                         + ",".join(known_chrom))

    
    # numpy array representation of sequence and convert to string
    chrom_node = genome_file.getNode("/%s" %chrom)
    np_seq = chrom_node[start-1:end]
    seq_str = str_from_nparray(np_seq)

    if strand == "-":
        # reverse-complement sequence
        seq_str = revcomp(seq_str)

    
  
    return seq_str

    # write header and sequence
    #sys.stdout.write(">%s:%d-%d(%s)\n" % (options.chrom, options.start,
                                        #options.end, options.strand))
    #sys.stdout.write("%s\n" % seq_str)

def specificity_calc(listIndex, infoList1,testguide,seq_h5):
    Shitlist = []
    #infoList1 now has everything from .sam file
    current_guide_length = len(testguide)
    while listIndex < len(infoList1):
            refSEQ = get_seq(infoList1[listIndex][0], int(infoList1[listIndex][2]), int(infoList1[listIndex][2]) + (current_guide_length-1) , infoList1[listIndex][1], seq_h5)            

            #call function that takes 
            #len1 = len(refSEQ)
            #len2 = len(guideRna[q])
            mismatchPos = []
            m = 0

            while m < len(refSEQ):
                if refSEQ[m] != testguide[m]:
                    mismatchPos.append(m)
                    m = m + 1
                else:
                    m = m +1


        #for each read that the guide allined to we need what position the mismatches occured at compared to the guide we had
        #using this information calculate SHit for each read that gRNA aligned at and using all the reads calculate Sguide
            if current_guide_length == 20:        
                W = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445,
                    0.508, 0.613, 0.851, 0.731, 0.828, 0.615, 0.804, 0.685, 0.583]
            elif current_guide_length == 21:
                #added a 0
                W = [0, 0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445,
                    0.508, 0.613, 0.851, 0.731, 0.828, 0.615, 0.804, 0.685, 0.583]
            elif current_guide_length == 19:
                    #removed a zero
                W = [0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445,
                    0.508, 0.613, 0.851, 0.731, 0.828, 0.615, 0.804, 0.685, 0.583]    

            d = 0
            l = (current_guide_length-1)


            k = 0
            if len(mismatchPos) == 0:

                #d = 0
                #SHit = (1 - W[mismatchPos[0]])*(1/((1-d*1.0/l)*4)+1)
                #Shitlist.append(SHit)

                listIndex = listIndex + 1
                
            if len(mismatchPos) == 3:
                d = ((mismatchPos[2] - mismatchPos[1]) + (mismatchPos[1] - mismatchPos[0])) / (len(mismatchPos) - 1 )

                a = (1 - W[mismatchPos[2]])*(1 - W[mismatchPos[1]])*(1 - W[mismatchPos[0]])
                b = (1.0/(((1-d*1.0/l)*4.0)+1))
                c = 1.0/9.0
                SHit = a*b*c

                Shitlist.append(SHit)
                listIndex = listIndex + 1
            if len(mismatchPos) == 2:
                d = ((mismatchPos[1] - mismatchPos[0])) / (len(mismatchPos) - 1 )

                a =  (1 - W[mismatchPos[1]])*(1 - W[mismatchPos[0]])
                b = (1/(((1-d*1.0/l)*4.0)+1))
                c = 1.0/4.0
                SHit = a*b*c
                Shitlist.append(SHit)

                listIndex = listIndex + 1
            if len(mismatchPos) == 1:

                d = 0
                SHit = (1 - W[mismatchPos[0]])*(1/((1-d*1.0/l)*4)+1)
                Shitlist.append(SHit)

                listIndex = listIndex + 1

            if len(mismatchPos) > 3:
                listIndex = listIndex + 1


    Sguide = 100 / (1 + sum(Shitlist))
    return Sguide

def finding_guides(seq, inStart, inEnd,gRNA_length):
    guideRna = []
    #has list with RNAs with proper direction
    guideRnaOutput = []
    location = []
    direction = []
    i = 0
    while i < (len(seq)-1):
        
        if seq[i] == "G":
            #print "HI"
            if seq[i+1] == "G":
                #print i
                #only if we are at least 23 bases in can we store a guide RNA, so check and store in list
                if i >= gRNA_length+1:
                    if seq[i-(gRNA_length+1)] == "G":
                        guideRna.append(seq[i-(gRNA_length+1):i-1])
                        
                        guideRnaOutput.append(seq[i-(gRNA_length+1):i-1])
                        #guideRnaPam.append(seq[i-21:i-1] + " " + seq[i-1:i+2])
                        start = (i - (gRNA_length+1)) + inStart
                        end = (i - 1) + inStart -1
                        location.append(str(start) + '-' + str(end))
                        direction.append('fwd')
        #increment counter
        i = i + 1




    ## add loop to look for CCN start once found take reverse complement
    w = 0
    forward = ""
    while w < (len(seq)-1):
        
        if seq[w] == "C":
            
            if seq[w+1] == "C":
                #print w
                #print len(seq)
                #only if we are at least 23 bases in can we store a guide RNA, so check and store in list
                if (w+(gRNA_length+3)) <= (len(seq)):
                    #add check here for N
                    
                    forward = seq[w+3: w+(gRNA_length+3)]
                    #print forward
                    CODE={'A':'T','T':'A','C':'G','G':'C','N':'N'} 
                    minus_seq=''
                    for c in forward:
                        minus_seq=minus_seq+CODE[c]
                        reverse_seq=minus_seq[::-1]
                    
                    #print reverse_seq
                    if reverse_seq[:1] == "G":
                        guideRnaOutput.append(reverse_seq)
                        guideRna.append(forward)
                        #guideRnaPam.append(seq[w:w+3] + " " + seq[w+3:w+23])
                        start = w+3+inStart
                        end = w+(gRNA_length+3)+inStart -1
                        location.append(str(start) + '-' + str(end))
                        direction.append('rev')
        #increment counter
        w = w + 1


    print guideRna
    #print guideRnaPam
    #print guideRnaOutput
    #print location
    #print direction
    b = 0
    while b < len(guideRna):
        if "N" in guideRna[b]:
            del guideRna[b]
            del guideRnaOutput[b]
            del location[b]
            del direction[b]
            b = b + 1
        else:
            b = b +1 
    #write to file
    print len(guideRna)
    return guideRna, guideRnaOutput, location, direction

def main():
    

    if len(sys.argv) < 2:
      sys.stderr.write("usage: %s <chrom> [<start> <end>]\n" % sys.argv[0])
      exit(2)

    inChrom = sys.argv[1]

    if len(sys.argv) > 2:
      inStart = int(sys.argv[2])
      inEnd = int(sys.argv[3])
    else:
      start = ""
      end = ""

    sys.stderr.write("%s %d %d\n" % (inChrom, inStart, inEnd))
    inStrand = "+"
    fastainput = open("Input_GATA3.fasta", 'w')
    seq_h5 = tables.openFile("/iblm/netapp/data1/external/GRC37/GRC37.h5", "r")
    
    inputRegion = get_seq(inChrom, inStart, inEnd, inStrand, seq_h5)
    #print inputRegion
    fastainput.write(">" + inChrom + '\n')
    fastainput.write(inputRegion)
    fastainput.close()
    f = open("Input_GATA3.fasta", "r")
    #declare chr start and end for region that you want to extract guides from
    
    file = f.readlines()

    
    #outfile = open('GuideRNAs5','w')
    #outfile2 = open('FQGuides5.fq' , 'w')
    outfile3 = open('Database_GATA3.txt', 'w')


    outfile3.write('Location' + '\t' + '\t' + '\t' + '\t' + 'gRNA' + '\t' + '\t' + '\t' + 'Direction' + '\t' + '#other perfect aligns' + '\t' +  'Specificity Score' + '\t' + '%GC Content' + '\t' + '#MM1'+ '\t' + '#MM2'+ '\t' + '#MM3' +  '\n')

    #can change seq and header to lists if there are more than one header in fasta file with several sequences
    #declare empty lists
    sequence = []
    header = []
 
    seq = ""
    header = ""
    GC = []

    #store all the headers in a list
    for f in file:

        if f.startswith('>'):
            header = header + f
            header = header.splitlines()[0]
        #get ride of new line charaters and spaces
        else:
            f = f.replace(" ", "")
            f = f.replace("\n", "")
            seq = seq + f

    #i = 0
    #make it all upper case, easier to parse
    seq = seq.upper()
    
    #print guideRna
    #call function to find all guide RNAs
    #guideRna, guideRnaOutput, location, direction = finding_guides(seq,inStart, inEnd)
    #need to add chromosome start and end position to file

    gRNA_length = 20
    guideRna1, guideRnaOutput1, location1, direction1 = finding_guides(seq,inStart, inEnd,gRNA_length)
    gRNA_length = 19
    guideRna2,guideRnaOutput2, location2, direction2 = finding_guides(seq,inStart, inEnd,gRNA_length)
    gRNA_length = 21
    guideRna3, guideRnaOutput3, location3, direction3 = finding_guides(seq,inStart, inEnd,gRNA_length)

    guideRna = guideRna1 + guideRna2 + guideRna3
    guideRnaOutput = guideRnaOutput1 + guideRnaOutput2 + guideRnaOutput3
    location = location1 + location2 + location3
    direction = direction1 + direction2 + direction3
    j = 0
    #outfile.write(header)
    while j<len(guideRna):

        totalgc = guideRnaOutput[j].count("G") + guideRnaOutput[j].count("C")

        gccontent = (totalgc / float(len(guideRna[j]))) * 100.00
 
        GC.append(gccontent)
        j = j+1


    import subprocess
    f2 = open("aln_GATA3.sam", "r")
    print "Sam file has been opened"
    import re
    
    infoListFull = []
    alnList = []
    q = 0
   
    i = 0
   
   
    
    for line in f2:
        
        #print line
        if line.startswith('@'):
            continue
    
        sys.stderr.write("We are on guide number %d\n" % q)

        alnList = re.findall('(chr[\dMXYIVLR]+),([+-])(\d+),\d+M,(\d)', line)

        #this give list for the qth GUIDE

        infoList1 = alnList
        print len(infoList1)
        #call calculation function

        
        if len(infoList1) != 0:
            #removes all reference sequences that match perfectly somewhere else == 0 mismatches

            numPerfAlns = 0
            mm1 = 0
            mm2 = 0
            mm3 = 0

            for i in range(len(infoList1)):
                if infoList1[i][3] == '0':
                    numPerfAlns = numPerfAlns + 1
                elif infoList1[i][3] == '1':
                    mm1 = mm1 + 1
                elif infoList1[i][3] == '2':
                    mm2 = mm2 + 1
                elif infoList1[i][3] == '3':
                    mm3 = mm3 + 1
                else:
                    # should not see these, but you never know with BWA
                    mm3 = mm3 + 1

            
            if numPerfAlns > 0:
                outfile3.write(header + ':' +  location[q] + '\t' + '\t' + guideRnaOutput[q] + '\t'  + direction[q] + '\t' + '\t' +  str(numPerfAlns) + '\t' + '\t' + '\t' + '1' + '\t' + '\t' + '\t' + str(GC[q]) + '\t' +  '\t' + str(mm1) + '\t' + str(mm2) + '\t' + str(mm3) + '\n')
                q = q +1
                #f2.next()
                continue


            #these thresholds will be changed using graphs from excel that show better values
            #if any of these conditions are true set specificity score to a low value approx 1.
            #Don't do any calculation and move onto next guide
            if mm1 > 5 or mm2 > 40 or mm3 > 500:
                outfile3.write(header + ':' +  location[q] + '\t' + '\t' + guideRnaOutput[q] + '\t'  + direction[q] + '\t' + '\t' +  str(numPerfAlns) + '\t' + '\t' + '\t' + '1' + '\t' + '\t' + '\t' + str(GC[q]) + '\t' +  '\t' + str(mm1) + '\t' + str(mm2) + '\t' + str(mm3) + '\n')
                q = q +1
                continue

            listIndex = 0
            
            testguide = guideRna[q]
            
            SGUIDE = specificity_calc(listIndex,infoList1,testguide,seq_h5)

            outfile3.write(header + ':' +  location[q] + '\t' + '\t' + guideRnaOutput[q] + '\t'  + direction[q] + '\t' + '\t' +  str(numPerfAlns) + '\t' + '\t' + '\t' + str(SGUIDE) + '\t' + '\t' + str(GC[q]) + '\t' +  '\t' + str(mm1) + '\t' + str(mm2) + '\t' + str(mm3) + '\n')
            q = q +1
    

        elif len(infoList1) == 0:
            
            outfile3.write(header + ':' +  location[q] + '\t' + '\t' + guideRnaOutput[q] + '\t'  + direction[q] + '\t' + '\t' +  '-' + '\t' + '\t' + '\t' + '1' + '\t' + '\t'+ '\t' + str(GC[q]) + '\t' +  '\t' + 'too many' + '\t' + 'too many' + '\t' + 'too many' + '\n')
            q = q +1
            
           
        #q = q +1
      
    print "End of program"
    seq_h5.close()


        
    
main()
