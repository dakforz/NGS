# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 10:47:03 2021

@author: User
"""

from Bio import SeqIO
import random
import sys

kmersdict = {} #this is the location that your kmer count file will be 

kmercountfile =sys.argv[1] 
readstocorrect =sys.argv[2]

try:
   kmersize = int(sys.argv[3]) #error handling to make sure number is valid
except ValueError:  
    print("Your kmersize should be an integar number!")
    sys.exit()
    
try:   
    threshold = int(sys.argv[4]) #error handling to make sure number is valid
except ValueError: 
    print("Your threshold should be an integar number!")
    sys.exit()

allseqlist = []
ATCG =['A','T','C','G'] #this is for fixing error, it will combine with random function in line 9
try:
    seqFile = open(sys.argv[1]) #read the kmer file generated by jellyfish
    seqIOObject = SeqIO.parse(seqFile,'fasta')
    
except IOError:   #check file exist
    print("Could not open the kmer count file, please check your spelling!") 
    sys.exit()
    
except ValueError:  #check file is in fasta
    print("The kmer count file should be fasta file!")
    sys.exit()
    
    
for kmer in seqIOObject: #this is for parsing and reading the kmers' seqenuce only
   seq = kmer.seq
   kmersdict[str(seq)] =int(kmer.description)
seqFile.close()

try:
    seqFile2 = open(sys.argv[2]) #read the file you want to correct
    seqIOObject2 = SeqIO.parse(seqFile2,'fastq')
    
except IOError:   #check file exist
    print("Could not open the file that you want to correct, please check your spelling!")
    sys.exit()
    
except ValueError:  #check file is in fasta
    print("The file that you want to correct should be fastq file!")
    sys.exit()
    

for record in seqIOObject2: #this read all the seqeunces and store them in a list
    seq2 = record.seq 
    allseqlist.append(str(seq2))
    
seqFile.close()

countseq = 0
checklist = ["0"] #this is used for confirm your sequence is correct or not, get 0 (False) is because now I used while true

for seq in allseqlist: 
    count = 0
    for i in range(len(seq)-1): #this is the first loop
        start = count - kmersize +1 #this give the kmer range start and end from that nucleotide (position)
        end = start + kmersize - 1
        
    
            
        count1 = 0
        
        while True:
            
          
            if "1"  in checklist or "0" not in checklist: #since if there is one kmer has high freqeuncy, we can say that is correct, 1 will generate when kmer coverage pass the threshold
                checklist = ["0"] #reset to zero, so next postion can keep running
                break
            
            else:
                checklist = [] #since we want a blank checklist now, it will just in this loop till the end 
                count1 = 0
                
                for x in range(start, end+1):
            
                    kmertest = allseqlist[countseq][start+count1 : start+kmersize+ count1] #this will be your kmer 
                  
                    if len(kmertest) == kmersize: #this line exist becaus first or following and last or previous postion have less kmers that include in that positon, if do not have this line, they will be no seqeunces or less sequence
         
                        if  kmertest not in kmersdict or kmersdict[kmertest] <= threshold: #check the kmer has passed the threshold or not
                            checklist.append("0") #if not append 0 (False)
                        else:
                            checklist.append("1")#if yes append 1 (True)
                       
                    else: #those not valid kmer just get 0 instantly, saved some space based on line 96
                        checklist.append("0")
                        
                    count1 += 1     #this will end until one nucleotide (postion) is done  
                    
                    
                if "1"  not in checklist and "0" in checklist:   #this step will just help fix the kmer, and go again while loop(line 81) until it success, next positon will then move (line 72)
                    allseqlist[countseq] = allseqlist[countseq][0:count] + random.choice(ATCG) + allseqlist[countseq][count+1:len(seq)] #even though random may have generated same character, it is just elimate lots of condition required to set 
                    checklist = ["0"]
                 
                           
                continue                        

        count += 1
    countseq +=1
    
outfile = open("outfile.fastq","w")    #write the output

seqFile2 = open("test_data.fastq","r") #open the orginal file that you want to correct

count2 = 3 #this is based on fastq format, i want seqeunces are in 4* so set it initially as 3
count3 = 1
for record in seqFile2:
    if count2/(4*count3) != 1: #this will only true if you are at seqeunces 
        outfile.write(record) #copy the information of your fastq format, except sequences
    else:
        outfile.write(allseqlist[count3-1]+"\n")  #since we store the sequences in list in order, we can just based on order copy to the file 
        count3 += 1
    count2 += 1
    
seqFile2.close() 
outfile.close() #close all the file


