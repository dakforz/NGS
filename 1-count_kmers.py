# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 23:50:48 2022

@author: dakforz

This program will convert sequences into kmers

Details: 
1.read any fasta or fastq file
2.break all the sequences into your desired kmer size 
3.count how many times each kmer appears in the file of sequencing read
4.create a table into result.txt that have all the kmer seqeunce and it's number of appearance

#sys.argv[1] input with filename (ex.sequence.fastq)
#sys.argv[2] input with interested kmer size (ex.17)
#sys.argv[3] input with filetype (fastq)



"""

from Bio import SeqIO
import sys

try:  #error handling for file's name and type 
    readFile = open(sys.argv[1],"r")   
    seqIOObject = SeqIO.parse(readFile, sys.argv[3])
    
except IOError:
    print("Can not open the file, Please re-check your spelling.")
    sys.exit()

except ValueError:
    print("The filetype should be in .fasta or .fastq.") 
    sys.exit()
 
try:  #error handling for kmer size's number
    kmer = int(sys.argv[2]) 
except ValueError:
    print("The kmer size should be an integar number.")
    sys.exit()



#kmers and kmers_counts shared the same index
kmers = []  #this will only store the kmer information
kmers_counts = [] #this will store the counts information

writeFile= open("result.txt","w")  #write the results into result.txt file 
space = " "  #this is for formating the file only
writeFile.writelines("%s %s %s\n" % ("sequence(s)",space*(kmer-4),"count(s)")) 

for record in seqIOObject:
    seq =record.seq #only getting sequence information
    n= 0  #this is for tracking kmer's number (n-k+1)
  
    while kmer+n <= int(len(seq)):  
        newstring = seq[0+n:kmer+n]
        
        if newstring in kmers:  #if the kmer has already existed, it will increase the counts
            location = int(kmers.index(newstring))
            existing_counts = kmers_counts[location] 
            new_counts = int(existing_counts) + 1 
            kmers_counts[location] = new_counts  
            
        else:
            kmers.append(newstring)
            kmers_counts.append('1')
        n += 1 
    
       
for index in range(len(kmers)): 
     writeFile.writelines("%s\t %s\n" % (kmers[index],kmers_counts[index]))


readFile.close() 
writeFile.close() 
