"""
Updated on 2022-04-02

@author: Dakforz Chang 
"""


from Bio import SeqIO
import sys #since this program use spyder, the argument is stored in run-configuration per file function

try:
 seqFile = open(sys.argv[1],"r")  #try for first argument - sequences.fastq
 seqIOObject = SeqIO.parse(seqFile, sys.argv[3]) #try for third argument - fastq 
 
except IOError:
 print("Could not open the file!") #confirm first argument is exist
 sys.exit()

except ValueError:
 print("The file should be fasta or fastq!") #if the file type is wrong, it will pop out 
 sys.exit()
 
try: 
 kmer = int(sys.argv[2]) #try second argument -17 here, if it is not a integar number the error will show
except ValueError:
 print("The kmers should be an integar number!")
 sys.exit()

kmers = [] #this store every kmer sequence (ex.ATCG TCGA..)
kmersnum = [] #this store every sequence repeat number (they are follow by kmer's index)

Filehandle= open("result.txt","w") #write the following graph in result file
space = " " #help with formating 
Filehandle.writelines("%s %s %s\n" % ("sequence(s)",space*(kmer-4),"count(s)")) #this is the title of the graph;space*(kmer-4) will locate in the center  

for record in seqIOObject: #ths read each line of sequences.fastq and then anaylze each line
    seq =record.seq #only get the sequence information
    n= 0 #ths is used for tracking kmer's number (n-k+1)
  
    while kmer+n <= int(len(seq)):  #this make sure each kmers do not exceed its assign number
        newstring = seq[0+n:kmer+n] #this specify the sequences's kmer using the range
        #print(kmer+n, len(seq),newstring)
        if newstring in kmers:  #this check for kmer replication
            location = int(kmers.index(newstring)) #this check the replication index location
            existingnum = kmersnum[location] #this read the existing replication number (since it may be larger than 1)
            newnum = int(existingnum) + 1 #the kmer number will add value 1 
            kmersnum[location] = newnum  #store/replace the above number into list 
            
        else: #if kmer do not replicate(exist before) it will write the kmer sequence in list, and assign value 1 to another list (first appearance)
            kmers.append(newstring)
            kmersnum.append('1')
        n += 1 #this will make sure our kmer push left
    
       
for writein in range(len(kmers)): #this line just write the list into file 
     Filehandle.writelines("%s %s %s\n" % (kmers[writein],space*(kmer-8),kmersnum[writein]))


seqFile.close() #close the references sequence file 
Filehandle.close() #close the graph file
       
  
    
