# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 03:34:12 2021

@author: User
"""

from Bio import SeqIO
import sys

inputMSA = sys.argv[1] 

try:  #check threshold is valid or not
    threshold = float(sys.argv[2])
 
except ValueError:  
    print("Your threshold should be float number under 1(100%)!")
    sys.exit()
    
try: #check file is valid or not
    seqFile = open(inputMSA, "r")   #parsing the file
    seqIOOject = SeqIO.parse(seqFile, "fasta")
    
except IOError:   #check file exist
    print("Could not open the file, please check your spelling!")
    sys.exit()
    
except ValueError:  #check file is in fasta
    print("The file should be fasta file!")
    sys.exit()
    

SampleList = [] #this will store all the file information respect to the position



count1 = 0
for record in seqIOOject: #get sequence oly
    seq = record.seq
    
    count2 = 0
    for nucleotide in seq: #this get each position data from the first seqeunces, it helps create a list till all the number of seqeunces are identified  
        if count1 == 0:
            SampleList.insert(count2, [nucleotide])
        
        else: #after seqeunces 1 (index0) it will just append the nucleotide to the respective position 
            SampleList[count2].append(nucleotide)
        
        count2 += 1
        
    count1 += 1


         
numberList = [] #this count ATCG number in list, so index 0 is A; 1 is T; 2 is C; 3 is G
count3 = 0

for value in SampleList:
  
    numberList.append([0,0,0,0]) #start with all 0
     
    for ATCG in value: #based on the ATCG, each number will use line 55 concept to store data
        
        if ATCG == "A":
           numberList[count3][0] += 1
        if ATCG == "T":
           numberList[count3][1] += 1
        if ATCG == "C":
           numberList[count3][2] += 1
        if ATCG == "G":
           numberList[count3][3] += 1
        
    count3 += 1
    

sumList = [] #this count each position valid data, since they are not always the same, it will affect the denominator

for total in numberList:
   sumList.append(sum(total))
   
possibilityList = [] #this count ATCG apparence possibility for each position

count4 = 0
for value in numberList:
    possibilityList.append([0,0,0,0]) #create blank possibility first 
    

    count5 = 0
    for number in value:
        try: #try to caculate the possiblity 
            possibilityList[count4][count5] = number/sumList[count4]
    
        except: #else will just consider as 0%
            
            possibilityList[count4][count5] = 0
        
        count5 += 1
    
    count4 += 1

resultpost = [] #this will store the position, nucleotide, possibility if they are passed the threeshold


count7 = 0
for possibility in possibilityList:
    
    count8 = 0
    for eachone in possibility:
        if eachone <= (1-threshold) and threshold <= eachone: #if this requirement is met, it will record it. Since it will affect both maximum and minimum value of the possibility, it requires in between the value
           Resultlist = count7
           
           resultpost.append([Resultlist,count8, eachone])
         
        count8 += 1
      
    count7 += 1
    
lastposition = 0
lastpercent =  0

count9 = 0
for details in resultpost:
   
    if lastposition == details[0]: #since they will be duplicate position; now I will just make the major allele in front of the minor allele
        
        if details[2] > lastpercent: 
           newresult = resultpost[count9] #if not store previously data, when overwrite it, it will just two duplicate data
           resultpost[count9] = resultpost[count9-1] #switch two postion
           resultpost[count9-1] = newresult
    
    lastposition = details[0]
    lastpercent = details[2]
    count9 += 1
 
lastposition = 0
identify = [] #this store the nucleotide (ATCG), this can use for identify Transversion or  Transition


print("|   Position  | Major allele(Freqeuncy)  |  Minor allele(Freqeuncy) |  Polymorphism type ") #table format
print("===========================================================================================")
for details in resultpost: #since now it is still in number form; switch to ATCG
    if details[1] == 0:
       details[1] = "A"
    if details[1] == 1:
       details[1] = "T"
    if details[1] == 2:
       details[1] = "C"
    if details[1] == 3:
       details[1] = "G"
     
    identify.append(details[1])
    
    if  lastposition != details[0]: #this is use for duplicate position handing
        print("\t",details[0]+1," \t\t\t",details[1],":",(round(details[2]*100)),"%",end = "") #print Position | Major allele(Freqeuncy) data
        #this number add 1 is because index 0 is postion 1

    else: #if it read to the duplicated position (the minor allele one), it will identify Polymorphism type and continue fill up the table 
        if "G" in identify and "A"in identify:
            identify = ["Transition"]
            
        elif "C" in identify and "T"in identify:
            identify = ["Transition"]
        
        else:
            identify = ["Transversion"]
        
        
        print(" \t\t\t\t\t ",details[1],":",round(details[2]*100),"%","\t\t\t\t",identify[0]) #differnet than previous is that, after this line finish, it will just print next line (since they will only be one major and minor allele)
        identify = [] #if not blank it, it will store infinite ATCG
    
    lastposition = details[0] #if not blank it, it will store infinite position that use for checking duplicate position happening

print("===========================================================================================")


    
    