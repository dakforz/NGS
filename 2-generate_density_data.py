"""
Created on Thu Oct 13 09:58:24 2021

@author: Dakforz Chang 400076414
"""
import matplotlib.pyplot as plt #for graph
import sys #for command line argument


vcfFile = sys.argv[1]

try:   #confirm it is a vcf file
    openFile = open(vcfFile,"r")    
    
except IOError:
    print("Could not open the file!")
    sys.exit()
except ValueError:
    print("The file should be vcf file!")
    sys.exit()
    
position = []
sample = []
totalpos = 0



try:  #confirm window size is integar number
  windowsize = int(sys.argv[2]) 
  
except ValueError:
  print("The window size should be an positive integar number!")
  sys.exit()
  
try:   #confirm increment value is integar number 
  incremement = int(sys.argv[3])
  incremement > 0
  
except ValueError:
  print("The incremement value should be an positive integar number!")
  sys.exit()
  
  

for record in openFile:  #this read the vcf file
    
    if '##' not in record: #this is used to filter out the desciption line, I keep #, if require sample name, it can be useful
    
        for line in record: #this read each line
            data = record.split('\t')  #this separate the line by tab
   
        count = 0 
    
        polydatainrow = [] #this store the polymorphism information
    
    
        for line in data: 

            if count == 1: #this is the position location
                position.append(line) #collect all the position information in "postion"
            
            if count >= 9: #this is the sample data (it can collect all the sample)
                if "1/1" in line: 
                    polydatainrow.append(1) #1/1 indicate polymorphism so +1
                else:
                    polydatainrow.append(0) #0/0 is the only other case, it indicate no polymorphism so +0
        
        
            count += 1
        
        sample.append(polydatainrow)  #(all sample at once) ex. [ [position1][position2]....positionlast];each position can occupy mutiple sample [position1]= [sample1,sample2]
        totalpos += 1 
        
    
del position[0] #right know the first information will be "POS" not necessary so delete (since we keep #)
position = [int(x) for x in position]
del sample[0] #this is the sample name, right know we don't use it, so delete it

#till this point all the vcf file information that required are stored in python

totalsample = count-9 + 1 #this is the total sample amount


xpos = 0 
ypos = 0

x = [0] #the first number of it will always be 0
y = [] #y will not be 0, since we have windowsize
total = []


#x axis number can be known simply by having incremement value and last position of DNA 
while int(x[xpos]) < position[-1]:
    x.append(incremement*(xpos+1)) #if added incremement value is not bigger than last position just keep add it and store in list 
    xpos += 1
    


#this is the step to find the position of y when x = 0

while position[ypos] < windowsize: #this will check if sample postion is less than windowsize
    count = 0
    
    for value in sample[ypos]: #if aboe is the case, t will store the polymorphism information respect to that position 
        
        if ypos == 0: #for the first sample, it requires to establish the index postion
            total.insert(count,[int(value)])
            count += 1
            
        else:
            total[count].append(int(value)) #other than first one, other can simply add the value
            count += 1

    
    ypos += 1
#till this point, all the sample will be stored in same list no longer separate by position

for value in range(totalsample-1):
    total[value] = sum(total[value]) #sum the number to get the y when x = 0; this is only one y value now
    


y.append(total) #make the above number in y, so it is clear
#print(y)



#now, below we want to do one thing to get y + 1*n, add the range that is going to cover further and substract the range that is no longer cover 
#for instance, 1000 windowsize, 10 incremement. x=0 y=polymorphism in range 0-1000. Go to x=10, we want 10-1010. So add range of polymorphism number between 1000-1010 and delete 0-10 range's polymorphism number

newy = y[0].copy() #without copy they will share the same list, newy will be treated to +/- and always be the recent y in further step 

initalcount = 0
lastcount = ypos #this is the last number of windowsize end (above)
#print(y[0])



incrementcount = 1 

for i in range(len(x)-2): #-2 is because we start at 1 and do not need last 
  
     windowsize += incremement #the whole process will stop until we reach x last positon
     
     deltotal = [] #the substrate part, it will be blank after each round
     
     
     while int(position[initalcount]) < incremement*incrementcount: #first find the part we want to substract
    
         deltotal.append(sample[initalcount].copy()) #if it is smaller than incremement*incrementcount (the increment count indicate shift to right)
         initalcount += 1
        
     if deltotal == []: #to prevent no polymorphism in the range and store nothing, we add blank 0 list into it
         deltotal = [[0 for x in range(totalsample-1)]]
            
        
         
     for value in deltotal:
         count = 0
             
         for num in value: #this two line badsically are making substract part work to each sample 
             newy[count] = (int(newy[count])-int(num))
             count += 1
        
         
                 
           
     incrementcount += 1
     
     addtotal = []   #same as  deltotal concept
     
     if  lastcount < (totalpos-2): #same as above, we want to find the range we want to add 

     
         #print(lastcount,windowsize)  
         
         while int(position[lastcount]) < windowsize:
               addtotal.append(sample[lastcount].copy())
               lastcount += 1  
               
               if lastcount == (totalpos-2): #this avoid that lastcount bigger than existing value and cause error, since it might happen due to in further step, we will only substract value not add value
                   y.append(newy.copy())
                   break
           
        
          
         if addtotal == []: #same as line 153
             addtotal = [[0 for x in range(totalsample-1)]]
             
         for value in addtotal: #same as line 158
             count = 0
             
             for num in value:
                 newy[count] = (newy[count]+int(num))
                 count += 1
                    
           
         y.append(newy.copy()) #this is just add the final newy.copy() in it, since newy is never reset, it will be treated to become next y 
             #print(newy)
             
     else:
         y.append(newy.copy()) #concept same as 181

        
     

   


#this definely can be treated with range function and line 77 to automatically generated the name or so 
#However, getting the color, and compared interested sample will be more convience with those small group

I = []
C = [] 
B = []
A = []
P = []
D = []	
F = []
E = []

for sample in y: #this is the step that organized the sample to different y list 
    
    I.append(sample[0])
    C.append(sample[1])
    B.append(sample[2])
    A.append(sample[3])
    P.append(sample[4])	 
    D.append(sample[5])
    F.append(sample[6])
    E.append(sample[7])
    
    
plt.plot(x, I, label = "I") #plot the graph by sample
plt.plot(x, C, label = "C") 
plt.plot(x, B, label = "B")
plt.plot(x, A, label = "A")
plt.plot(x, P, label = "P")
plt.plot(x, D, label = "D")
plt.plot(x, F, label = "E")
plt.plot(x, E, label = "F")

plt.xlabel('Chromosome Position')    
plt.ylabel('Number of Polymorphism')
plt.title('Chromosome 2 Polymorphism Density Plot')
plt.yscale("log") #y with 10* scale
plt.legend() #top right corner information
plt.show()
    
 
