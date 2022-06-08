"""
Created on Sun May  1 00:20:54 2022

@author: dakfo

This program will create polymorphism density plot based on vcf file

Details: 
1.read vcf file
2.separate the vcf file into chromosome position and its' sample polymorphism
2.Count the number of polymorphisms based on input windowsize and increment value 
3.create a polymorphism density plot 

#sys.argv[1] input with filename (ex.chr02.vcf)
#sys.argv[2] input with interested windowsize (ex.1000000)
#sys.argv[3] input with increment value (100000)

"""

import matplotlib.pyplot as plt 
import sys 


vcfFile = sys.argv[1]

try:   #error handling for file's name and type 
    readFile = open(vcfFile,"r")  
      
except IOError:
    print("Can not open the file, Please re-check your spelling.")
    sys.exit()
    
except ValueError:
    print("The file type should be in .vcf.")
    sys.exit()
    
    
try:   #error handling for windowsize
  windowsize = int(sys.argv[2]) 
  
except ValueError:
  print("The window size should be an positive integar number!")
  sys.exit()
  
try:  #error handling for incremement value
  incremement = int(sys.argv[3])
  incremement > 0
  
except ValueError:
  print("The incremement value should be an positive integar number!")
  sys.exit()
  
  
position = [] #this store all the chromosome position information
sample = []  #this store polymorphism data as array 
rows = 0

samplename = []

for record in readFile: 
    
    if '##' not in record: #this line filter out the file desciption line if existed
    
        for line in record: 
            line_data = record.split('\t') #separated each line elements
   
        line_elements_counts = 0 
    
        polymorphism_count_row_samples = [] #this store the polymorphism information
    
    
        for elements in line_data: 

            if line_elements_counts == 1: #this is the position data 
                position.append(elements) 
            
            if line_elements_counts >= 9: #this is the sample data (it can collect all the sample)
                if "1/1" in elements: 
                    polymorphism_count_row_samples.append(1) #1/1 indicate polymorphism existed
                else:
                    if rows == 0:
                        samplename.append(elements)
                    else:
                        polymorphism_count_row_samples.append(0) 
               
        
            line_elements_counts += 1
        
        sample.append(polymorphism_count_row_samples)
 
        rows += 1 
        
del position[0] # columns name 

del sample[0] #first column [will be []]

position = [int(x) for x in position]


#till this point all the vcf file information that required are stored in python

totalsample = line_elements_counts-9 + 1 #this is the total sample amount


x_index = 0  
y_index = 0

x = [0] #x axis; the first number will always be 0
y = [] #y axis
total_polymorphisms_x_0 = [] #all the polymorphism values when x= 0; sum(total_polymorphisms_x_0) = y[0]


#create x axis 
while int(x[x_index]) < position[-1]:
    x.append(incremement*(x_index+1))
    x_index += 1
    


#Find y value when x = 0

while position[y_index] < windowsize: 
   
    count = 0
    
    for value in sample[y_index]: 
        
        if y_index == 0: #first sample require an inital value
            total_polymorphisms_x_0.insert(count,[int(value)])
            count += 1
            
        else:
            total_polymorphisms_x_0[count].append(int(value)) 
            count += 1

    y_index+= 1

for value in range(totalsample-1):
    total_polymorphisms_x_0[value] = sum(total_polymorphisms_x_0[value]) 


y.append(total_polymorphisms_x_0) #first y value; other will based on y[0] for caculation 



#The above function add the range that is going to cover and substract the range that is no longer cover 
#for instance, 1000 windowsize, 10 incremement. x=0 y= (polymorphisms in range 0-1000). When it is x=10, the range will turn to 10-1010. So add range of polymorphism number between 1000-1010 and delete 0-10 range's polymorphism 
handle_y_value = y[0].copy() #without copy they will share the same list, handle_y_value will be treated to +/- and always be the recent y in further step 

initalcount = 0
lastcount = y_index #this is the last count of inital windowsize 


incrementcount = 1 

for i in range(len(x)-2): #-2 is because x start with index = 1 and do not need last x position  
  
     windowsize += incremement #the whole process will stop until x > last position
     
     deltotal = [] #the substrate part, it will be blank after each round
     
     while int(position[initalcount]) < incremement*incrementcount: #substract the value that is no longer cover
    
         deltotal.append(sample[initalcount].copy()) 
         initalcount += 1
        
     if deltotal == []: #to prevent no polymorphism in the range and store nothing, we add blank 0 list into it
         deltotal = [[0 for x in range(totalsample-1)]]
            
        
     for value in deltotal:
         count = 0
             
         for num in value: #this add all the number in the substrate part and store in "handle_y_value"
             handle_y_value[count] = (int(handle_y_value[count])-int(num))
             count += 1
           
     incrementcount += 1
     
     
     addtotal = []  #the add part, it will be blank after each round
     
     if  lastcount < (rows-2): #the concept is same as substrate part 
         
         while int(position[lastcount]) < windowsize:
               addtotal.append(sample[lastcount].copy())
               lastcount += 1  
               
               if lastcount == (rows-2): #this avoid that lastcount can be bigger than existing value and cause error
                   y.append(handle_y_value.copy())
                   break
           
        
         if addtotal == []: 
             addtotal = [[0 for x in range(totalsample-1)]]
             
         for value in addtotal: 
             count = 0
             
             for num in value: 
                 handle_y_value[count] = (handle_y_value[count]+int(num))
                 count += 1
                    
    
         y.append(handle_y_value.copy()) #this is just add copy of each final handle_y_value.copy() in y, since handle_y_value is never reset, it will be treated to gnerate next y value 
             
     else:
         y.append(handle_y_value.copy()) 
     



#this is plotting graph part; until this point, y store all sample polymorphism data 

plt.xlabel("Chromosome Position")
plt.ylabel("Number of Polymorphism")
plt.yscale("log") #y with 10* scale
plt.title("Chromosome 2 Polymorphism Density Plot")

for i in range(len(y[0])):
    plt.plot(x,[pt[i] for pt in y],label = samplename[i])
plt.legend()
plt.show()

