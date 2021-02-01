import math
from PIL import Image, ImageDraw, ImageFont


#This method retrieves all of the coding sequences from the genomes as well as the names of the genes
def getIntervals(filename,intervals,names):
    file1=open(filename,"r")
    for line in file1.readlines():
        if line.find('CDS') != -1:      #find a coding sequence,parse the numbers, and then add the interval to the list
            arr=line.split()
            if 'complement' in arr[1]:
                start=(int)(arr[1][arr[1].index("(")+1:arr[1].index(".")])
                end=(int)(arr[1][arr[1].index(".")+2:-1])
                intervals.append([start,end,0])         #the third coordinate indicates the offset layer
            else:
                start=(int)(arr[1][0:arr[1].index(".")])
                end=(int)(arr[1][arr[1].index(".")+2:])
                intervals.append([start,end,0])
        if line.find('/gene') !=-1:             #if a gene name is found, parse it
            add=line[line.index("\"")+1:-2]
            if add not in names:
                names.append(add)

def is_overlapping(x1,x2,y1,y2):            #check if two intervals are overlapping
    return max(x1,y1) <= min(x2,y2)

intervals=[]            #numeric intervals for coding sequences
names=[]                #names of genes
colors=[                                #potential list of colors for coding sequences, can add more if needed
    "red","blue","green","orange","purple","yellow"
]
nameofgenome=""             #organism name
size=0                      #size of genome in bases
w, h = 800, 800             #size of image

#these two variables control location and size of the circle
shape1=300
shape2=500


shape=[(shape1,shape1),(shape2,shape2)]
radius=(shape2-shape1)/2
center=(shape1+shape2)/2        #center coordinate of circle
# creating new Image object
img = Image.new("RGB", (w, h),color="black")

# create rectangle image
img1 = ImageDraw.Draw(img)

#draw base circle for genome
img1.arc(shape, start=270, end=630, fill="white")
img1.arc((shape1+1,shape1+1,shape2-1,shape2-1),start=270, end=630, fill="white")
img1.arc((shape1-1,shape1-1,shape2+1,shape2+1),start=270, end=630, fill="white")

#evenly divide the circle into 6 parts, put 6 marks around the circle
points=[]
points.append(((int)(math.cos(math.pi/2-math.pi/3)*radius+center),(int)(center-math.sin(math.pi/2-math.pi/3)*radius)))
points.append(((int)(math.cos(math.pi/2-2*math.pi/3)*radius+center),(int)(center-math.sin(math.pi/2-2*math.pi/3)*radius)))
points.append(((int)(math.cos(math.pi/2-4*math.pi/3)*radius+center),(int)(center-math.sin(math.pi/2-4*math.pi/3)*radius)))
points.append(((int)(math.cos(math.pi/2-5*math.pi/3)*radius+center),(int)(center-math.sin(math.pi/2-5*math.pi/3)*radius)))
img1.line([(points[0][0]-(int)(math.cos(math.pi/6)*10),points[0][1]+(int)(math.sin(math.pi/6)*10)),(points[0][0]+(int)(math.cos(math.pi/6)*10),points[0][1]-(int)(math.sin(math.pi/6)*10))],fill="white",width=0)
img1.line([(points[1][0]-(int)(math.cos(11*math.pi/6)*10),points[1][1]+(int)(math.sin(11*math.pi/6)*10)),(points[1][0]+(int)(math.cos(11*math.pi/6)*10),points[1][1]-(int)(math.sin(11*math.pi/6)*10))],fill="white",width=0)
img1.line([(points[2][0]-(int)(math.cos(7*math.pi/6)*10),points[2][1]+(int)(math.sin(7*math.pi/6)*10)),(points[2][0]+(int)(math.cos(7*math.pi/6)*10),points[2][1]-(int)(math.sin(7*math.pi/6)*10))],fill="white",width=0)
img1.line([(points[3][0]-(int)(math.cos(5*math.pi/6)*10),points[3][1]+(int)(math.sin(5*math.pi/6)*10)),(points[3][0]+(int)(math.cos(5*math.pi/6)*10),points[3][1]-(int)(math.sin(5*math.pi/6)*10))],fill="white",width=0)
img1.line([(center,shape1-10),(center,shape1+10)],fill="white",width=0)
img1.line([(center,shape2-10),(center,shape2+10)],fill="white",width=0)

#parse the name of the organism and number of bases
filename="Genome.gb"
file=open(filename,"r")

for line in file.readlines():
    if line.find('ORGANISM') != -1:
        nameofgenome=line[line.index("M")+3:-1]

    if line.find('source') != -1:
        arr=line.split()
        size=(int)(arr[1][arr[1].index(".")+2:])
        break

#font sizes
fnt = ImageFont.truetype('/Library/Fonts/Arial.ttf', 24)
fntsmall= ImageFont.truetype('/Library/Fonts/Arial.ttf', 16)
fntxsmall=ImageFont.truetype('/Library/Fonts/Arial.ttf',12)

img1.text((center-100,20),nameofgenome, font=fnt,fill="white")
img1.text((center-100,60),"Number of bases: "+(str)(size),font=fnt,fill="white")

#put size/6 at each of the 6 markers to give a sense of scale to the map
img1.text((center-5,shape1+20),(str)(0),font=fntsmall,fill="white")
img1.text((points[0][0]-(int)(math.cos(math.pi/6)*15)-20,points[0][1]+(int)(math.sin(math.pi/6)*15)),(str)((int)(size/6)),font=fntsmall,fill="white")
img1.text((points[1][0]-(int)(math.cos(11*math.pi/6)*15)-20,points[1][1]+(int)(math.sin(11*math.pi/6)*15)-20),(str)((int)(size/3)),font=fntsmall,fill="white")
img1.text((center-20,shape2-30),(str)((int)(size/2)),font=fntsmall,fill="white")
img1.text((points[2][0]-(int)(math.cos(7*math.pi/6)*15),points[2][1]+(int)(math.sin(7*math.pi/6)*15)-20),(str)((int)(2*size/3)),font=fntsmall,fill="white")
img1.text((points[3][0]-(int)(math.cos(5*math.pi/6)*15),points[3][1]+(int)(math.sin(5*math.pi/6)*15)),(str)((int)(5*size/6)),font=fntsmall,fill="white")

#get the coding sequences and gene names
getIntervals(filename,intervals,names)

#drawing the coding sequences on the map
for i in range(0,len(intervals)):
    
    for k in range(0,i):
        if intervals[i][2]==intervals[k][2]:            #check which offset layer has room for the current coding sequence
            if is_overlapping(intervals[i][0],intervals[i][1],intervals[k][0],intervals[k][1]):
                intervals[i][2]=intervals[i][2]+1
                intervals[0:i].sort(key=lambda lis: lis[2])     #sort previous intervals based on offset
                #print(intervals)

    #calculate start angle by getting start of coding sequence and seeing where it is relative to the whole genome
    startangle=(int)((intervals[i][0]/size)*360)
    startangle+=270

    #take the length of the genome and see how it is relative to the size of the whole genome in order to calculate angle difference
    angledif=(int)(((intervals[i][1]-intervals[i][0])/size)*360)

    #plot the coding sequence on the correct layer based on the offset
    x=10*(intervals[i][2]+1)
    box=[(shape1-20-x,shape1-20-x),(shape2+20+x,shape2+20+x)]
    img1.arc(box,start=startangle,end=startangle+angledif,fill=colors[i])

    #make an entry on the key for the current gene
    img1.rectangle([(720,20+i*40),(740,40+i*40)],fill=colors[i],outline=colors[i])
    img1.text((745,30+i*40),names[i],font=fntxsmall,fill="white")
#display the image
img.show()
