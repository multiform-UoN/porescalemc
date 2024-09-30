#!/usr/bin/python
import sys

f=open(sys.argv[1])

regions=0
for line in f:
	if "*Number of regions: " in line:
		regions=line[(line.find('*Number of regions: ')+20):-1]
#print(regions)
if regions==0:
    quit()


print("cellSet saved new cellToCell region1"+"\n")

regs=int(regions)
for i in range (regs-2):
	print("cellSet saved add cellToCell region"+str(i+2)+"\n")

print("cellSet saved invert"+"\n")
print("cellSet saved subset"+"\n")
print("quit")
