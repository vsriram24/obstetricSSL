# Author: Vivek Sriram
# Last revised: 08/18/2021
# ------------------------------------------------
# Function: Given a directory of comorbidity data,
#			create a corresponding comorbidity map

# Import statements
import csv
import argparse
import ast
import os
import sys
import glob

inputFiles = "/Users/viveksrm/Desktop/ukbbFemaleComorbidityData/*"
allComorbiditiesPath = "/Users/viveksrm/Desktop/ukbbFemaleComorbidities.csv"

#comorbidityMap = set()
#comorbidityMissingLinks = set()
comorbiditiesAll = set()

filesToBeParsed = glob.glob(inputFiles)
if(len(filesToBeParsed) == 0):
	sys.exit('ERROR: No input files have been provided. Check the accuracy of file names for the --input-files argument')

i = 1
for name in filesToBeParsed:
	print("On phenotype " + str(i) + " out of " + str(len(filesToBeParsed)))
	with open(name, 'r') as f:
		phenotype1 = name.split("/")[5].split(".csv")[0].split("PheCode_")[1] + "_ld"
		firstLine = next(f)

		fileLines = f.readlines()[1:]
		for line in fileLines:
			lineStripped = line.strip()
			lineAsArray = lineStripped.split(",")

			phenotype2 = lineAsArray[1].split("PheCode_")[1]
			preDecimalOfPhenotype2 = phenotype2.split(".")[0]
			if(len(preDecimalOfPhenotype2) == 1):
				phenotype2 = "00" + phenotype2
			elif(len(preDecimalOfPhenotype2) == 2):
				phenotype2 = "0" + phenotype2
			phenotype2 = phenotype2 + "_ld"

			coOccur = lineAsArray[3]
			phiValue = lineAsArray[5]
			pValue = lineAsArray[6]
			#print(phenotype1)
			#print(phenotype2)
			#print(pValue)

			if(pValue != ""):
				currentEdge = [phenotype1, phenotype2, coOccur, phiValue, pValue]

				comorbiditiesAll.add(tuple(currentEdge))

				#if(float(pValue) < float(0.05) and float(coOccur) > 1 and float(phiValue) > 0):
				#	comorbidityMap.add(tuple(currentEdge))
				#else:
				#	comorbidityMissingLinks.add(tuple(currentEdge))
	
	i = i+1

print("Finished identifying edge map and missing link map")


#with open(edgeMapOutputPath, 'w') as outfile:
#	outfile.write("phenotype1\tphenotype2\tcoOccur\tphiValue\tpValue\n")
#	for elem in comorbidityMap:
#		outfile.write(str(elem[0]) + "\t" + str(elem[1]) + "\t" + str(elem[2]) + "\t" + str(elem[3]) + "\t" + str(elem[4]) + "\n")
#print("Finished outputting " + str(len(comorbidityMap)) + " edges to the edge map.")


#with open(missingLinkMapOutputPath, 'w') as outfile2:
#	outfile2.write("phenotype1\tphenotype2\tcoOccur\tphiValue\tpValue\n")
#	for elem in comorbidityMissingLinks:
#		outfile2.write(str(elem[0]) + "\t" + str(elem[1]) + "\t" + str(elem[2]) + "\t" + str(elem[3]) + "\t" + str(elem[4]) + "\n")
#print("Finished outputting " + str(len(comorbidityMissingLinks)) + " edges to the missing link map.")


with open(allComorbiditiesPath, 'w') as outfile3:
	outfile3.write("phenotype1\tphenotype2\tcoOccur\tphiValue\tpValue\n")
	for elem in comorbiditiesAll:
		outfile3.write(str(elem[0]) + "\t" + str(elem[1]) + "\t" + str(elem[2]) + "\t" + str(elem[3]) + "\t" + str(elem[4]) + "\n")
print("Finished outputting " + str(len(comorbiditiesAll)) + " edges to the all comorbidities map.")



