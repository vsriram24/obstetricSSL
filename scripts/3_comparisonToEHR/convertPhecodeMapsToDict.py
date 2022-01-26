# Import statements
import pandas as pd


##############################################
# Read in our mappings
phecodeToICD10_directory = "/Users/viveksrm/Desktop/Phecode_map_v1_2_icd10_beta.csv"
phecodeToICD10 = pd.read_csv(phecodeToICD10_directory, sep = ",", dtype = "str")
phecodeToICD10 = phecodeToICD10.dropna()
numRows_ICD10 = phecodeToICD10['ICD10'].count()

phecodeToICD9_directory = "/Users/viveksrm/Desktop/phecode_icd9_map_unrolled.csv"
phecodeToICD9 = pd.read_csv(phecodeToICD9_directory, sep = ",", dtype = "str")
phecodeToICD9.dropna()

icd10Dict = {}
for ind in phecodeToICD10.index: 
	currentPhecode = phecodeToICD10['PHECODE'][ind]
	currentICD10 = phecodeToICD10['ICD10'][ind]
	currentICD10AsArray = currentICD10.split(".")
	newICD10 = ""
	for elem in currentICD10AsArray:
		newICD10 = newICD10 + elem

	# We want to ignore phenotypes that have two numbers after the decimal place
	currentPhecodeAsArray = currentPhecode.split(".")
	if(len(currentPhecodeAsArray) > 1):
		if(len(currentPhecodeAsArray[1]) < 2):
			if(currentPhecode in icd10Dict):
				icd10Dict[currentPhecode].append(newICD10)
			else:
				icd10Dict[currentPhecode] = [newICD10]

	if(len(currentPhecodeAsArray) == 1):
		if(currentPhecode in icd10Dict):
			icd10Dict[currentPhecode].append(newICD10)
		else:
			icd10Dict[currentPhecode] = [newICD10]
        
icd9Dict = {}
for ind in phecodeToICD9.index: 
	currentPhecode = phecodeToICD9['phecode'][ind]
	currentICD9 = phecodeToICD9['icd9'][ind].strip(".")
	currentICD9AsArray = currentICD9.split(".")
	newICD9 = ""
	for elem in currentICD9AsArray:
		newICD9 = newICD9 + elem

	currentPhecodeAsArray = currentPhecode.split(".")
	if(len(currentPhecodeAsArray) > 1):
		if(len(currentPhecodeAsArray[1]) < 2):
			if(currentPhecode in icd9Dict):
				icd9Dict[currentPhecode].append(newICD9)
			else:
				icd9Dict[currentPhecode] = [newICD9]

	if(len(currentPhecodeAsArray) == 1):
		if(currentPhecode in icd9Dict):
			icd9Dict[currentPhecode].append(newICD9)
		else:
			icd9Dict[currentPhecode] = [newICD9]

with open('/Users/viveksrm/Desktop/phecode10Map.csv', 'w') as f:
    for key in icd10Dict.keys():
        f.write("%s\t%s\n"%(key,icd10Dict[key]))
        
with open('/Users/viveksrm/Desktop/phecode9Map.csv', 'w') as f:
    for key in icd9Dict.keys():
        f.write("%s\t%s\n"%(key,icd9Dict[key]))