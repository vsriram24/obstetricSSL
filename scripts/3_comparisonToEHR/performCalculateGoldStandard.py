import subprocess
import pandas as pd

# Read in the GBSSL output file to get a list of phenotypes
gbsslResults = pd.read_csv("/Users/viveksrm/Desktop/pregnancy_gbssl_outcomes/gbsslScores_653.csv", dtype=str)
phenotypesToQuery = gbsslResults['phenotype'].tolist()

# Read in the phecode9Map.csv file
phecodeToICD9File = pd.read_csv("/Users/viveksrm/Documents/GitHub/peGBSSL/data/phecode9Map.csv", delimiter='\t', header=0, dtype=str)
phecodeToICD9Dict = {}
# Read in the phecode10Map.csv file
phecodeToICD10File = pd.read_csv("/Users/viveksrm/Documents/GitHub/peGBSSL/data/phecode10Map.csv", delimiter='\t', header=0, dtype=str)
phecodeToICD10Dict = {}

# Make dictionaries out of phecodeToICD10Map and phecodeToICD9Map
for index, row in phecodeToICD9File.iterrows():
    phecodeToICD9Dict[row.iloc[0]] = row.iloc[1]

for index, row in phecodeToICD10File.iterrows():
	phecodeToICD10Dict[row.iloc[0]] = row.iloc[1]

pvalues = []
i = 1
# For each row in our GBSSL output file, get the phecode
for elem in phenotypesToQuery:
	phenotypeLength = len(elem)
	elem = elem[0:phenotypeLength-3]
	print("On phenotype " + elem + ", " + str(i) + " out of " + str(len(phenotypesToQuery)))
	# Find the corresponding ICD-9 codes from phecode9Map.csv
	phecode9List = phecodeToICD9Dict.get(elem, '')
	if(phecode9List != ''):
		phecode9List = phecode9List[1:len(phecode9List)-1]

	# Find the corresponding ICD-10 codes from phecode10Map.csv
	phecode10List = phecodeToICD10Dict.get(elem, '')
	if(phecode10List != ''):
		phecode10List = phecode10List[1:len(phecode10List)-1]

	# ['<command_to_run>', '<path_to_script>', 'arg1' , 'arg2', 'arg3', 'arg4']
	# Input = ukbbEHR csv file, phecode0list, phecode10list
	# Use our R script to get the p-value
	cmd = ['/usr/local/bin/RScript', '/Users/viveksrm/Documents/GitHub/peGBSSL/scripts/calcPValue_HDP.R', '/Users/viveksrm/Desktop/ukbbEHRFeather/ukbbEHR_653.feather', phecode9List, phecode10List]
	out = subprocess.check_output(cmd)
	outAsString = out.decode("utf-8")
	outAsArray = outAsString.split(" ")
	finalElem = outAsArray[len(outAsArray)-1]
	print(finalElem)
	print(len(finalElem))
	if(len(finalElem) == 1):
		print("No p-value from test")
		currentPValue = 1
	#currentPValue = float(out)
	else:
		print("Calculating pvalue")
		currentPValue = float(finalElem)
	print("Got a p-value of " + str(currentPValue))
	pvalues.append(currentPValue)
	i = i+1

gbsslResults['ehrSignificance'] = pvalues
gbsslResults.to_csv("/Users/viveksrm/Desktop/coAssociationResults_653.csv")

