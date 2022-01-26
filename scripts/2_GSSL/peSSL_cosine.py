#####################################
# Title: peSSL_cosine.py
# Authors: Vivek Sriram, Yonghyun Nam 
# Last revised: 06/07/21
# -----------------------------------
# Function: Run graph-based SSL on input PheWAS summary data to calculate disease similarity
# 			for a chosen index phenotype(s). Includes cosine similarity normalization of the
#			association matrix

# Import statements
import numpy as np
import pandas as pd
import os
import sklearn
import sklearn.metrics
import statistics
import math



##############################################
# 0. GET OUR INPUT DATA AND DISTANCE ARRAY
# Read in our raw data
rawDataDirectory = "/Users/viveksrm/Documents/GitHub/peGBSSL/data/rawPheWASData/neg4RawData_tenths_noSymptoms/"

# CREATE OUR ASSOCIATION MATRIX
# The rows of the matrix will be the phenotypes, and the columns will be the SNPs
# Initialize a dictionary
phenotypeToSNPDictionary = {}
uniquePhenotypes = []
uniqueSNPs = []

i = 1
numFiles = len(os.listdir(rawDataDirectory))

for filename in os.listdir(rawDataDirectory):
	print("Processing file " + filename + ", " + str(i) + " out of " + str(numFiles))
	currentFile = pd.read_csv(rawDataDirectory + "/" + filename, sep = " ")
	phenotypeName = os.path.splitext(filename)[0]

	snpsForThisPhenotype = 	currentFile['ID'].tolist()
	if(len(snpsForThisPhenotype) > 0):
		uniquePhenotypes.append(phenotypeName)
		for elem in snpsForThisPhenotype:
			if(elem not in uniqueSNPs):
				uniqueSNPs.append(elem)
            
	phenotypeToSNPDictionary[str(phenotypeName)] = snpsForThisPhenotype
	i = i+1

# Now we have a dictionary of phenotypes and their corresponding SNPs
# We also have a list of unique phenotypes and SNPs
# Initialize a dataframe that is the right dimensions: number of rows = number of phenotypes, number of columns = number of SNPs
# Label all the rows and columns according to our uniquePhenotypes and uniqueSNPs lists
associationMatrix = pd.DataFrame(np.zeros((len(uniquePhenotypes), len(uniqueSNPs))),
 					columns=uniqueSNPs, 
 					index=uniquePhenotypes)

for phen in phenotypeToSNPDictionary:
	currentListOfSnps = phenotypeToSNPDictionary[phen]
	for snp in currentListOfSnps:
		associationMatrix[snp][phen] = 1.0


associationMatrix['rs111398676'].sum()
droplist = [i for i in associationMatrix if associationMatrix[i].sum()<=1]
associationMatrix.drop(droplist,axis=1,inplace=True)


associationMatrix.to_hdf("/Users/viveksrm/Desktop/assocMat.h5", key='stage', mode='w')

associationMatrix = pd.read_hdf("/Users/viveksrm/Desktop/assocMat.h5")

W_w = sklearn.metrics.pairwise.cosine_similarity(associationMatrix)
W_u = np.where(W_w > 0, 1, 0)



np.save("/Users/viveksrm/Desktop/W_w.npy", W_w)

DistW = 1-W_w
DistU = 1-W_u


##############################################
# 1. LABEL PROPAGATION FUNCTION
# Label_Propagation function - takes as input the number of nodes,
#	an inverse Laplace, the indices of group 1, and the indices
#	of group -1
def Label_Propagation(N_nodes, inv_Lap, labeled_pos, labeled_neg):
	# Initialize np lists of length equaling the number of nodes
    f = np.zeros(N_nodes)
    y = np.zeros(N_nodes)

    # Set the labels in our list y
    y[labeled_pos] = 1;
    y[labeled_neg] = -1;

    # Multiple the inverse Laplace with our list of labels
    f = inv_Lap@y
    f = f.transpose()
    
    return (f)


##############################################
# 2. CALCULATE SSL MATRICES
# Initialize values for gamma and mu
gamma = 1
gamma = gamma**2
mu = 100
N_nodes = len(uniquePhenotypes)

# Calculate the matrix for SSL given weighted distance matrix
I_ = np.eye(N_nodes)
#W_w = np.exp(-DistW/gamma)
D_w = np.diag(sum(W_w))
L_w = D_w - W_w
SSL_w = np.linalg.inv(I_ + mu*L_w)

# Calculate the matrix for SSL given unweighted distance matrix
#W_u = np.exp(-DistU/gamma)
D_u = np.diag(sum(W_u))
L_u = D_u - W_u
SSL_u = np.linalg.inv(I_ + mu*L_u)



aMatRowNames = list(associationMatrix.index)
D_w_pandas = pd.DataFrame(D_w)
D_w_pandas.columns = aMatRowNames
D_w_pandas.index = aMatRowNames
D_w_pandas.to_csv("/Users/viveksrm/Desktop/Dw.csv")



##############################################
# 3. PROCESS INPUT DATA
# Make lists that correspond to positive labels in our data for 642 and 642.1
class_snp_642 = []
class_snp_6421 = []
class_snp_642and6421 = []
for elem in uniquePhenotypes:
	if(elem == "642"):
		class_snp_642.append(1)
	else:
		class_snp_642.append(0)
	
	if(elem == "642.1"):
		class_snp_6421.append(1)
	else:
		class_snp_6421.append(0)

for (item1, item2) in zip(class_snp_642, class_snp_6421):
    class_snp_642and6421.append(item1+item2)

# Calling DataFrame constructor after zipping 
# both lists, with columns specified 
nodeMap = pd.DataFrame(list(zip(uniquePhenotypes, class_snp_642, class_snp_6421, class_snp_642and6421)), 
               columns =['phenotypeName', '642Label', '6421Label', '642and6421Label']) 

# Convert the input csv file to a numpy array
data = np.array(nodeMap)
# Label y is column 3 - it's a binary class of +1/-1. Unlabeled nodes have a label of 0
y_642 = data[:,1]
y_6421 = data[:,2]
y_642and6421 = data[:,3]


# Get the list of indices where our label is 1
idx_labeled_pos_642 = np.where(y_642 == 1)[0]
# Get the list of indices where our label is -1
idx_labeled_neg_642 = np.where(y_642 == -1)[0]  
# Get the list of indices where our label is known
idx_labeled_642 = np.where(y_642 != 0)[0]

# Get the list of indices where our label is 1
idx_labeled_pos_6421 = np.where(y_6421 == 1)[0]
# Get the list of indices where our label is -1
idx_labeled_neg_6421 = np.where(y_6421 == -1)[0]  
# Get the list of indices where our label is known
idx_labeled_6421 = np.where(y_6421 != 0)[0]

# Get the list of indices where our label is 1
idx_labeled_pos_642and6421 = np.where(y_642and6421 == 1)[0]
# Get the list of indices where our label is -1
idx_labeled_neg_642and6421 = np.where(y_642and6421 == -1)[0]  
# Get the list of indices where our label is known
idx_labeled_642and6421 = np.where(y_642and6421 != 0)[0]



##############################################
# 4. PERFORM LABEL PROPAGATION
def performLabelPropagation(phenotypeId, idsLabeledPos, idsLabeledNeg, idsLabeled, nnodes, sslw, sslu, data):
	print("Performing label propagation on phenotype " + phenotypeId)

	# Apply the Label_Propagation function to come up with label predictions for all nodes
	f_w = Label_Propagation(nnodes, sslw, idsLabeledPos, idsLabeledNeg)
	print("Predictions for all nodes based upon weighted distance matrix:")
	print(f_w)

	f_u = Label_Propagation(nnodes, sslu, idsLabeledPos, idsLabeledNeg)
	print("Predictions for all nodes based upon unweighted distance matrix:")
	print(f_u)

	if(sum(f_w-f_u) == 0):
	    print("Note – predictions are identical between our weighted and unweighted distance matrices")

	print("True labels:")
	data = np.array([int(a) for a in data])
	print(data)
	print()
	return([f_w, f_u])

out_snp_642 = performLabelPropagation("642", idx_labeled_pos_642, idx_labeled_neg_642, idx_labeled_642, N_nodes, SSL_w, SSL_u, y_642)
out_snp_6421 = performLabelPropagation("642.1", idx_labeled_pos_6421, idx_labeled_neg_6421, idx_labeled_642, N_nodes, SSL_w, SSL_u, y_6421)
out_snp_642and6421 = performLabelPropagation("642 and 642.1", idx_labeled_pos_642and6421, idx_labeled_neg_642and6421, idx_labeled_642and6421, N_nodes, SSL_w, SSL_u, y_642and6421)



##############################################
# 5. NORMALIZE LABEL PROPAGATION RESULTS
def normalizeAndProduceResults(outputScores, outputFilePath, yValues, originalDiseases):
	scores = outputScores
	
	final = pd.DataFrame()
	final['phenotype'] = nodeMap['phenotypeName']
	final['score'] = scores
	final['trueLabel'] = yValues

	# Drop the original disease from the dataframe
	for originalDisease in originalDiseases:
		final = final[final.phenotype != originalDisease]

	# Normalize the scores
	sigma = statistics.pstdev(final['score'])
	normalized =  [float(i)/max(final['score']) for i in final['score']]
	transformed = [a/(1+math.exp(-a/sigma)) for a in normalized]
	final['score'] = transformed

	final.to_csv(outputFilePath)

normalizeAndProduceResults(out_snp_642[0], "/Users/viveksrm/Desktop/out_642_w.csv", y_642, ["642"])
normalizeAndProduceResults(out_snp_6421[0], "/Users/viveksrm/Desktop/out_6421_w.csv", y_6421, ["642.1"])
normalizeAndProduceResults(out_snp_642and6421[0], "/Users/viveksrm/Desktop/out_642and6421_w.csv", y_642and6421, ["642","642.1"])


