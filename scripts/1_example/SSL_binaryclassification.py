# Import statements
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn import metrics
import random

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
# 2. PREPARE INPUT DATA
# Read in the input sample data
data = pd.read_csv("sample_data.csv", sep=",")
# Convert the input csv file to a numpy array
data = np.array(data)

# X includes columns 2 onward
X = data[:,2:]
# Label y is column 1 - it's a binary class of +1/-1. Unlabeled nodes have a label of 0
y = data[:,1]

# Get the list of indices where our label is 1
idx_labeled_pos = np.where(y == 1)
idx_labeled_pos = idx_labeled_pos[0]
# Get the list of indices where our label is -1
idx_labeled_neg = np.where(y == -1)  
idx_labeled_neg = idx_labeled_neg[0]
# Get the list of indices where our label is known
idx_labeled = np.where(y != 0)
idx_labeled = idx_labeled[0]

# Scaler transforms features by scaling each feature to a given range
Scaler = MinMaxScaler()
# Get the min and max of X, and scale all elements to the range (0, 1)
scaled_data = Scaler.fit(X)
scaled_data = Scaler.transform(X)
X = scaled_data


##############################################
# 3. CALCULATE SSL MATRIX
# Get the number of nodes
N_nodes = len(X)

# Compute the distance between all pairs of rows of our data
Dist = metrics.pairwise.euclidean_distances(X) 

# Initialize values for gamma and mu
gamma = 1
gamma = gamma^2
mu = 0.5

# Calculate the matrix for SSL
I_ = np.eye(N_nodes)
W = np.exp(-Dist/gamma)
D = np.diag(sum(W))
L = D - W

SSL = np.linalg.inv(I_ + mu*L)
#print(SSL.shape)


##############################################
# 4. PERFORM LABEL PROPAGATION USING A RANDOM SUBSET OF LABELED DATA
# In our worst case, we have just 1 node labeled positively and 1 labeled negatively
N_pos = 5
N_neg = 5

# Randomly sample N positive and negative nodes from our positively and negatively labeled nodes
pos_test = random.sample(range(len(idx_labeled_pos)-1), N_pos)
neg_test = random.sample(range(len(idx_labeled_neg)-1), N_neg)
set_positive = idx_labeled_pos[pos_test]
set_negative = idx_labeled_neg[neg_test]

# Make a union of our two sets to figure out a list of labeled node indices
labeled = np.union1d(set_positive, set_negative)
# Get the indices for all labeled data that weren't in our random sample
test_set = np.setdiff1d(idx_labeled, labeled)

# Apply the Label_Propagation function to come up with label predictions for all nodes
f_ = Label_Propagation(N_nodes, SSL, set_positive, set_negative)

print(f_)

# Evaluate accuracies from our model
fpr, tpr, thresholds = metrics.roc_curve(y[test_set], f_[test_set], pos_label=1)
print(metrics.auc(fpr, tpr))
