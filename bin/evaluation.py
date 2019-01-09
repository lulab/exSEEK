import numpy as np
from sklearn.utils.linear_assignment_ import linear_assignment
from sklearn.mixture import GaussianMixture as GMM
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

def convert_label_to_int(sample_class):
    classes, counts = np.unique(sample_class, return_counts=True)
    classes = np.argmax(sample_class.reshape((-1, 1)) == classes.reshape((1, -1)), axis=1)
    return classes

def unsupervised_clustering_accuracy(y, y_pred):
    assert len(y_pred) == len(y)
    u = np.unique(np.concatenate((y, y_pred)))
    n_clusters = len(u)
    mapping = dict(zip(u, range(n_clusters)))
    reward_matrix = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    for y_pred_, y_ in zip(y_pred, y):
        if y_ in mapping:
            reward_matrix[mapping[y_pred_], mapping[y_]] += 1
    cost_matrix = reward_matrix.max() - reward_matrix
    ind = linear_assignment(cost_matrix)
    return sum([reward_matrix[i, j] for i, j in ind]) * 1.0 / y_pred.size, ind

def uca_score(X, y, prediction_algorithm='knn'):
    #X_log = np.log2(X + 0.001).T
    X_scale = StandardScaler().fit_transform(X)
    # convert string labels to integer labels
    unique_classes = np.unique(y)
    labels = np.zeros(y.shape)
    for i, c in enumerate(unique_classes):
        labels[y == c] = i
    
    cluster_num = np.unique(y).shape[0]
    if prediction_algorithm == 'knn':
        labels_pred = KMeans(cluster_num, n_init=200).fit_predict(X_scale)  
    elif prediction_algorithm == 'gmm':
        gmm = GMM(cluster_num)
        gmm.fit(X_scale)
        labels_pred = gmm.predict(X_scale)
    labels_int = convert_label_to_int(labels)
    score = unsupervised_clustering_accuracy(labels_int, labels_pred)[0]
    return score