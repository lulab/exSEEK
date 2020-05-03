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

def knn_score(X, y, K=10, n_shuffle=None):
    from sklearn.neighbors import NearestNeighbors

    N = X.shape[0]
    assert K < N
    def knn_fractions(X, y):
        nn = NearestNeighbors(K)
        nn.fit(X)
        distances, indices = nn.kneighbors(X, K + 1)
        neighbor_classes = np.take(y, indices[:, 1:])
        return np.sum(neighbor_classes == y[:, np.newaxis], axis=1)
    
    def expected_fractions(X, y):
        stats = np.zeros((n_shuffle, N))
        for i in range(n_shuffle):
            y = np.random.permutation(y)
            stats[i] = knn_fractions(X, y)
        return stats.mean(axis=0)
    
    classes, class_sizes = np.unique(y, return_counts=True)
    classes = np.argmax(y.reshape((-1, 1)) == classes.reshape((1, -1)), axis=1)
    class_sizes = np.take(class_sizes, classes)
    # expected fraction
    mean_r = K/(N - 1)*class_sizes
    observed_r = knn_fractions(X, y)
    #mean_r = expected_fractions(X, y)
    max_r = np.minimum(K, class_sizes)
    #print(observed_r, mean_r, max_r)
    scores = (observed_r - mean_r)/(max_r - mean_r)
    return scores.mean()