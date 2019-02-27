from knn import NearestNeighborExperiment, NearestNeighborWithPCAExperiment

'''Format: setting, # folds, # principal components'''
hlg_knn = NearestNeighborWithPCAExperiment("HLG", 5, 100)
hlg_knn.run(0, write_log=True)
