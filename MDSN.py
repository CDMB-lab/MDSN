import math

import matlab
import matlab.engine
import pandas as pd
import numpy as np
import scipy.io as sio
from datetime import datetime
import os
from itertools import combinations
import commontools as comm

engine = matlab.engine.start_matlab()
# dataPath = 'output/high_order'
dataFile = 'threshold_7_0.5_h0.40_EDM-1_01'
maxOrder = 3
t_threshold = 2
CLUSTERING_NUM = 5
# [0.4, 0.8]
NEIGHBOR_RANGE = [0, 1]
SPAN = 3

def threshold(scores, idx):
    # Score = multiSURFScores(:, 2)';
    # num = length(Score);
    # SlidingWindow = max(min(10, num / 10), 5);
    # Score2 = zeros(1, num - 2);
    score = np.transpose(scores[:, idx])
    num = score.shape[0]
    sliding_window = max(min(10, num / 10), 5)
    score2 = np.zeros(num - 2)
    for i in range(num - 2):
        score2[i] = score[i + 2] - 2 * score[i + 1] + score[i]

    mean = np.average(score2)
    std = np.std(score2)
    outliers = np.where(np.logical_and(score2 <= (mean + std / 2), score2 >= (mean - std / 2)) == False)[0]
    num_outliers = outliers.shape[0]
    topN = num_outliers
    for i in range(num_outliers - 1):
        if (outliers[i + 1] - outliers[i]) > sliding_window:
            topN = outliers[i] + 1
            break

    topN = max(min([topN, 100]), 25)
    return topN


def association(dataFile, mtsf_scores, topN, maxOrder):
    # global engine
    snps = np.transpose(mtsf_scores[0:topN + 1, 0])
    tmp = sio.loadmat(dataFile + '.mat')
    pts = tmp['pts']
    clas = tmp['class']
    values = np.empty((maxOrder - 1, 1), dtype=object)
    for i in range(2, maxOrder + 1):
        values[i - 2, 0] = np.array(list(combinations(snps, i)), dtype=float)
        print('{} {}'.format(i, values[i - 2, 0].shape[0]))

    for i in range(maxOrder - 1):
        [num, order] = values[i, 0].shape
        # values[i, 0] = np.insert(values[i, 0], order, values=np.zeros(num), axis=1)
        # for j in range(num):
        # values[i, 0][j, order] = engine.MutualInformation(
        #     matlab.double(pts.tolist()),
        #     matlab.double(clas.tolist()),
        #     matlab.double(list(values[i, 0][j, 0:2])))
        mutualInfo = engine.mutualInfoTool(num, order,
                                           matlab.double(pts.tolist()),
                                           matlab.double(clas.tolist()),
                                           matlab.double(values[i, 0].tolist()))
        values[i, 0] = np.insert(values[i, 0], order, values=np.array(mutualInfo).transpose(), axis=1)
        a = np.sort(values[i, 0][:, order])[::-1]
        b = np.argsort(values[i, 0][:, order])[::-1]
        values[i, 0][:, 0:order] = values[i, 0][b, 0:order]
        values[i, 0][:, order] = a

    outFile = dataFile + '_Association.mat'
    sio.savemat(outFile, {'values': values, 'topN': topN})
    return values


def transCytoscape(dataFile, Edges, factor):
    # edges_txt = Edges
    Edges = np.insert(Edges, 2, values=np.zeros((3, Edges.shape[0])), axis=1).astype(int)
    for i in range(2):
        for j in range(factor.shape[1]):
            Edges[np.where(Edges[:, i] == factor[0, j]), i + 3] = 1

    Edges[:, 2] = 1
    # Edges = Edges.astype(int)
    # write down Edges
    with open(dataFile + '_cytoscape.txt', 'w') as f:
        f.write('source\t')
        f.write('target\t')
        f.write('weights\t')
        f.write('SouAttri\t')
        f.write('TarAttri\n')

        for i in range(Edges.shape[0]):
            for j in range(Edges.shape[1] - 1):
                f.write('{}\t'.format(Edges[i, j]))
            f.write('{}\n'.format(Edges[i, Edges.shape[1] - 1]))


def transAdjacentMatrix(Edges):
    vertexTable = np.unique(Edges)
    numIndex = vertexTable.shape[0]
    vertexTable = np.insert(vertexTable.reshape(numIndex, 1), 1, values=np.array(range(numIndex)), axis=1).astype(int)
    adjacentMatrix = np.zeros((numIndex, numIndex), dtype=int)
    for i in range(numIndex):
        loc = np.where(Edges == vertexTable[i, 0])
        # loc = np.insert(np.transpose(loc_tpl[0]).reshape(loc_tpl[0].shape[0], 1), 1, values=loc_tpl[1], axis=1)
        Edges[loc[0], loc[1]] = vertexTable[i, 1]

    for i in range(Edges.shape[0]):
        adjacentMatrix[Edges[i, 0], Edges[i, 1]] = 1
        adjacentMatrix[Edges[i, 1], Edges[i, 0]] = 1

    return adjacentMatrix, vertexTable


def pruneMatrix(adjacentMatrix, vertexTable):
    vertex = list(vertexTable[:, 1])
    numVertex = adjacentMatrix.shape[0]
    # vertex = list(range(numVertex))
    numEdge = np.sum(adjacentMatrix) / 2
    maxDensity = 2 * numEdge / (numVertex * (numVertex - 1))
    if numVertex < maxOrder:
        mtx_pruned = adjacentMatrix
        tbl_pruned = vertexTable[:, 0]
        return mtx_pruned, tbl_pruned
    while len(vertex) >= maxOrder:
        current_adj = adjacentMatrix[vertex][:, vertex]
        degrees = np.sum(current_adj, axis=0)
        k = min(degrees)
        current_numVertex = len(vertex)
        current_numEdge = np.sum(current_adj) / 2
        density = 2 * current_numEdge / (current_numVertex * (current_numVertex - 1))
        if density >= maxDensity:
            mtx_pruned = current_adj
            tbl_pruned = vertexTable[vertex, 0]
        delVertex = np.where(degrees == k)[0]
        vertex = [vertex[x] for x in range(current_numVertex) if x not in delVertex]
    return mtx_pruned, tbl_pruned


def getDegreeAndDensity(matrix, vtxTable):
    mtx, tbl = pruneMatrix(matrix, vtxTable)
    degrees = np.sum(mtx, axis=0)
    degree_avg = np.average(degrees)
    numVertex = mtx.shape[0]
    numEdge = np.sum(mtx) / 2
    density = 2 * numEdge / (numVertex * (numVertex - 1))
    return degree_avg, density


def vertexWeighting(adjacentMatrix):
    numVertex = adjacentMatrix.shape[0]
    v_scores = np.zeros((1, numVertex))
    v_weights = np.ones((1, numVertex))
    e_weights = np.zeros((numVertex, numVertex))
    for i in range(numVertex):
        neighbors = np.where(adjacentMatrix[i, :] == 1)[0]
        neighbors_idx = np.insert(neighbors, 0, values=np.array([i]))
        neighbors_idx_table = np.insert(neighbors_idx.reshape(neighbors_idx.shape[0], 1), 1,
                                        values=np.array(range(neighbors_idx.shape[0])), axis=1)
        neighbors_adj = adjacentMatrix[neighbors_idx][:, neighbors_idx]
        # mtx_pruned = pruneMatrix(neighbors_adj)
        degree_avg, density_max = getDegreeAndDensity(neighbors_adj,
                                                      neighbors_idx_table)  # 与pruneMatrix一同改写成两参数，adj_mtx idx_table
        v_scores[0, i] = degree_avg * density_max

    t = t_threshold
    while t > 0:
        for i in range(numVertex):
            neighbors = np.where(adjacentMatrix[i, :] == 1)[0]
            for j in neighbors:
                u = np.where((adjacentMatrix[j, :] == 1) & (adjacentMatrix[i, :] == 1))[0]
                e_weights[i, j] = v_weights[0, i] * v_scores[0, i] + v_weights[0, j] * v_scores[0, j] + np.sum(
                    v_weights[0, u] * v_scores[0, u])

        for i in range(numVertex):
            neighbors = np.where(adjacentMatrix[i, :] == 1)[0]
            v_weights[0, i] = np.sum(e_weights[i, neighbors])

        t = t - 1
    return v_weights


def cutCommunities(communities, snpTable, communitiesTag, seed, span):
    global NEIGHBOR_RANGE
    delta = abs(NEIGHBOR_RANGE[1] - NEIGHBOR_RANGE[0])
    down = NEIGHBOR_RANGE[0] + (communitiesTag-1)*0.05*delta
    up = NEIGHBOR_RANGE[1] - (communitiesTag-1)*0.05*delta
    span = span - 1
    if span <= 0:
        return snpTable
    neighbors = np.where((communities[seed, :] == 1) & (snpTable[:, 2] == 0))[0]
    # snpTable[neighbors[0:math.ceil(len(neighbors)*0.05)], 2] = communitiesTag
    snpTable[neighbors[math.ceil(len(neighbors)*down):math.ceil(len(neighbors)*up)], 2] = communitiesTag
    for i in range(neighbors.shape[0]):
        snpTable = cutCommunities(communities, snpTable, communitiesTag, neighbors[i], span)
    return snpTable


def getEpistasis(cluster, clusterTable):
    epistaticNetwork, epistasis = pruneMatrix(cluster, clusterTable)
    return epistaticNetwork, epistasis


def detectingCommunities(weights, vertexTable, adjacentMatrix, factor):
    global CLUSTERING_NUM
    # a = np.sort(weights)[0, ::-1]
    index = np.argsort(weights)[0, ::-1]
    snps = vertexTable[index, 0]
    communities = adjacentMatrix[index][:, index]
    snpTable = np.zeros((snps.shape[0], 3), dtype=int)
    snpTable[:, 0] = snps
    snpTable[:, 1] = np.array(range(snps.shape[0]))
    meta_snpTable = snpTable * 1

    table_list = []
    community_num = CLUSTERING_NUM
    communitiesTag = 1
    while communitiesTag <= community_num:
        seed = np.where(snpTable[:, 2] == 0)[0][communitiesTag-1]
        # communitiesTag = communitiesTag + 1
        snpTable[seed, 2] = communitiesTag
        # span = 3
        snpTable = cutCommunities(communities, snpTable, communitiesTag, seed, SPAN)
        communitiesTag = communitiesTag + 1
        # index_sorted = np.argsort(snpTable[:, 2], kind='mergesort')
        # snpTable = snpTable[index_sorted, :]
        table_list.append(snpTable)
        snpTable = meta_snpTable * 1

    # numCom = np.max(snpTable[:, 2])
    numCom = community_num
    epistasis = np.empty((numCom, 1), dtype=object)
    epistaticNetwork = np.empty((numCom, 1), dtype=object)
    locations = np.empty((numCom, 1), dtype=object)
    count = np.zeros((1, numCom), dtype=int)
    for i in range(numCom):
        tag = np.where(table_list[i][:, 2] == i + 1)[0]
        cluster = communities[tag][:, tag]
        clusterTable = np.insert(table_list[i][tag, 0].reshape(tag.shape[0], 1), 1, values=np.array(range(tag.shape[0])),
                                 axis=1)
        epistaticNetwork[i, 0], epistasis[i, 0] = getEpistasis(cluster, clusterTable)
        locations[i, 0] = comm.ismember(factor[0], epistasis[i, 0])
        count[0, i] = np.where(locations[i, 0] > -1)[0].size

    print('====Communities====')
    print('factor: ', factor)
    for i in range(numCom):
        print('Commnunity:')
        print(epistasis[i, 0])
        print('Locations:')
        print(locations[i, 0])
        print(count[0, i])
    if numCom == 1:
        epistasis = epistasis[0, 0]
        epistaticNetwork = epistaticNetwork[0, 0]
        locations = locations[0, 0]
        count = count[0, 0]


    return epistasis, epistaticNetwork, locations, count
    # return epistasis, epistaticNetwork


def run(filename, t_max=2, max_order=3):
    # global engine
    # dataFile = 'threshold_7_0.5_h0.40_EDM-1_10'
    # maxOrder = 3
    global dataFile
    global t_threshold, maxOrder
    t_threshold = t_max
    # dataPath = data_path
    dataFile = filename
    maxOrder = max_order

    tic = datetime.now()
    print('==========multiSURF==========')
    if os.path.exists(dataFile + '_multiSURF.mat') is False:
        engine.multiSURF(dataFile)

    mlts_mat = sio.loadmat(dataFile + '_multiSURF.mat')
    mtsf_scores = mlts_mat['multiSURFScores']
    mtsf_locations = mlts_mat['locations_multiSURF']
    factor = mlts_mat['factor']
    toc = datetime.now()
    print('Elapsed time: {} s'.format((toc - tic).total_seconds()))

    print('==========Threshold==========')
    topN = threshold(mtsf_scores, 1)
    print('threshold: {}'.format(topN))
    print('locations: ')
    print(mtsf_locations)
    toc2 = datetime.now()
    print('Elapsed time: {} s'.format((toc2 - toc).total_seconds()))

    print('==========Association==========')
    # values = []
    # tmp = sio.loadmat(dataFile + '.mat')
    # factor = tmp['factor']
    # topN = tmp['pts'].shape[1]
    # mtsf_scores = np.array(range(1, topN+1), dtype=float).reshape(topN, 1)

    if os.path.exists(dataFile + '_Association.mat') is False:
        values = association(dataFile, mtsf_scores, topN, maxOrder)
    else:
        association_mat = sio.loadmat(dataFile + '_Association.mat')
        values = association_mat['values']

    toc3 = datetime.now()
    print('Elapsed time: {} s'.format((toc3 - toc2).total_seconds()))

    print('==========EdgeMatrix==========')
    numOrder = values.shape[0]
    topN_edges = np.zeros((numOrder, 1), dtype=int)
    for i in range(numOrder):
        topN_edges[i, 0] = threshold(values[i, 0], i + 2)
    Edges = values[0, 0][0:topN_edges[0, 0], 0:-1]
    for i in range(1, numOrder):
        size = Edges.shape[0]
        for j in range(topN_edges[i, 0]):
            Edges = np.insert(Edges, size, values=np.array(list(combinations(values[i, 0][j, 0:-1], 2))), axis=0)

    Edges = np.unique(np.sort(Edges, axis=1), axis=0).astype(int)
    # Edges = Edges.astype(int)
    TopCombinations = np.empty((numOrder, 1), dtype=object)
    for i in range(numOrder):
        TopCombinations[i, 0] = values[i, 0][0:topN_edges[i, 0], 0:-1]
    toc4 = datetime.now()
    print('Elapsed time: {} s'.format((toc4 - toc3).total_seconds()))

    print('==========TransCytoscape==========')
    transCytoscape(dataFile, Edges, factor)
    toc5 = datetime.now()
    print('Elapsed time: {} s'.format((toc5 - toc4).total_seconds()))

    print('==========TransAdjacentMatrix==========')
    adjacentMatrix, vertexTable = transAdjacentMatrix(Edges)
    toc6 = datetime.now()
    print('Elapsed time: {} s'.format((toc6 - toc5).total_seconds()))

    print('==========HSIDE==========')
    # VertexWeighting
    print('====VertexWeighting====')
    weights = vertexWeighting(adjacentMatrix)
    toc7 = datetime.now()
    print('Elapsed time: {} s'.format((toc7 - toc6).total_seconds()))
    print('====DetectingCommunities====')
    # weights = vertexWeighting(adjacentMatrix)
    epistasis, epistaticNetwork, locations, count = detectingCommunities(weights, vertexTable, adjacentMatrix, factor)

    outFile = dataFile + '_Epistasis.mat'
    toc8 = datetime.now()
    totalTime = toc8 - tic
    sio.savemat(outFile, {'epistasis': epistasis, 'epistaticNetwork': epistaticNetwork, 'factor': factor,
                          'locations': locations, 'count': count, 'time': totalTime.total_seconds()})
    print('Elapsed time: {} s'.format((toc8 - toc7).total_seconds()))
    print('Total time: {}s'.format((toc8 - tic).total_seconds()))


if __name__ == '__main__':
    # engine = matlab.engine.start_matlab()
    # engine.multiSURF('threshold_7_0.5_h0.40_EDM-1_10.mat')
    run('additive_6_0.5_h0.05_EDM-1_01')

    # run('output/high_order_attr1000/additive_7_0.5_h0.40_EDM-1/additive_7_0.5_h0.40_EDM-1_06')
