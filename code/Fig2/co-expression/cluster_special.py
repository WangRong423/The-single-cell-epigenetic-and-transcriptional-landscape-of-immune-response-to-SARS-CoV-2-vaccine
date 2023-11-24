import os
import gc
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
from sklearn import metrics
from sklearn.cluster import MeanShift
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import Normalizer
from kneed import KneeLocator

warnings.filterwarnings('ignore')
annotations = ['Eryth/Platelet CD14+ Mono new']
# annotations = ['Naïve CD4+ T cells', 'Naïve CD8+ T cells', 'NK/NKT', 'Plasmablasts/Memrary B', 'Naïve B',
#                'Intermediate B', 'CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'Plasmacytoid DC', 'Treg', 
#                 'CD8+ Tem', 'Eryth/Platelet', 'MAIT']
 
for celltype in annotations:
    startTime = time.time()
    print('Start processing {}.'.format(celltype))
    
    celltype = celltype.replace(' ', '_')
    celltype = celltype.replace('/', '_')
    # os.makedirs('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/'.format(celltype))

    slope = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/getSlope/Slope_{}_ED.csv'.format(celltype))
    zero = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/getZero/ZeroPT_{}_ED.csv'.format(celltype))
    zero = zero.fillna(0)
    zero_newcol = {}
    for i in range(len(zero.columns)):
        y = zero.iloc[:, i]
        newcol = zero.columns[i] + '_flip'
        zero_newcol[newcol] = y
    zero_newcol_df = pd.DataFrame(zero_newcol)
    zero = pd.concat([zero, zero_newcol_df], axis=1)
    # cv = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/getCV/all_{}_NED.csv'.format(celltype))
    # cv = cv.fillna(0)
    # cv_newcol = {}
    # for i in range(len(cv.columns)):
    #     y = cv.iloc[:, i]
    #     y = -y
    #     newcol = cv.columns[i] + '_flip'
    #     cv_newcol[newcol] = y
    # cv_newcol_df = pd.DataFrame(cv_newcol)
    # cv = pd.concat([cv, cv_newcol_df], axis=1)
        
    gex1 = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/getGEX/byPeople/P1_{}.csv'.format(celltype))
    gex1 = gex1.fillna(0)
    gex1_newcol = {}
    for i in range(1,len(gex1.columns)):
        y = gex1.iloc[:,i]
        ymax = max(y)
        ymin = min(y)
        yflip = (ymax + ymin) - y
        newcol = gex1.columns[i] + '_flip'
        gex1_newcol[newcol] = yflip
    gex1_newcol_df = pd.DataFrame(gex1_newcol)
    gex1 = pd.concat([gex1, gex1_newcol_df], axis=1)
    
    gex2 = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/getGEX/byPeople/P2_{}.csv'.format(celltype))
    gex2 = gex2.fillna(0)
    gex2_newcol = {}
    for i in range(1,len(gex2.columns)):
        y = gex2.iloc[:,i]
        ymax = max(y)
        ymin = min(y)
        yflip = (ymax + ymin) - y
        newcol = gex2.columns[i] + '_flip'
        gex2_newcol[newcol] = yflip
    gex2_newcol_df = pd.DataFrame(gex2_newcol)
    gex2 = pd.concat([gex2, gex2_newcol_df], axis=1)
    
    gex3 = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/getGEX/byPeople/P3_{}.csv'.format(celltype))
    gex3 = gex3.fillna(0)
    gex3_newcol = {}
    for i in range(1,len(gex3.columns)):
        y = gex3.iloc[:,i]
        ymax = max(y)
        ymin = min(y)
        yflip = (ymax + ymin) - y
        newcol = gex3.columns[i] + '_flip'
        gex3_newcol[newcol] = yflip
    gex3_newcol_df = pd.DataFrame(gex3_newcol)
    gex3 = pd.concat([gex3, gex3_newcol_df], axis=1)
    
    gex4 = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/getGEX/byPeople/P4_{}.csv'.format(celltype))
    gex4 = gex4.fillna(0)
    gex4_newcol = {}
    for i in range(1,len(gex4.columns)):
        y = gex4.iloc[:,i]
        ymax = max(y)
        ymin = min(y)
        yflip = (ymax + ymin) - y
        newcol = gex4.columns[i] + '_flip'
        gex4_newcol[newcol] = yflip
    gex4_newcol_df = pd.DataFrame(gex4_newcol)
    gex4 = pd.concat([gex4, gex4_newcol_df], axis=1)

    geneNames = gex1.columns[1:]
    delete_list = []
    for i in geneNames:
        if sum(gex1[i]) == 0 or sum(gex2[i]) == 0 or sum(gex3[i]) == 0 or sum(gex4[i]) == 0:
            delete_list.append(i)
    gex1 = gex1.drop(delete_list, axis=1)    
    gex2 = gex2.drop(delete_list, axis=1) 
    gex3 = gex3.drop(delete_list, axis=1) 
    gex4 = gex4.drop(delete_list, axis=1) 
    slope = slope.drop(delete_list, axis=1) 
    zero = zero.drop(delete_list, axis=1)   

    merge = pd.concat([slope, zero], axis=0)
    merge = pd.DataFrame(merge.values.T, index=merge.columns)

    num_genes = len(gex1.columns) - 1  
    num_timepoints = len(gex1)
    alldim = len(merge.columns)
    dim_slope_origin = int((len(slope) + 1) / 2)
    dim_slope_gap = len(slope) - dim_slope_origin
    dim_zero = len(zero)
    
    print('='*25 + 'The summary of the data being processed.' + '='*25)
    print('num_genes: ', num_genes)
    print('num_timepoints: ', num_timepoints)
    print('dim_slope_origin: ', dim_slope_origin)
    print('dim_slope_gap: ', dim_slope_gap)
    print('dim_zero: ', dim_zero)
    print('alldim: ', alldim)
    print('='*88)

    del slope,zero
    gc.collect()

    slope_origin = merge.iloc[:, 0:dim_slope_origin]
    slope_origin = slope_origin.copy(deep=True)
    for i in range(len(slope_origin)):
        for j in range(len(slope_origin.columns)):
            if slope_origin.iloc[i,j] > 0:
                slope_origin.iloc[i,j] = 2
            elif slope_origin.iloc[i,j] < 0:
                slope_origin.iloc[i,j] = -2

    slope_origin['sum'] = slope_origin.apply(lambda x: abs(x).sum(), axis=1)
    merge['pre_cluster'] = ''
    pre_cluster_idx = len(merge.columns) - 1
    #print('pre_cluster_idx: ', pre_cluster_idx)
    sum_idx = len(slope_origin.columns) - 1
    #print('sum_idx: ', sum_idx)
    for i in range(len(merge)):
        if slope_origin.iloc[i,sum_idx] == 2 or slope_origin.iloc[i,sum_idx] == 0:
            merge.iloc[i,pre_cluster_idx] = 0
        else:
            merge.iloc[i,pre_cluster_idx] = 1
    # pre_counts = pd.value_counts(merge['pre_cluster'])
    # x = pre_counts.index
    # y = pre_counts
    # plt.figure(figsize=(20,10))
    # plt.title('{}: The numbers of zero and nonzero.'.format(celltype), fontsize=15)
    # plt.bar(x, y)
    # for a,b,i in zip(x,y,range(len(x))):
    #     plt.text(a, b+0.1, "%.0f"%y[len(x)-1-i], ha='center', fontsize=15)
    # plt.savefig('/home/special/user/huhuajie/clusters_final/cluster/{}/RemoveZero.png'.format(celltype))
    # plt.close()

    print('='*25 + 'Start first-layer clustering.' + '='*25)
    os.makedirs('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/'.format(celltype))
    df1 = pd.DataFrame(columns = range(0,dim_slope_origin))
    for i in range(len(merge)):
        if merge.iloc[i, pre_cluster_idx] == 1:
            df1 = df1.append(slope_origin.iloc[i,:sum_idx])

    X = df1
    ms = MeanShift(bandwidth=2, n_jobs = -1).fit(X)
    labels = ms.labels_
    n_clusters = len(set(labels))
    num_layer1 = n_clusters
    print("Estimated number of clusters: %d" % n_clusters)

    df1['cluster'] = labels

    cents = []
    for i in range(n_clusters):
        sum = 0
        cluster_temp = df1.loc[df1['cluster'] == i]
        length = len(cluster_temp.columns) - 1
        for j in range(len(cluster_temp)):
            temp = cluster_temp.iloc[j,0:length]
            sum += temp

        cents.append(sum/len(cluster_temp))

    Separation = []
    Cohesion = []
    for i in range(n_clusters):
        inter_sum = 0
        intra_sum = 0
        for j in range(n_clusters):
            inter_temp = np.sqrt(np.sum((cents[j] - cents[i])**2))
            inter_sum += inter_temp

        if len(cents) == 1:
            Separation.append(0)
        else:
            Separation.append(inter_sum/(len(cents)-1))

        cluster_temp = df1.loc[df1['cluster'] == i]
        length = len(cluster_temp.columns) - 1
        for k in range(len(cluster_temp)):
            intra_temp = np.sqrt(np.sum((cluster_temp.iloc[k,0:length] - cents[i])**2))
            intra_sum += intra_temp

        Cohesion.append(intra_sum/len(cluster_temp))

    merge['layer1 clustering results'] = ''
    for i in range(len(df1)):
        merge.loc[df1.index[i],'layer1 clustering results'] = labels[i]

    plt.figure(figsize=(20,10))
    counts2 = pd.value_counts(labels)
    x = counts2.index
    y = counts2
    plt.title('{}: The number of clusters is {}.'.format(celltype, n_clusters), fontsize=15)
    plt.bar(x, y)
    plt.savefig('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/layer1_NumofClusters.png'.format(celltype))
    plt.close()

    results_clusters = pd.DataFrame(counts2)
    results_clusters.columns=['numbers']
    results_clusters['Cohesion'] = ''
    results_clusters['Separation'] = ''
    for i in range(len(Separation)):
        results_clusters.loc[i,'Cohesion'] = Cohesion[i]
        results_clusters.loc[i,'Separation'] = Separation[i]
    results_clusters.to_csv('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/layer1_results_of_Clustering.csv'.format(celltype))
    print('='*24 + 'First-layer clustering was Down.' + '='*24)


    print('='*25 + 'Start second-layer clustering.' + '='*25)
    merge['layer2 clustering results'] = ''
    num_layer2 = []
    for i in range(num_layer1):
        os.makedirs('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/{}/'.format(celltype, i))
        df2 = merge.loc[merge['layer1 clustering results'] == i]
        df2 = df2.iloc[:,dim_slope_origin:(dim_slope_origin + dim_slope_gap)]
        df2 = df2.copy(deep=True)
        for j in range(len(df2)):
            for k in range(len(df2.columns)):
                if df2.iloc[j,k] > 0:
                    df2.iloc[j,k] = 2
                elif df2.iloc[j,k] < 0:
                    df2.iloc[j,k] = -2

        X = df2
        if len(df2) <= 3:
            labels = [0 for a in range(len(df2))]
        else:
            ms = MeanShift(bandwidth=2, n_jobs=-1).fit(X)
            labels = ms.labels_

        n_clusters = len(set(labels))
        num_layer2.append(n_clusters)
        # print("Estimated number of clusters: %d" % n_clusters)

        df2['cluster'] = labels

        cents = []
        for g in range(n_clusters):
            sum = 0
            cluster_temp = df2.loc[df2['cluster'] == g]
            length = len(cluster_temp.columns) - 1
            for h in range(len(cluster_temp)):
                temp = cluster_temp.iloc[h,0:length]
                sum += temp

            cents.append(sum/len(cluster_temp))

        Separation = []
        Cohesion = []
        for l in range(len(cents)):
            inter_sum = 0
            intra_sum = 0
            for m in range(len(cents)):
                inter_temp = np.sqrt(np.sum((cents[m] - cents[l])**2))
                inter_sum += inter_temp

            if len(cents) == 1:
                Separation.append(0)
            else:
                Separation.append(inter_sum/(len(cents)-1))

            cluster_temp = df2.loc[df2['cluster'] == l]
            length = len(cluster_temp.columns) - 1
            for n in range(len(cluster_temp)):
                intra_temp = np.sqrt(np.sum((cluster_temp.iloc[n,0:length] - cents[l])**2))
                intra_sum += intra_temp

            Cohesion.append(intra_sum/len(cluster_temp))

        for j in range(len(df2.index)):
            merge.loc[df2.index[j],'layer2 clustering results'] = labels[j]

        plt.figure(figsize=(20,10))
        counts2 = pd.value_counts(labels)
        x = counts2.index
        y = counts2
        plt.title('{}: The number of clusters is {}.'.format(celltype, n_clusters), fontsize=15)
        plt.bar(x, y)
        plt.savefig('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/{}/layer2_NumofClusters.png'.format(celltype, i))
        plt.close()

        results_clusters = pd.DataFrame(counts2)
        results_clusters.columns=['numbers']
        results_clusters['Cohesion'] = ''
        results_clusters['Separation'] = ''
        for k in range(len(Separation)):
            results_clusters.loc[k,'Cohesion'] = Cohesion[k]
            results_clusters.loc[k,'Separation'] = Separation[k]
        results_clusters.to_csv('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/{}/layer2_Results_of_Clustering.csv'.format(celltype, i))

        # print('Layer2: No.' + str(i) + ' was down!')
    print('='*24 + 'Second-layer clustering was down.' + '='*24)

    print('='*25 + 'Start third-layer clustering.' + '='*25)
    merge['layer3 clustering results'] = ''
    num_layer3 = []
    for i in range(num_layer1):
        num_layer3_temp = []
        for j in range(num_layer2[i]):
            df3 = merge.loc[(merge['layer1 clustering results'] == i) & (merge['layer2 clustering results'] == j)]
            os.makedirs('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/{}/{}/'.format(celltype, i, j))
            df3 = df3.iloc[:,dim_slope_gap:(dim_slope_origin + dim_slope_gap)]
            df3 = df3.copy(deep=True)

            X = Normalizer().fit_transform(df3)

            if len(df3) > 3:
                nn = NearestNeighbors(n_neighbors=2, n_jobs=-1).fit(X)
                distances, idx = nn.kneighbors(X)
                distances = np.sort(distances, axis=0)
                distances = distances[:,1]
                x = np.arange(len(distances))
                y = distances
                kneedle = KneeLocator(x, y, S=1.0, curve="convex", direction="increasing")
                # print('The elbow is: ', kneedle.elbow_y)

            if (len(df3) <= 3) or (kneedle.elbow_y == 0) or (kneedle.elbow_y == 'None') or (kneedle.elbow_y == None):
                labels = [0 for k in range(len(df3))]
            else:
                bandwidth = kneedle.elbow_y
                ms = MeanShift(bandwidth=bandwidth, n_jobs=-1).fit(X)
                labels = ms.labels_


            n_clusters = len(set(labels))
            num_layer3_temp.append(n_clusters)
            # print("Estimated number of clusters: %d" % n_clusters)

            df3['cluster'] = labels

            cents = []
            for g in range(n_clusters):
                sum = 0
                cluster_temp = df3.loc[df3['cluster'] == g]
                length = len(cluster_temp.columns) - 1
                for h in range(len(cluster_temp)):
                    temp = cluster_temp.iloc[h,0:length]
                    sum += temp

                cents.append(sum/len(cluster_temp))

            Separation = []
            Cohesion = []
            for l in range(len(cents)):
                inter_sum = 0
                intra_sum = 0
                for m in range(len(cents)):
                    inter_temp = np.sqrt(np.sum((cents[m] - cents[l])**2))
                    inter_sum += inter_temp

                if len(cents) == 1:
                    Separation.append(0)
                else:
                    Separation.append(inter_sum/(len(cents)-1))

                cluster_temp = df3.loc[df3['cluster'] == l]
                length = len(cluster_temp.columns) - 1
                for n in range(len(cluster_temp)):
                    intra_temp = np.sqrt(np.sum((cluster_temp.iloc[n,0:length] - cents[l])**2))
                    intra_sum += intra_temp

                Cohesion.append(intra_sum/len(cluster_temp))

            for a in range(len(df3.index)):
                merge.loc[df3.index[a],'layer3 clustering results'] = labels[a]


            plt.figure(figsize=(20,10))
            counts3 = pd.value_counts(labels)
            x = counts3.index
            y = counts3
            plt.title('{}: The number of clusters is {}.'.format(celltype, n_clusters), fontsize=15)
            plt.bar(x, y)
            plt.savefig('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/{}/{}/layer3_NumofClusters.png'.format(celltype, i, j))
            plt.close()

            results_clusters = pd.DataFrame(counts3)
            results_clusters.columns=['numbers']
            results_clusters['Cohesion'] = ''
            results_clusters['Separation'] = ''
            for k in range(len(Separation)):
                results_clusters.loc[k,'Cohesion'] = Cohesion[k]
                results_clusters.loc[k,'Separation'] = Separation[k]
            results_clusters.to_csv('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/{}/{}/layer3_Results_of_Clustering.csv'.format(celltype, i, j))
            # print('Layer3: No.' + str(j) + ' in cluster ' + str(i) + '(Layer2) was down!')

        num_layer3.append(num_layer3_temp)
    print('='*25 + 'Third-layer clustering was down.' + '='*25)

    merge['layer4 clustering results'] = ''
    merge['The_cohesion_of_the_last_layer'] = ''
    merge['The_separation_of_the_last_layer'] = ''
    num_layer4 = []
    for i in range(num_layer1):
        for j in range(num_layer2[i]):
            for k in range(num_layer3[i][j]):
                df4 = merge.loc[(merge['layer1 clustering results'] == i) & (merge['layer2 clustering results'] == j) & (merge['layer3 clustering results'] == k)]
                os.makedirs('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/{}/{}/{}/'.format(celltype, i, j, k))
                df4 = df4.iloc[:,(dim_slope_origin + dim_slope_gap):(dim_slope_origin + dim_slope_gap + dim_zero)]
                df4 = df4.copy(deep=True)

                X = Normalizer().fit_transform(df4)

                if len(df4) > 3:
                    nn = NearestNeighbors(n_neighbors=2, n_jobs=-1).fit(X)
                    distances, idx = nn.kneighbors(X)
                    distances = np.sort(distances, axis=0)
                    distances = distances[:,1]
                    x = np.arange(len(distances))
                    y = distances
                    kneedle = KneeLocator(x, y, S=1.0, curve="convex", direction="increasing")
                    # print('The elbow is: ', kneedle.elbow_y)

                if (len(df4) <= 3) or (kneedle.elbow_y == 0) or (kneedle.elbow_y == 'None') or (kneedle.elbow_y == None):
                    labels = [0 for b in range(len(df4))]
                else:
                    bandwidth = kneedle.elbow_y
                    ms = MeanShift(bandwidth=bandwidth, n_jobs = -1).fit(X)
                    labels = ms.labels_

                n_clusters = len(set(labels))
                num_layer3_temp.append(n_clusters)
                # print("Estimated number of clusters: %d" % n_clusters)

                df4['cluster'] = labels

                cents = []
                for g in range(n_clusters):
                    sum = 0
                    cluster_temp = df4.loc[df4['cluster'] == g]
                    length = len(cluster_temp.columns) - 1
                    for h in range(len(cluster_temp)):
                        temp = cluster_temp.iloc[h,0:length]
                        sum += temp

                    cents.append(sum/len(cluster_temp))

                Separation = []
                Cohesion = []
                for l in range(len(cents)):
                    inter_sum = 0
                    intra_sum = 0
                    for m in range(len(cents)):
                        inter_temp = np.sqrt(np.sum((cents[m] - cents[l])**2))
                        inter_sum += inter_temp

                    if len(cents) == 1:
                        Separation.append(0)
                    else:
                        Separation.append(inter_sum/(len(cents)-1))

                    cluster_temp = df4.loc[df4['cluster'] == l]
                    length = len(cluster_temp.columns) - 1
                    for n in range(len(cluster_temp)):
                        intra_temp = np.sqrt(np.sum((cluster_temp.iloc[n,0:length] - cents[l])**2))
                        intra_sum += intra_temp

                    Cohesion.append(intra_sum/len(cluster_temp))

                for a in range(len(df4.index)):
                    merge.loc[df4.index[a],'layer4 clustering results'] = labels[a]
                    merge.loc[df4.index[a],'The_cohesion_of_the_last_layer'] = Cohesion[labels[a]]
                    merge.loc[df4.index[a],'The_separation_of_the_last_layer'] = Separation[labels[a]]

                plt.figure(figsize=(20,10))
                counts4 = pd.value_counts(labels)
                x = counts4.index
                y = counts4
                plt.title('{}: The number of clusters is {}.'.format(celltype, n_clusters), fontsize=15)
                plt.bar(x, y)
                plt.savefig('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/{}/{}/{}/layer4_NumofClusters.png'.format(celltype, i, j, k))
                plt.close()

                results_clusters = pd.DataFrame(counts4)
                results_clusters.columns=['numbers']
                results_clusters['Cohesion'] = ''
                results_clusters['Separation'] = ''
                for b in range(len(Separation)):
                    results_clusters.loc[b,'Cohesion'] = Cohesion[b]
                    results_clusters.loc[b,'Separation'] = Separation[b]
                results_clusters.to_csv('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/1/{}/{}/{}/layer4_Results_of_Clustering.csv'.format(celltype, i, j, k))

                # print('Layer4: No.' + str(k) + ' in cluster ' + str(j) + '(Layer3) was down!')


    x = np.arange(num_timepoints)
    for i in range(1, num_genes + 1):
        geneName = merge.index[i-1]

        y1 = gex1[geneName]
        y2 = gex2[geneName]
        y3 = gex3[geneName]
        y4 = gex4[geneName]
        y_upper, y_lower, y = [], [], []
        plt.figure(figsize = (6, 3))
        plt.xlim((-0.5,9.5))
        plt.xticks(x)
        for k in range(len(y1)):
          y_tmp = (y1[k] + y2[k] + y3[k] + y4[k]) / 4
          ymax = (max(y1[k], y2[k], y3[k], y4[k]))
          ymin = (min(y1[k], y2[k], y3[k], y4[k]))
          y_upper.append(ymax - y_tmp)
          y_lower.append(y_tmp - ymin)
          y.append(y_tmp)
        plt.errorbar(x, y, yerr = [y_lower, y_upper], fmt='o-', ecolor='r', color='b', elinewidth=1, capsize=4)
        path = '/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/'.format(celltype)
        for j in [alldim, alldim + 1, alldim + 2, alldim + 3, alldim + 4]:
            if merge.iloc[i-1,j] != None:
                path = path + str(merge.iloc[i-1,j]) + '/'
        if not os.path.exists(path):
            os.makedirs(path)
        if geneName == 'THRA1/BTR':
          geneName = 'THRA1_BTR'
        if geneName == 'THRA1/BTR_flip':
          geneName = 'THRA1_BTR_flip'
        plt.savefig('{}{}.png'.format(path, geneName))
        plt.close()
        if i % 1000 == 0 or i == num_genes:
            print(str(i) + '/' + str(num_genes) + ' was down.')

    merge.to_csv('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/Results_{}.csv'.format(celltype,celltype))

    overTime = time.time()
    spendtime = (overTime - startTime) / 60
    print('{} was Down! Time taken: {} mins'.format(celltype, spendtime))
