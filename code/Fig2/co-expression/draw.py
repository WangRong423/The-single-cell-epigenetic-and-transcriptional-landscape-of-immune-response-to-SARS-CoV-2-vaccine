import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import numpy as np
import time
import os

changecolor = colors.Normalize(vmin=0, vmax=1.0)

annotations = ['Naïve CD4+ T cells', 'Naïve CD8+ T cells', 'NK/NKT', 'Plasmablasts/Memrary B', 'Naïve B',
               'Intermediate B', 'CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'Plasmacytoid DC', 'Treg', 
                'CD8+ Tem', 'MAIT']
# annotations = ['Eryth/Platelet CD14+ Mono new']

for celltype in annotations:
    startTime = time.time()
    print('Start processing {}.'.format(celltype))
    celltype = celltype.replace(' ', '_')
    celltype = celltype.replace('/', '_')

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

    ref = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/filtered/filtered_{}.csv'.format(celltype))
    ref = ref.iloc[:,1:]
    results = ref['clustering results']
    geneNames = ref['gene']

    size = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/getZero/ZeroPT_{}.csv'.format(celltype))
    size = size.fillna(0)
    size_newcol ={}
    for i in range(0,len(size.columns)):
        zero_flip = size.iloc[:,i]
        newcol = size.columns[i] + '_flip'
        size_newcol[newcol] = zero_flip
    size_newcol_df = pd.DataFrame(size_newcol)
    size = pd.concat([size, size_newcol_df], axis=1)

    cv = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/getCV/CV_{}.csv'.format(celltype))
    cv = cv.fillna(0)
    cv_newcol = {}
    for i in range(0,len(cv.columns)):
        cv_flip = cv.iloc[:,i]
        newcol = cv.columns[i] + '_flip'
        cv_newcol[newcol] = cv_flip
    cv_newcol_df = pd.DataFrame(cv_newcol)
    cv = pd.concat([cv, cv_newcol_df], axis=1)

    x = [num for num in range(len(gex1))]
    x_label = [0, 1, 3, 6, 14, 30, 31, 33, 36, 44]

    print('num of genes:{}'.format(len(geneNames)))
    for i in range(len(geneNames)):
        if i % 50 == 0 or i == len(geneNames) - 1:
            print(i)
        geneName = geneNames[i]
        path = '/database/huhuajie/findCoexpressionGenes/finalResults/draw/{}/{}/'.format(celltype, results[i])
        if not os.path.exists(path):
            os.makedirs(path)

        data1 = gex1[geneName]
        data2 = gex2[geneName].iloc[0:10]
        data3 = gex3[geneName]
        data4 = gex4[geneName]
        y_mean, y = [], []
        for j in range(len(data1)):
            y_mean.append((data1[j] + data2[j] + data3[j] + data4[j]) / 4)
            y.append([data1[j], data2[j], data3[j], data4[j]])
        size_tmp = np.array(size[geneName]) * 200
        cv_tmp = np.array(cv[geneName])

        if geneName == 'THRA1/BTR':
            geneName = 'THRA1_BTR'
        if geneName == 'THRA1/BTR_flip':
            geneName = 'THRA1_BTR_flip'

        plt.figure(figsize=(12, 11))

        plt.subplot(2, 1, 1)
        plt.xlim((-0.5,9.5))
        plt.xlabel('Timepoints')
        plt.ylabel('Gene Expression')
        plt.plot(x, y_mean, linewidth=2.5, alpha=0.6,)
        plt.scatter(x, y_mean, marker='o', s=size_tmp, c=cv_tmp, cmap='rainbow', alpha=0.8)
        plt.boxplot(y, sym='o', showmeans=True, meanline=True, meanprops={'linestyle':'-', 'color':'orange'}, medianprops={'linestyle':''}, positions=x)
        plt.colorbar(label="Coefficient of Variation", )
        plt.xticks(x, x_label)
        plt.title('Celltype:{}   Gene:{}'.format(celltype, geneName))
        for s in [50, 100, 200]:
            plt.scatter([], [], c='k', alpha=0.3, s=s,
                        label=str(int(s/2)) + '%')
        plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Non-zero Percentage', loc='best')

        plt.subplot(2, 1, 2)
        plt.xlim((-0.5,9.5))
        plt.xlabel('Timepoints')
        plt.ylabel('Gene Expression')
        plt.plot(x, data1, linewidth=2.5, label = 'P1')
        plt.plot(x, data2, linewidth=2.5, label = 'P2')
        plt.plot(x, data3, linewidth=2.5, label = 'P3')
        plt.plot(x, data4, linewidth=2.5, label = 'P4')
        plt.xticks(x, x_label)
        plt.legend()
        plt.title('Celltype:{}   Gene:{}'.format(celltype, geneName))

        plt.savefig('{}/{} among 4 people.pdf'.format(path, geneName))
        plt.close()

    overTime = time.time()
    spendtime = (overTime - startTime) / 60
    print('{} was Down! Time taken: {} mins'.format(celltype, spendtime))
