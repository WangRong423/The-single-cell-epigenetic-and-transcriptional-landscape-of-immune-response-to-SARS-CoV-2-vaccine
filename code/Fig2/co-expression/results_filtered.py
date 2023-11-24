import operator
import pandas as pd

# annotations = ['Naïve CD4+ T cells', 'Naïve CD8+ T cells', 'NK/NKT', 'Plasmablasts/Memrary B', 'Naïve B',
#                'Intermediate B', 'CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'Plasmacytoid DC', 'Treg', 
#                 'CD8+ Tem', 'Eryth/Platelet', 'MAIT']
annotations = ['Eryth/Platelet CD14+ Mono new']

for annotation in annotations:
    annotation = annotation.replace(' ', '_')
    annotation = annotation.replace('/', '_')
    data = pd.read_csv('/database/huhuajie/findCoexpressionGenes/finalResults/cluster/{}/Results_{}.csv'.format(annotation, annotation))
    alldim = len(data.columns)

    data.insert(loc=alldim-2,column='clustering results',value='')
    data.insert(loc=alldim-1,column='the numbers of genes',value='')
    data = data.rename(columns={'Unnamed: 0':'gene',
                                'The_cohesion_of_the_last_layer':'The cohesion of the last layer',
                                'The_separation_of_the_last_layer':'The separation of the last layer'})

    for i in range(len(data)):
        if data.iloc[i,alldim-7] != 0:
            data.iloc[i,alldim-2] = '{}-{}-{}-{}-{}'.format(int(data.iloc[i,alldim-7]),int(data.iloc[i,alldim-6]),int(data.iloc[i,alldim-5]),int(data.iloc[i,alldim-4]),int(data.iloc[i,alldim-3]))
        else:
            data.iloc[i,alldim-2] = 0
            
    counts = data.groupby(['clustering results']).size()
    for i in range(len(data)):
        data.iloc[i,alldim-1] = counts[data.iloc[i,alldim-2]]
        
    genename = []
    num = []
    for i in range (int(data['layer1 clustering results'].max()) + 1):
        temp_genename = data[data['layer1 clustering results'] == i]['gene'].to_list()
        temp_num = 0
        for j in range(len(temp_genename)):
            if temp_genename[j].find('flip') != -1:
                temp_genename[j] = temp_genename[j].split('_')[0]
                temp_num += 1
        temp_genename.sort()
        num.append(temp_num)
        genename.append(temp_genename)

    delete_list = []
    for i in range(len(genename)):
        for j in range(len(genename)):
            if (i != j) & (operator.eq(genename[i],genename[j])):
                if num[i] > num[j]:
                    delete_list.append(i)
                else:
                    delete_list.append(j)
    delete_list = list(set(delete_list))

    for i in range(len(delete_list)):
        delete_gene = data[data['layer1 clustering results'] == delete_list[i]]
        data = data.drop(delete_gene.index)
        
    data = data.sort_values(by=['the numbers of genes', 'clustering results'], ascending=[False, False])
    data = data.dropna()
    data.to_csv('/database/huhuajie/findCoexpressionGenes/finalResults/filtered/results_{}.csv'.format(annotation))
    newdata = data[data['the numbers of genes'] >= 4]
    newdata = newdata.iloc[:, [0, 32, 33, 34, 35]]
    newdata.to_csv('/database/huhuajie/findCoexpressionGenes/finalResults/filtered/filtered_{}.csv'.format(annotation))
    print('{} was down.'.format(annotation))
