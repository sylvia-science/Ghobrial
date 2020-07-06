library('vroom')
file1 = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/AllSamplesDownsample/Batch_Kit/SNN_Umap/connectivities_old.csv'
connect_kit = vroom(file1, delim = ",",col_names = F)

file2 = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/AllSamplesDownsample/Batch_Sample/SNN_Umap/connectivities_old.csv'
connect_sample = vroom(file2, delim = ",",col_names = F)

all((connect_kit[1:10000,1:10000] == connect_sample[1:10000,1:10000]))

all((colSums(connect_kit) == colSums(connect_sample) ))

all((connect_kit == connect_sample))
