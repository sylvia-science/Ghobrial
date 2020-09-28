import signatureanalyzer as sa
import pandas as pd

# ---------------------
# RUN SIGNATURE ANALYZER
# ---------------------

celltype = 'Stem_Cell'
celltype = ''
path = "/home/sujwary/Desktop/scRNA/Data/NMF/Harmony_AllSamples_Sample_Kit_" + celltype
input_matrix = pd.read_csv(path + ".tsv", sep='\t')
sa.run_matrix(matrix=input_matrix, outdir= path + '_phi1_alpha10/', nruns=100,verbose=True,plot_results=False,K0=20,tolerance=1e-6,objective='gaussian', max_iter=30000, phi=1.0, a=10.0)
 


text_file = open("/home/sujwary/Desktop/scRNA/Data/NMF/PrePostEOTNBM_MT15_All_phi1_alpha10/Output.txt", "w")
text_file.write("Finished")
text_file.close()



# ---------------------
# LOADING RESULTS
# ---------------------
import pandas as pd

H = pd.read_hdf('nmf_output.h5', 'H')
W = pd.read_hdf('nmf_output.h5', 'W')
Hraw = pd.read_hdf('nmf_output.h5', 'Hraw')
Wraw = pd.read_hdf('nmf_output.h5', 'Wraw')
feature_signatures = pd.read_hdf('nmf_output.h5', 'signatures')
markers = pd.read_hdf('nmf_output.h5', 'markers')
cosine = pd.read_hdf('nmf_output.h5', 'cosine')
log = pd.read_hdf('nmf_output.h5', 'log')

# Output for each run may be found at...
Hrun1 = pd.read_hdf('nmf_output.h5', 'run1/H')
Wrun1 = pd.read_hdf('nmf_output.h5', 'run1/W')
# etc...

# Aggregate output information for each run
aggr = pd.read_hdf('nmf_output.h5', 'aggr')

# ---------------------
# PLOTTING
# ---------------------
sa.pl.marker_heatmap(...)
sa.pl.signature_barplot(...)
sa.pl.stacked_bar(...)
sa.pl.k_dist(...)
sa.pl.consensus_matrix(...)