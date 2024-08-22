import pandas as pd
panglaoDB = pd.read_csv("../resources/PanglaoDB_markers_27_Mar_2020.tsv", sep='\t')
cellMaker = pd.read_excel("../resources/Cell_marker_Seq.xlsx")