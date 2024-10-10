import pandas as pd
import pkg_resources
panglao_path = pkg_resources.resource_filename('pairpot', 'resources/PanglaoDB_markers_27_Mar_2020.tsv')
print("PanglaoDB path:"+panglao_path)
panglaoDB = pd.read_csv(panglao_path, sep='\t')

marker_path=pkg_resources.resource_filename('pairpot','resources/Cell_marker_Seq.xlsx')
print("CellMarker path:"+marker_path)
cellMaker = pd.read_excel(marker_path)