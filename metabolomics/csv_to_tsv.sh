# sed -E 's/("([^"]*)")?,/\2\t/g' data/SOURCE2022_Serum_Data_MZM13_N2_filtered_calour.csv > data/SOURCE2022_Serum_Data_MZM13_N2_filtered_calour.tsv
# sed -E 's/("([^"]*)")?,/\2\t/g' data/SOURCE2022_Stool_Data_MZM13_N2_filtered_calour.csv > data/SOURCE2022_Stool_Data_MZM13_N2_filtered_calour.tsv


sed -E 's/("([^"]*)")?,/\2\t/g' data/SOURCE2022_Serum_Data_MZM13_N2_filtered_v2.csv > data/SOURCE2022_Serum_Data_MZM13_N2_filtered_v2.tsv
sed -E 's/("([^"]*)")?,/\2\t/g' data/SOURCE2022_Stool_Data_MZM13_N2_filtered_v2.csv > data/SOURCE2022_tool_Data_MZM13_N2_filtered_v2.tsv

