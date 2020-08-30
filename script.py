import pandas as pd

# df = pd.read_csv('am_flows.csv')
# df[['ID1', 'AB_Time', 'BA_Time', 'AB_Flow', 'BA_Flow']].set_index('ID1').to_csv('times_and_flows.csv')
# try:
#     calculated = pd.read_csv('calculated_file.csv')
# except:
#     print("Calculated file missing")
#     exit(1)
# extras = list(set(df['ID1']) - set(calculated['ID1']))
# updated_truth = df[~df['ID1'].isin(extras)]
# reformatted = dict()
# for idx, row in calculated.iterrows():
#     if row['ID1'] not in reformatted:
#         reformatted[row['ID1']] = dict()
#     reformatted[row['ID1']]['AB_Time' if row['Is_AB'] else 'BA_Time'] = row['time']
# final_df = pd.DataFrame.from_dict(reformatted, orient='index')
# final_df.fillna(0, inplace=True)
# final_df.index.name = 'ID1'
# updated_truth = df.set_index('ID1')
# updated_truth.index.name = 'ID1'
# updated_truth.fillna(0, inplace=True)
# error = (updated_truth[['AB_Time', 'BA_Time']] - final_df).abs()
# for col in error.columns:
#     print(error[col].describe())

df = pd.read_csv('../../Downloads/am_inputMtx.csv', header=False)
print(df.columns)