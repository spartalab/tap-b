import pandas as pd
import matplotlib.pyplot as plt

tap_output = pd.read_csv('translated_tap_flows.csv')
truth = pd.read_csv('am_flows.csv').drop(columns=['AB_Time', 'BA_Time']).fillna(0)
extras = list(set(truth['ID1']) - set(tap_output['ID1']))
updated_truth = truth[~truth["ID1"].isin(extras)]

error = tap_output.set_index('ID1').subtract(updated_truth.set_index('ID1'), axis='index').abs()
print(len(error))
for col in error.columns: 
    print(error[col].describe()) 

error.hist(column=["Tot_Flow"], bins=2000, color='#607c8e')
plt.show()