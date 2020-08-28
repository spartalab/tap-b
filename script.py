import pandas as pd

df = pd.read_csv('am_flows.csv')
df[['ID1','AB_Time', 'BA_Time', 'AB_Flow', 'BA_Flow']].set_index('ID1').to_csv('times_and_flows.csv')