import pandas as pd
import itertools

column_names = ["_HNWDAINC123", "_HBWDAINC123", "_HNWDAINC4andNHBDANONWRK", "_HBWDAINC4andNHBDAWRK", 
				"_HNWSRINC123", "_HBWSRINC123", "_HNWSRINC4andNHBSRNONWRK", "_HBWSRINC4andNHBSRWRK",
				"_MedTruck", "_HvyTruck"]
AB_PREFIX = 'AB_Flow'
BA_PREFIX = 'BA_Flow'

AB_Column_names = [AB_PREFIX+name for name in column_names]
BA_Column_names = [BA_PREFIX+name for name in column_names]

translator = pd.read_csv('__translator.csv')
truth = pd.read_csv('am_flows.csv').drop(columns=['AB_Time', 'BA_Time'])

ids = set(truth['ID1'])
unaccounted_total_flow = 0
tapb_only_links = 0

def process_total_flow(line):
	return float(line.split(' ')[1])

def process_class_flows(line, prefix):
	return {prefix + column_names[idx]: val for idx, val in enumerate(map(lambda x: float(x), line.split(' ')[:-1]))}

data = dict.fromkeys(ids)
with open('low_gap_flows_tapb.txt') as tapb_flows:
	for num, line in enumerate(tapb_flows):
		idx = num // 2
		try:
			translation = translator[translator['ID_tapb'] == idx].iloc[0]
		except:
			print("Couldn't find idx ", idx, " in translation table")
			continue
		is_AB = translation['AB']
		prefix = AB_PREFIX if is_AB else BA_PREFIX
		try:
			tap_output_row = truth[truth['ID1'] == translation['ID_TransCAD']].iloc[0]
		except:
			print("Couldn't find ID in am_flows", translation['ID_TransCAD'])
			tapb_only_links += 1
			if num % 2 == 0:
				print("Total flow for link: ", process_total_flow(line))
				unaccounted_total_flow += process_total_flow(line)
			else:
				print("Class flows for link: ", process_class_flows(line, prefix))
			continue
		try: 
			ids.remove(translation['ID_TransCAD'])
		except:
			pass
		if data[translation['ID_TransCAD']] is None:
			data[translation['ID_TransCAD']] = dict()
		
		row = data[translation['ID_TransCAD']]
		if num % 2 == 0:
			total_flow = process_total_flow(line)
			row[prefix] = total_flow
		else:
			class_flows = process_class_flows(line, prefix)
			row.update(class_flows)

data = {k: v for k, v in data.items() if v is not None}
tap_output = pd.DataFrame.from_dict(data, orient='index')
tap_output.fillna(0, inplace=True)
tap_output['Tot_Flow'] = tap_output[AB_PREFIX] + tap_output[BA_PREFIX]
tap_output.index.name = 'ID1'
print("Number of Links in am_flows but not in tapB output: ", len(ids))
print("Number of links in tapB output/translation but not in am_flows", tapb_only_links)
print("Total flow in lines in tapB output but not in am_flows ", unaccounted_total_flow)

truth.set_index('ID1', inplace=True)
tap_output.fillna(0, inplace=True)

tap_output.to_csv("translated_tap_flows.csv")