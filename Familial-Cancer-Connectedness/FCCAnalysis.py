import os
import csv
import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import time
import json
from scipy.stats import wilcoxon


DOMAIN = 'https://labs.tib.eu/sdm/clarify_kpi3/'

KG = 'https://labs.tib.eu/sdm/clarify-kg-7-1/sparql'


dtype = [('group', '|S256'), ('CancerType', '|S256'), ('Biomarker', '|S256'), ('value', float)]
dtypeNew = [('group', '|S256'), ('CancerType', '|S256'), ('value', float)]

def current_milli_time():
    return round(time.time() * 1000)

def totalPopulation(input_data,endpoint):
	query="""SELECT count(DISTINCT ?ehr1)  as ?num \n
	    WHERE {
		?ehr a <http://clarify2020.eu/vocab/LCPatient>. 
		?ehr <http://clarify2020.eu/vocab/has_LC_SLCG_ID> ?ehr1 . }"""
	sparql = SPARQLWrapper(endpoint)
	sparql.setQuery(query)
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()
	return int(results["results"]["bindings"][0]["num"]["value"])

def query_generation(input_data,endpoint):
	where_clause = {"Biomarker":"OPTIONAL {?ehr1 <http://clarify2020.eu/vocab/hasBio> ?Biomarker.}.",
					"Relapse":"OPTIONAL {?ehr1 <http://clarify2020.eu/vocab/hasProgressionRelapse> ?Relapse.}.",
					"Stages":"OPTIONAL {?ehr1 <http://clarify2020.eu/vocab/hasDiagnosis> ?o1 . \n ?o1 <http://clarify2020.eu/vocab/hasDiagnosisStage> ?Stages.}.",
					"Tumor":"OPTIONAL {?ehr1 <http://clarify2020.eu/vocab/hasTumorHistology> ?Tumor.}.",
					"FamilyDegree":"""?ehr1 <http://clarify2020.eu/vocab/hasFamilyHistory> ?o .\n ?o <http://clarify2020.eu/vocab/familyRelationDegree> ?familyType .\n ?o  <http://clarify2020.eu/vocab/hasFamilyCancerType> ?CancerType .""",
					"FamilyRelationship":"""?ehr1 <http://clarify2020.eu/vocab/hasFamilyHistory> ?o .\n ?o <http://clarify2020.eu/vocab/familyType> ?familyType .\n ?o  <http://clarify2020.eu/vocab/hasFamilyCancerType> ?CancerType .""",}
	select_clause = {"Biomarker":"?Biomarker",
					"Relapse":"?Relapse",
					"Stages":"?Stages",
					"Tumor":"?Tumor"}
	query_select_clause = "SELECT DISTINCT ?ehr1 ?familyType ?CancerType"
	query_where_clause="""WHERE {
		?ehr a <http://clarify2020.eu/vocab/LCPatient>. 
		?ehr <http://clarify2020.eu/vocab/has_LC_SLCG_ID> ?ehr1 . \n"""
	# print(input_data)
	if "Age" in input_data["Input"]["IndependentVariables"]:
		query_where_clause = query_where_clause + "?ehr1  <http://clarify2020.eu/vocab/age> ?age.\n"
		query_select_clause = query_select_clause + "  ?age "
	if "Gender" in input_data["Input"]["IndependentVariables"]:
		query_where_clause = query_where_clause + "?ehr1 <http://clarify2020.eu/vocab/sex> ?gender. FILTER (regex(?gender,\"" + input_data["Input"]["IndependentVariables"]["Gender"] + "\"))\n"
	if "SmokingHabits" in input_data["Input"]["IndependentVariables"]:
		query_where_clause = query_where_clause + "?ehr1 <http://clarify2020.eu/vocab/hasSmokingHabit> ?smoking. FILTER (regex(?smoking,\"" + input_data["Input"]["IndependentVariables"]["SmokingHabits"] + "\"))\n"
	if "FamilyType" in input_data["Input"]["IndependentVariables"].keys() and "FamilyDegree" == input_data["Input"]["IndependentVariables"]["FamilyType"]:
		query_where_clause = query_where_clause + where_clause["FamilyDegree"] + "\n"
	if "FamilyType" in input_data["Input"]["IndependentVariables"].keys() and "FamilyRelationship" == input_data["Input"]["IndependentVariables"]["FamilyType"]:
		query_where_clause = query_where_clause + where_clause["FamilyRelationship"] + "\n"

	for variable in input_data["Input"]["DependentVariables"]:
		if variable != "CancerType":
			query_select_clause = query_select_clause + select_clause[variable] + " "
			query_where_clause= query_where_clause + where_clause[variable] + " \n"

	query_where_clause = query_where_clause[:-1] + "}"
	sparql_query = query_select_clause + " " + query_where_clause
	# print(sparql_query)

	sparql = SPARQLWrapper(endpoint)
	sparql.setQuery(sparql_query)
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()
	return results["results"]["bindings"]


def sort_dic(dic):
	temp_dic = {}
	sorted_names = sorted(dic.keys(), key=str.lower)
	for name in sorted_names:
		temp_dic[name] = ""
	return temp_dic


def dependent_values(data,variables):
	unique_values = {}
	for variable in variables:
		unique_values[variable] = {}
		for row in data:
			if variable in row:
				unique_value = row[variable]["value"].replace("http://clarify2020.eu/entity/","")
				if unique_value not in unique_values[variable]:
					unique_values[variable][unique_value] = ""
		unique_values[variable] = sort_dic(unique_values[variable])
	return unique_values


def divided_by_age(data,age):
	divided_data = {"older":[],"younger":[]}
	for row in data:
		if int(row["age"]["value"]) <= age:
			divided_data["younger"].append(row)
		else:
			divided_data["older"].append(row)
	return divided_data


def projected_data(data,variable,values):
	pro_data = {}
	#if variable == "familyType":
	#	pro_data["No_Relative"] = []
	for value in values:
		pro_data[value] = []
		for row in data:

			if variable in row:
				if value in row[variable]["value"]:
					pro_data[value].append(row["ehr1"]["value"])
			#else:
			#	if variable == "familyType":
			#		pro_data["No_Relative"].append(row["ehr1"]["value"])
	return(pro_data)


def fcc_ternary_set(list1, list2, list3):
	set1 = set(list1)
	set2 = set(list2)
	set3 = set(list3)
	intersection = set1.intersection(set2.intersection(set3))
	union = (len(set1) + len(set2) + len(set3)) - len(intersection)
	if len(intersection) == 0.0 and union == 0.0:
		return 0.0
	else:
		return float(len(intersection)) / union


def fcc_binary(list1, list2):
	"""A function for finding the similarity between two binary vectors"""
	intersection = []
	list3 = list(list2)
	for value in list1:
		if value in list3:
			intersection.append(value)
			list3.remove(value)
	union = (len(list1) + len(list2)) - len(intersection)
	if len(intersection) == 0.0 and union == 0.0:
		return 0.0
	else:
		return float(len(intersection)) / union

def fcc_ternary(list1, list2, list3):
	"""A function for finding the similarity between three  vectors"""
	intersection = []
	list41 = list(list2)
	list42 = list(list3)
	for value in list1:
		if (value in list41) and (value in list42):
			intersection.append(value)
			list41.remove(value)
			list42.remove(value)
	union = (len(list1) + len(list2) + len(list3)) - len(intersection)
	if len(intersection) == 0.0 and union == 0.0:
		return 0.0
	else:
		return float(len(intersection)) / union

def overlap_ternary(list1, list2, list3):
	intersection = []
	list41 = list(list2)
	list42 = list(list3)
	for value in list1:
		if (value in list41) and (value in list42):
			intersection.append(value)
			list41.remove(value)
			list42.remove(value)
	union = min(len(list1), len(list2), len(list3))
	if len(intersection) == 0.0 and union == 0.0:
		return 0.0
	else:
		return float(len(intersection)) / union


def average(df):
	groups = np.unique(df['group'])
	variables = np.unique(df['variable'])

	df_avgs = np.empty(shape=0, dtype=dtype)

	for group in groups:
		df_group = df[(df['group'] == group) & (df['variable'] != b'No_Relative')]
		avg = np.average(df_group['value'])
		res = np.array([(group.decode('utf-8'), 'COLUMN AVERAGE', avg)], dtype=dtype)
		df_avgs = np.append(df_avgs, res, axis=0)

	for variable in variables:
		df_variable = df[df['variable'] == variable]
		avg = np.average(df_variable['value'])
		res = np.array([(str('ROW AVERAGE'), str(variable.decode('utf-8')), avg)], dtype=dtype)
		df_avgs = np.append(df_avgs, res, axis=0)

	return np.append(df, df_avgs, axis=0)


def fccMultiset(input_data):
	if not os.path.exists("static"):
		os.mkdir("static")
	random_id = str(current_milli_time())

	results = query_generation(input_data, KG)

	if "FamilyType" in input_data["Input"]["IndependentVariables"].keys() and input_data["Input"]["IndependentVariables"]["FamilyType"] == "FamilyRelationship":
		paramFamily = {"Father":"", "Mother":"", "Brother":"", "Sister":"", "Daugther":"", "Son":"", "Uncle":"", "Nephew":"","Grandfather":"", "Grandmother":"", "Aunt":"", "Niece":"", "Female_Cousin":"", "Male_Cousin":"", "Great_Grandmother":"", "Half_Brother":"", "Great_Grandfather":"", "Grandson":"", "Granddaugther":""}
	else:
		paramFamily = {"Male_First_Degree":"", "Female_First_Degree":"", "Male_Second_Degree":"", "Female_Second_Degree":"", "Male_Third_Degree":"", "Female_Third_Degree":""}

	if "Gender" in input_data["Input"]["IndependentVariables"].keys():
		gender=input_data["Input"]["IndependentVariables"]["Gender"]
	else:
		gender=""

	metric=input_data["Input"]["IndependentVariables"]["Metric"]
	
	if "SmokingHabits" in input_data["Input"]["IndependentVariables"].keys():
		smoking=input_data["Input"]["IndependentVariables"]["SmokingHabits"]
	else:
		smoking=""

	dic_file = {}
	dic_file["populationSizes"]= {}
	dic_file["populationSizes"]["cohortSize"] = len(set([row["ehr1"]["value"] for row in results]))
	dic_file["populationSizes"]["allPopulationSize"] = totalPopulation(input_data, KG)
	dependent_variables = dependent_values(results,input_data["Input"]["DependentVariables"])
	if "Age" in input_data["Input"]["IndependentVariables"]:
		divided_data = divided_by_age(results,int(input_data["Input"]["IndependentVariables"]["Age"]))
		older_family_data = projected_data(divided_data["older"],"familyType",paramFamily)
		younger_family_data = projected_data(divided_data["younger"],"familyType",paramFamily)
		dic_file["populationSizes"]["cohortOldPatients"] = len(set([row1["ehr1"]["value"] for row1 in divided_data["older"]]))
		dic_file["populationSizes"]["cohortYoungPatients"] = len(set([row1["ehr1"]["value"] for row1 in divided_data["younger"]]))
		older_data_cancerType = projected_data(divided_data["older"],"CancerType",dependent_variables["CancerType"])
		younger_data_cancerType = projected_data(divided_data["younger"],"CancerType",dependent_variables["CancerType"])
		older_data_biomarker = projected_data(divided_data["older"],"Biomarker",dependent_variables["Biomarker"])
		younger_data_biomarker = projected_data(divided_data["younger"],"Biomarker",dependent_variables["Biomarker"])
		older_cancerType_biomaker_distribution = []
		younger_cancerType_biomaker_distribution = []
		df_older_patients_ALK = np.empty(shape=0, dtype=dtypeNew)
		df_older_patients_EGFR = np.empty(shape=0, dtype=dtypeNew)
		df_older_patients_Others = np.empty(shape=0, dtype=dtypeNew)
		df_older_patients = np.empty(shape=0, dtype=dtype)
		for var1 in older_data_cancerType:
			for var2 in older_data_biomarker:
				for member in older_family_data:
					if "Metric" in input_data["Input"]["IndependentVariables"].keys() and input_data["Input"]["IndependentVariables"]["Metric"] == "Overlap":
						value = overlap_ternary(np.array(older_data_cancerType[var1]), np.array(older_data_biomarker[var2]), np.array(older_family_data[member]))
					elif "Metric" in input_data["Input"]["IndependentVariables"].keys() and input_data["Input"]["IndependentVariables"]["Metric"] == "FCCSet":
						value = fcc_ternary_set(np.array(older_data_cancerType[var1]), np.array(older_data_biomarker[var2]), np.array(older_family_data[member]))
					else:
						value = fcc_ternary(np.array(older_data_cancerType[var1]), np.array(older_data_biomarker[var2]), np.array(older_family_data[member]))
					older_cancerType_biomaker_distribution.append(value)
					res = np.array([(str(var1), str(var2), str(member), value)], dtype=dtype)
					df_older_patients = np.append(df_older_patients, res, axis=0)
					if str(var2)=="ALK":
						resALK = np.array([(str(member),str(var1), value)], dtype=dtypeNew)
						df_older_patients_ALK = np.append(df_older_patients_ALK, resALK, axis=0)
					elif str(var2)=="EGFR":
						resEGFR = np.array([(str(member),str(var1), value)], dtype=dtypeNew)
						df_older_patients_EGFR = np.append(df_older_patients_EGFR, resEGFR, axis=0)
					else:
						resOthers = np.array([(str(member),str(var1), value)], dtype=dtypeNew)
						df_older_patients_Others = np.append(df_older_patients_Others, resOthers, axis=0)
			#df_older_patients = average(df_older_patients)


		
		df_ALKOld = pd.DataFrame(df_older_patients_ALK)
		df_ALKOld = df_ALKOld.pivot('group', 'CancerType', 'value')
		x_labels_alk_older = [lbl.decode('utf8') for lbl in np.unique(df_older_patients_ALK['CancerType'])]
		y_labels_alk_older = [lbl.decode('utf8') for lbl in np.unique(df_older_patients_ALK['group'])]
		df_EGFROld = pd.DataFrame(df_older_patients_EGFR)
		df_EGFROld = df_EGFROld.pivot('group', 'CancerType', 'value')
		x_labels_egfr_older = [lbl.decode('utf8') for lbl in np.unique(df_older_patients_EGFR['CancerType'])]
		y_labels_egfr_older = [lbl.decode('utf8') for lbl in np.unique(df_older_patients_EGFR['group'])]
		fig = plt.figure(figsize=(22, 12))

		plt.subplot(121)  
		plt.title('Heatmap EGFR Older ' + gender + '  ' + smoking + ' Patients--' + metric)
		sns.heatmap(df_EGFROld, xticklabels=x_labels_egfr_older, yticklabels=y_labels_egfr_older, fmt='.2f', square=True, cmap='Greens')

		plt.subplot(122)   
		plt.title('Heatmap ALK Older ' + gender + '  ' + smoking + ' Patients --' + metric)
		sns.heatmap(df_ALKOld, xticklabels=x_labels_alk_older, yticklabels=y_labels_alk_older, fmt='.2f', square=True, cmap='Blues')

		plt.tight_layout()
		plt.savefig("static/" + "OlderEGFR-OlderALK" + gender + smoking + metric + random_id + ".png")


		df_younger_patients = np.empty(shape=0, dtype=dtype)
		df_younger_patients_ALK = np.empty(shape=0, dtype=dtypeNew)
		df_younger_patients_EGFR = np.empty(shape=0, dtype=dtypeNew)
		df_younger_patients_Others = np.empty(shape=0, dtype=dtypeNew)
		for var1 in younger_data_cancerType:
			for var2 in younger_data_biomarker:
				for member in younger_family_data:
					if "Metric" in input_data["Input"]["IndependentVariables"].keys() and input_data["Input"]["IndependentVariables"]["Metric"] == "Overlap":
						value = overlap_ternary(np.array(younger_data_cancerType[var1]), np.array(younger_data_biomarker[var2]), np.array(younger_family_data[member]))
					elif "Metric" in input_data["Input"]["IndependentVariables"].keys() and input_data["Input"]["IndependentVariables"]["Metric"] == "FCCSet":
						value = fcc_ternary_set(np.array(younger_data_cancerType[var1]), np.array(younger_data_biomarker[var2]), np.array(younger_family_data[member]))
					else:
						value = fcc_ternary(np.array(younger_data_cancerType[var1]), np.array(younger_data_biomarker[var2]), np.array(younger_family_data[member]))
						if (str(var2)=="EGFR" and str(var1)=="Breast" and str(member)=="Male_First_Degree"):
							print(np.array(younger_data_cancerType[var1]))
							print(np.array(younger_data_biomarker[var2]))
							print(np.array(younger_family_data[member]))
							print(value)
					younger_cancerType_biomaker_distribution.append(value)
					res = np.array([(str(member),str(var1), str(var2), value)], dtype=dtype)
					df_younger_patients = np.append(df_younger_patients, res, axis=0)

					if str(var2)=="ALK":
						resALK = np.array([(str(member),str(var1), value)], dtype=dtypeNew)
						df_younger_patients_ALK = np.append(df_younger_patients_ALK, resALK, axis=0)
					elif str(var2)=="EGFR":
						resEGFR = np.array([(str(member),str(var1), value)], dtype=dtypeNew)
						df_younger_patients_EGFR = np.append(df_younger_patients_EGFR, resEGFR, axis=0)
					else:
						resOthers = np.array([(str(member),str(var1), value)], dtype=dtypeNew)
						df_younger_patients_Others = np.append(df_younger_patients_Others, resOthers, axis=0)

			
	
		df_ALK1 = pd.DataFrame(df_younger_patients_ALK)
		df_ALK1 = df_ALK1.pivot('group', 'CancerType', 'value')
		x_labels_alk_younger = [lbl.decode('utf8') for lbl in np.unique(df_younger_patients_ALK['CancerType'])]
		y_labels_alk_younger = [lbl.decode('utf8') for lbl in np.unique(df_younger_patients_ALK['group'])]
		df_EGFR2 = pd.DataFrame(df_younger_patients_EGFR)
		df_EGFR2 = df_EGFR2.pivot('group', 'CancerType', 'value')
		x_labels_egfr_younger = [lbl.decode('utf8') for lbl in np.unique(df_younger_patients_EGFR['CancerType'])]
		y_labels_egfr_younger = [lbl.decode('utf8') for lbl in np.unique(df_younger_patients_EGFR['group'])]
		fig = plt.figure(figsize=(22, 12))

		plt.subplot(121)   
		plt.title('Heatmap EGFR Younger ' + gender + '  ' + smoking + ' Patients--' + metric)
		sns.heatmap(df_EGFR2, xticklabels=x_labels_egfr_younger, yticklabels=y_labels_egfr_younger, fmt='.2f', square=True, cmap='OrRd')

		plt.subplot(122)   
		plt.title('Heatmap ALK Younger ' + gender + '  ' + smoking + ' Patients--' + metric)
		sns.heatmap(df_ALK1, xticklabels=x_labels_alk_younger, yticklabels=y_labels_alk_younger,  fmt='.2f', square=True, cmap='BuPu')

		plt.tight_layout()
		plt.savefig("static/" + "YoungerEGFR-YoungerALK" + gender + smoking + metric + random_id + ".png")


						
		#df_younger_patients = average(df_older_patients)

	

		with open("static/" + "CancerType_Biomarker" + random_id + "_older.csv", "w") as temp_csv:
			writer = csv.writer(temp_csv, quoting=csv.QUOTE_ALL)
			writer.writerow(["group", "CancerType", "Biomarker", "value"])
			for row in df_older_patients:
				writer.writerow([row['group'].decode('utf-8'), row["CancerType"].decode('utf-8'), row["Biomarker"].decode('utf-8'),row['value']])
		with open("static/" + "CancerType_Biomarker" + random_id + "_younger.csv", "w") as temp_csv:
			writer = csv.writer(temp_csv, quoting=csv.QUOTE_ALL)
			writer.writerow(["group", "CancerType","Biomarker", "value"])
			for row in df_younger_patients:
				writer.writerow([row['group'].decode('utf-8'), row["CancerType"].decode('utf-8'), row["Biomarker"].decode('utf-8'), row['value']])
		dic_file["CancerType_Biomarker"] = {}
		dic_file["CancerType_Biomarker"]["files"] = [DOMAIN+"static/" + "CancerType_Biomarker" + random_id + "_older.csv", DOMAIN+"static/" + "CancerType_Biomarker" + random_id + "_younger.csv"]

			#stat, p = wilcoxon(older_cancerType_biomaker_distribution, younger_cancerType_biomaker_distribution,alternative='two-sided')
			# print('Statistics=%.3f, p=%.9f' % (stat, p))
		dic_file["CancerType_Biomarker"]["p_value"] = []

	else:
		family_data = projected_data(results,"familyType",paramFamily)
		cancerType_data = projected_data(results,"CancerType",dependent_variables["CancerType"])
		biomarker_data  = projected_data(results,"Biomarker",dependent_variables["Biomarker"])
		df_patients = np.empty(shape=0, dtype=dtype)
		df_patients_ALK = np.empty(shape=0, dtype=dtypeNew)
		df_patients_EGFR = np.empty(shape=0, dtype=dtypeNew)
		df_patients_Others = np.empty(shape=0, dtype=dtypeNew)
		for var1 in cancerType_data:
			for var2 in biomarker_data:
				for member in family_data:
					if "Metric" in input_data["Input"]["IndependentVariables"].keys() and input_data["Input"]["IndependentVariables"]["Metric"] == "Overlap":
						value = overlap_ternary(np.array(cancerType_data[var1]), np.array(biomarker_data[var2]), np.array(family_data[member]))
					elif "Metric" in input_data["Input"]["IndependentVariables"].keys() and input_data["Input"]["IndependentVariables"]["Metric"] == "FCCSet":
						value = fcc_ternary_set(np.array(cancerType_data[var1]), np.array(biomarker_data[var2]), np.array(family_data[member]))
					else:
						value = fcc_ternary(np.array(cancerType_data[var1]), np.array(biomarker_data[var2]), np.array(family_data[member]))
					res = np.array([(str(member),str(var1), str(var2), value)], dtype=dtype)
					df_patients = np.append(df_patients, res, axis=0)
					if str(var2)=="ALK":
						resALK = np.array([(str(member),str(var1), value)], dtype=dtypeNew)
						df_patients_ALK = np.append(df_patients_ALK, resALK, axis=0)
					elif str(var2)=="EGFR":
						resEGFR = np.array([(str(member),str(var1), value)], dtype=dtypeNew)
						df_patients_EGFR = np.append(df_patients_EGFR, resEGFR, axis=0)
					else:
						resOthers = np.array([(str(member),str(var1), value)], dtype=dtypeNew)
						df_patients_Others = np.append(df_patients_Others, resOthers, axis=0)
		

		df_ALK = pd.DataFrame(df_patients_ALK)
		df_ALK = df_ALK.pivot('group', 'CancerType', 'value')
		x_labels_alk = [lbl.decode('utf8') for lbl in np.unique(df_patients_ALK['CancerType'])]
		y_labels_alk = [lbl.decode('utf8') for lbl in np.unique(df_patients_ALK['group'])]
		df_EGFR = pd.DataFrame(df_patients_EGFR)
		df_EGFR = df_EGFR.pivot('group', 'CancerType', 'value')
		x_labels_egfr = [lbl.decode('utf8') for lbl in np.unique(df_patients_EGFR['CancerType'])]
		y_labels_egfr = [lbl.decode('utf8') for lbl in np.unique(df_patients_EGFR['group'])]
		fig = plt.figure(figsize=(22, 12))

		plt.subplot(121)   
		plt.title('Heatmap EGFR ' + gender + '  ' + smoking + ' Patients--' + metric)
		sns.heatmap(df_EGFR, xticklabels=x_labels_egfr, yticklabels=y_labels_egfr, fmt='.2f', square=True, cmap='OrRd')

		plt.subplot(122)   
		plt.title('Heatmap ALK ' + gender + '  ' + smoking + ' Patients--' + metric)
		sns.heatmap(df_ALK, xticklabels=x_labels_alk, yticklabels=y_labels_alk, fmt='.2f', square=True, cmap='BuPu')

		plt.tight_layout()
		plt.savefig("static/" + "EGFR-ALK" + gender + smoking + metric + random_id + ".png")


		#df_patients = average(df_patients)

		with open("static/" + "CancerType_Biomarker" + random_id + ".csv", "w") as temp_csv:
			writer = csv.writer(temp_csv, quoting=csv.QUOTE_ALL , lineterminator='\n')
			writer.writerow(["group", "CancerType","Biomarker", "value"])
			for row in df_patients:
				writer.writerow([row['group'].decode('utf-8'), row["CancerType"].decode('utf-8'), row["Biomarker"].decode('utf-8'), row['value']])
		dic_file["CancerType_Biomarker"] = {}
		dic_file["CancerType_Biomarker"]["files"] = [DOMAIN+"static/" +  "CancerType_Biomarker" + random_id + ".csv"]
		dic_file["CancerType_Biomarker"]["p_value"] = []
	return dic_file


if __name__ == '__main__':
	# query_generation("input.json","https://labs.tib.eu/sdm/clarify-kg-6-1/sparql")
	with open("input.json", "r") as input_file:
		input_data = json.load(input_file)

	print(fccMultiset(input_data))
