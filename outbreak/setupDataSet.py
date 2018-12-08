import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np
import os
from datetime import datetime
from datetime import timedelta
from datetime import date
import string 

import random
import csv

parser = argparse.ArgumentParser(
    description='This script formats the raw data present in the Raw_Data directory to match \
    the date of the workshop. The script will output a full line list (./DataSet/Anonymised_line_list.csv),\
     a subseted line list (/DataSet/Line_list.csv") for the 40 analyzed samples, \
     and corresponding fasta files.')


parser.add_argument('current_date', metavar='current_date', nargs='+',
                    help='The date of the workshop. The most recent sample will taken \
                    1 day prior to this. It should be dd-mm-yyyy format')

parser.add_argument('--seed', dest='seed', action='store',type=int,
                    help='The radom seed used in generating the sample ids and names (if desired)')

# parser.add_argument('--names', dest='name_file', action='store',
#                     help='a csv file containing names to be used in generating an nonanymized csv file\
#                     This current is not available')

args = parser.parse_args()
if args.seed:
    seed=args.seed
else:
        seed=np.random.randint(100000)

print("Using seed: %s" % seed)
random.seed(seed)
np.random.seed(seed)

dir_path = os.path.dirname(os.path.realpath(__file__))

full_fasta_file = dir_path+ "/Raw_Data/full_sequence_set.fasta"
subset_fasta_file = dir_path+ "/DataSet/all_sample_sequences.fasta"
prior_to_workshop_file = dir_path+"/DataSet/prior_to_workshop_sequences.fasta"
raw_line_list_file = dir_path+"/Raw_Data/fullLineList.csv"
expected_sequences_file = dir_path+"/DataSet/expected_sequences.fasta"
###########################   process line list    #######################################
raw_line_list = pd.read_csv(raw_line_list_file, sep=",")
last_sequenced_case = 366 # case366 is the most recently sampled isolate - it is one of the sequenced samples

last_sample_day = raw_line_list.iloc[last_sequenced_case, 2]

raw_line_list["days ago"] = (last_sample_day - raw_line_list["sampleTime"])
sampling_dates_backwards = raw_line_list['days ago']

# Set the end date : ie the date of sampling for case 352
ddmmyyy=args.current_date[0].split("-")
end_date = datetime(int(ddmmyyy[0]),int(ddmmyyy[1]),int(ddmmyyy[2])) 

date_list2 = []

for i in sampling_dates_backwards:
    difference = timedelta(i)
    new_date = end_date - difference
    date_notime = new_date.strftime("%Y-%m-%d") #Get rid of time from the variable, as this is not usually included in line listings
    date_list2.append(date_notime)
    
df_backwards_time = pd.DataFrame(date_list2)
df_backwards_time.columns = ["Date of sampling"]


parsed_line_list = pd.DataFrame(df_backwards_time["Date of sampling"])
#parsed_line_list.columns="Date of sampling"
parsed_line_list["Id"] = raw_line_list["Id"] #So that the contacts can be given names


parsed_line_list["Contact"] = raw_line_list["parentId"]



#with_contacts = new_names.assign(Id = assigned.values) #Replacese "case0" with the appropriate name

parsed_line_list.columns = ["Date of sampling", "Case id", "Contact"]
lstoflst = []
clean_list = []

for l in parsed_line_list["Contact"][1:]: #Missing out first row which has a none option in the contact column
    my_string = str(l) #From a tuple
    thing = my_string.replace("'", "").replace(",", "").replace(")", "").replace("(", "")
    
    lstoflst.append(thing)

r = pd.Series(lstoflst)

parsed_line_list["Contact"][1:] = r

#tidy = with_contacts
death_list = [] 

case_fatality_rate = 70

for i in range(0, len(raw_line_list["Id"])):
    prob = np.random.uniform(0,100)
    if prob > case_fatality_rate:
        death_list.append("Alive")
    else:
        death_list.append("Dead")

death_series = pd.Series(death_list)

parsed_line_list["Outcome"] = death_series
raw_line_list["onset days ago"] = (last_sample_day - raw_line_list["onset"])
onset_dates_backwards = raw_line_list['onset days ago']

date_list3 = []

for i in onset_dates_backwards:
    
    difference = timedelta(i)
    new_date = end_date - difference
    
    date_notime = new_date.strftime("%Y-%m-%d") #Get rid of time from the variable, as this is not usually included in line listings
    
    date_list3.append(date_notime)
    
df_onset_time = pd.DataFrame(date_list3)

df_onset_time.columns = ["Date of onset"]
parsed_line_list["Date of onset"] = df_onset_time
parsed_line_list = parsed_line_list[[ "Case id", "Date of onset", "Date of sampling", "Contact", "Outcome"]]


present_loc2 = []
## Adding locations - could eventually be in raw line list
with open(dir_path+"/Raw_Data/location2.csv", 'r') as f:
    next(f)
    for l in f:
        toks = l.strip().split(",")
        present_loc2.append(toks[0])
loc_list = []
loc_b = []

for i in parsed_line_list["Case id"]:
    #print(i)
    if i in present_loc2:
        loc_list.append("Location B")
    else:
        loc_list.append("Location A")


    
parsed_line_list["Location"] = pd.Series(loc_list)
#true_list = tidy_now

#true_list.to_csv(dir_path+"/DataSet/true_list.csv", na_rep = "NA", index = False)


#dataframe = pd.read_csv(line_list, sep=",")

## Adding random id and lab id
id_list = []

for i in range(0, len(raw_line_list)):
    i = ''.join(random.choice(string.ascii_uppercase + string.digits) for i in range(6))
    id_list.append(i)

#print(id_list)

id_list = pd.DataFrame(id_list)
id_list.columns = ["Id number"]

parsed_line_list["Id number"] = id_list["Id number"]

parsed_line_list.columns = ["Case id", "Date of symptom onset", "Date of sampling", "Contact", "Outcome", "Location", "Id number"]

parsed_line_list = parsed_line_list[["Case id", "Id number", "Date of symptom onset", "Date of sampling", "Contact", "Outcome", "Location"]]

headers = ["Case Id", "Id", "Date of infection", "Date of sampling", "Contact", "Outcome", "Location"]

parsed_line_list.columns = headers



# if args.name_file:

#     namefile = pd.read_csv(dir_path+"/Raw_data/Name_files/"+args.name_file) #Makes a dataframe of all hobbit names
#     names = namefile.sample(n=len(raw_line_list)) #Takes a random sample of those names, n=number of cases. Output is a dataframe
#     new_names = names.reset_index(drop= True) #So that they can be matched by index rather than case number (which are arbitrary)
#     #parsed_line_list["First"]=
#     #parsed_line_list["Surname"]=




parsed_line_list.to_csv(dir_path+"/DataSet/Anonymised_line_list.csv", na_rep = "NA", index = False)







# Sampled cases
case_id_list = ['case10', 'case13', 'case14', 'case15', 'case16', 'case18', 'case19', 'case20', 'case21', 'case22', 'case23', 'case24', 'case26', 'case29', 'case32', 'case28', 'case39', 'case44', 'case45', 'case49', 'case67', 'case77', 'case81', 'case95', 'case98', 'case111', 'case121', 'case123', 'case126', 'case132', 'case142', 'case145', 'case149', 'case155', 'case178', 'case196', 'case260', 'case321', 'case352', 'case370','case366','case339']

# row_list = []

# with open(dir_path+"/DataSet/Anonymised_line_list.csv", 'r') as f:
#     next(f)
#     for l in f:
#         tokens = l.strip("\n").split(",")
#         case = tokens[0]
#         for i in case_id_list:
#             if case == i:
#                 row_list.append(tokens)


subset_samples=parsed_line_list.loc[parsed_line_list['Case Id'].isin(case_id_list)]
subset_samples=subset_samples.sort_values(by=["Date of sampling"])
lab_ids = []
options=[0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,2,3,3,3,5,7]
current_index=1
for i in range(0,len(case_id_list)):
    pick =np.random.randint(0,len(options)-1)
    lab_ids.append("PHR-"+str(current_index+options[pick]))
    current_index += options[pick]+1

subset_samples["Sample Id"]=lab_ids

# The contacts here should be PHR labels and only if the contact is in 
# the subset
subsetted_contacts=[]

for index, row in subset_samples.iterrows():
    if row["Contact"] in case_id_list:
        contact_lab_id = subset_samples.loc[subset_samples["Case Id"]==row["Contact"],"Sample Id"]
        subsetted_contacts.append(contact_lab_id.values[0])
    else:
        subsetted_contacts.append("Unknown")

subset_samples.drop(columns=["Contact"])
subset_samples["Contact"] = subsetted_contacts
column_order = ['Sample Id', 'Id', 'Date of infection', 'Date of sampling', 'Contact', 'Outcome', 'Location', 'Case Id']
subset_samples=subset_samples[column_order]

subset_samples.to_csv(dir_path+"/DataSet/Line_list.csv", na_rep = "NA", index = False)

###########################   process fasta files   #######################################


fasta_sequences = SeqIO.parse(open(full_fasta_file),'fasta')
sequenced = ["case260","case321","case352","case370","case366","case399"]
with open(subset_fasta_file,"w") as sampled_seqs:
    with open(prior_to_workshop_file,"w") as prior_to_workshop_seqs:
        with open(expected_sequences_file,"w") as expected_sequences:
            for sequence in SeqIO.parse(open(full_fasta_file),'fasta'):
    #             # Write the sampled sequences to file 
                if sequence.id in case_id_list:
                    metaData = subset_samples.loc[subset_samples['Case Id'] == sequence.id]
                    if sequence.id=="case366":
                        outbreak_C = [ record for record in SeqIO.parse(open(full_fasta_file),'fasta') if record.id=="OB_C"]
                        nucleotide_sequence=outbreak_C[0].seq
                    else:
                        nucleotide_sequence = sequence.seq
                    sampled_seqs.write(">{}|{}|{}\n{}\n".format(metaData['Sample Id'].values[0], metaData["Location"].values[0],metaData["Date of sampling"].values[0], nucleotide_sequence))
                    # Write the sampled sequences to file exp
                    if sequence.id not in sequenced:
                        prior_to_workshop_seqs.write(">{}|{}|{}\n{}\n".format(metaData['Sample Id'].values[0], metaData["Location"].values[0],metaData["Date of sampling"].values[0], nucleotide_sequence))
                    else:
                        expected_sequences.write(">{}|{}|{}\n{}\n".format(metaData['Sample Id'].values[0], metaData["Location"].values[0],metaData["Date of sampling"].values[0], nucleotide_sequence))

