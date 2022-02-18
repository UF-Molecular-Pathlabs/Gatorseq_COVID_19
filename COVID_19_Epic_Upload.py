import requests 
from requests.exceptions import HTTPError 
import time 
import ast
import sys
import yaml
import numpy
import pandas as pd
import hl7
import os
import datetime
from shutil import copyfile
from shutil import move
import csv
import traceback
import mysql.connector
from filelock import FileLock
from pathlib import Path
import os.path
import re
import math
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt

print("Run start time: ", str(datetime.datetime.now()) + "\n")

script_path = os.path.dirname(os.path.abspath( __file__ ))
parent_path = os.path.abspath(os.path.join(script_path, '..'))
comment_files = {"MCDC1RNA" : parent_path + "/COVID_19/Comments_CDC_Clinical_MCDC1RNA.txt", "MCDC1DIRECTPROT" : parent_path + "/COVID_19/Comments_CDC_Screening_MCDC1DIRECTPROT.txt"}

CONFIG_FILE = parent_path +"/linux_gatorseq.config.yaml"
config_dict=dict()
with open(CONFIG_FILE, 'r') as stream:
    try:
        config_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

CODE_ENV=script_path.split('/')[-3]
USER_NAME=os.environ['USER']

def replace_env(strname):
    strname=strname.replace("USER_NAME",USER_NAME).replace("CODE_ENV",CODE_ENV)
    return strname

COVID_19_EPIC_UPLOAD_TABLE = replace_env(config_dict['COVID_19_EPIC_UPLOAD_TABLE'])
COVID_19_TEST_INPUT_FOLDER = replace_env(config_dict['COVID_19_TEST_INPUT_FOLDER']) #'G:\DRL\Molecular\Assays\PGX\PGX_Beaker_Interface' 
COVID_19_TEST_SAMPLE_LOG = replace_env(config_dict['COVID_19_TEST_SAMPLE_LOG'])
CONFIG_TOKENS_FILE = parent_path  + "/" + config_dict['CONFIG_TOKENS_FILE']
MIRTH_GATORSEQ = config_dict['MIRTH_GATORSEQ']
if CODE_ENV=='ProdEnv':
    MIRTH_GATORSEQ += '/PROD'
    COVID_19_EPIC_UPLOAD_TABLE = replace_env(config_dict['COVID_19_EPIC_UPLOAD_TABLE_PROD'])
    COVID_19_TEST_INPUT_FOLDER = replace_env(config_dict['COVID_19_TEST_INPUT_FOLDER_PROD']) 
else:
    MIRTH_GATORSEQ += '/TEST'



config_token_dict=dict()
with open(CONFIG_TOKENS_FILE, 'r') as stream:
    try:
        config_token_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

if CODE_ENV == "ProdEnv":
    MYSQL_HOST = config_dict['PROD_MYSQL_HOST']
    MYSQL_USERNAME = config_dict['PROD_MYSQL_USERNAME']
    MYSQL_PASSWORD = config_token_dict['PROD_MYSQL_PASSWORD']
    MYSQL_DATABASE = config_dict['PROD_MYSQL_DATABASE']

if CODE_ENV == "DevEnv":
    MYSQL_HOST = config_dict['DEV_MYSQL_HOST']
    MYSQL_USERNAME = config_dict['DEV_MYSQL_USERNAME']
    MYSQL_PASSWORD = config_token_dict['DEV_MYSQL_PASSWORD']
    MYSQL_DATABASE = config_dict['DEV_MYSQL_DATABASE']


def check_folders_exist():
    #if not os.path.isfile(GATOR_SEQ_SAMPLE_INPUT_FILE):
        #sys.exit("ERROR: Does not have access to following folder: " + GATOR_SEQ_SAMPLE_INPUT_FILE + "\n") 

    if not os.path.isdir(MIRTH_GATORSEQ):
        sys.exit("ERROR: Does not have access to following folder: " + MIRTH_GATORSEQ + "\n") 

    if not os.path.isdir(COVID_19_TEST_INPUT_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + COVID_19_TEST_INPUT_FOLDER + "\n") 


check_folders_exist()
def format_for_unity(x):
    if(x<10):
        return "0"+str(x)
    else:
        return str(x)

def get_current_formatted_date():
    currentDT = datetime.datetime.now()
    if currentDT:
        data = format_for_unity(currentDT.year) +format_for_unity(currentDT.month) +format_for_unity(currentDT.month)+format_for_unity(currentDT.hour)+format_for_unity(currentDT.minute)+format_for_unity(currentDT.second)
        
        return data
    return str(currentDT)

def get_ticks(dt):
    return (dt - datetime.datetime(1, 1, 1)).total_seconds() * 10000000

class hl7update:
    
    def __init__(self, h):
        self.h = h

    def update_msh_segment(self):
        if self.h and self.h['MSH']:
            for msh_segment in self.h['MSH']:
                if msh_segment:
                    msh_segment[7] = get_current_formatted_date()
                    msh_segment[8] = ''
                    msh_segment[9][0][0] = 'ORU'
                    msh_segment[9][0][1] = 'R01'
                    msh_segment[10] = get_ticks(datetime.datetime.now())

    def update_orc_segment(self):
        if self.h and self.h['ORC']:
            for orc_segment in self.h['ORC']:
                orc_segment[1] = 'RE'

    def update_obr_segment(self):
        if self.h and self.h['OBR']:
            for obr_segment in self.h['OBR']:
                obr_segment[22] = get_current_formatted_date()
                obr_segment[25] = 'P'
                obr_segment[27] = '^^^^^R^^'

    def update_obx_segment(self):
        if self.h and self.h['OBX']:
            for obx_segment in self.h['OBX']:
                obx_segment[2] = 'ST'
                obx_segment[11] = 'P'
                obx_segment[14] = get_current_formatted_date()
                if(len(obx_segment)==19):
                    obx_segment.append(obx_segment[14])
                elif(len(obx_segment)>=19):
                    obx_segment[19] = obx_segment[14]
   
    def update_obx_seg_containing_gene(self, sample):
        updates = 0
        temp_obx = self.h[:]
        l = len(self.h)
        for i in range(l):
            del temp_obx[l-i-1]
        new_obx_index = 1
        for obxSegment in self.h['OBX']:
            if obxSegment[3][0][1][0] == "SARS-COV-2, NAA":
                obxSegment[5][0] = sample.result
                obxSegment[1] = new_obx_index
                obxSegment[18][0] = "Genome"
                new_obx_index +=1 
                temp_obx.append(obxSegment) 
            
        h_t = self.h[:]
        l = len(self.h)
        for i in range(l):
            del h_t[l-i-1]
        for i in range(len(self.h)):
            if(self.h[i][0][0]!="OBX"):
                h_t.append(self.h[i])
        h_t.extend(temp_obx)
        return h_t

    def get_first_obx_index(self):
        idx = 0
        for seg in self.h:
            if seg[0][0] == 'OBX':
                return idx
            idx += 1
        return -1
    # Assuming insertion is just above first OBX segment
    def update_comments(self, comments):
        comments_arr = comments.split("\n")
        obx_idx = self.get_first_obx_index()
        if (obx_idx == -1):
            print("OBX segment not found, so appending it")
            obx_idx = len(self.h) - 1 
        i=1
        for comment in comments_arr:
            self.h.append('NTE|{}|L|{}'.format(i,comment))
            obx_idx += 1
            i += 1


def create_connection():
    conn = None
    try:
        conn = mysql.connector.connect(
            host=MYSQL_HOST,
            user=MYSQL_USERNAME,
            passwd=MYSQL_PASSWORD,
            database=MYSQL_DATABASE
        )
    except:
        print(traceback.format_exc())
    return conn
SQL_CONNECTION = create_connection()


def updateRowInDatabase(sample, PLMO, MRN, ptName, ptSex, ptAge, ordDept, excelFileName):
    cur = SQL_CONNECTION.cursor()

    updateSql = "UPDATE "+ COVID_19_EPIC_UPLOAD_TABLE +" set PLMO_Number= %s, MRN = %s, PATIENT_NAME = %s, PATIENT_SEX = %s, PATIENT_AGE = %s, ORDERING_DEPARTMENT = %s, EPIC_UPLOAD_TIMESTAMP = %s WHERE CONTAINER_ID = %s and SOURCE_EXCEL_FILE = %s;"
    #print(sample.completeSampleName, type(sample.name), "!!!!!!",(sample.name, PLMO, get_current_formatted_date(), sample.nCoV_N1, sample.nCoV_N2, sample.nCoV_N3, sample.RP, sample.result))
    cur.execute(updateSql, (PLMO[0], MRN, ptName, ptSex, ptAge, ordDept, str(datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S")), sample.name, excelFileName))
    SQL_CONNECTION.commit()
    cur.close()


# arguments are <class Sample> sample, str PLMO
def addRowInDatabase(sample, PLMO, MRN, ptName, ptSex, ptAge, ordDept, excelFileName):
    cur = SQL_CONNECTION.cursor()

    findsql = "SELECT * from " + COVID_19_EPIC_UPLOAD_TABLE + " where CONTAINER_ID = %s and SOURCE_EXCEL_FILE = %s;"
    cur.execute(findsql, (sample.name, excelFileName))
    rows = cur.fetchall()
    # print(rows)
    # print("-----------------")
    if len(rows) > 0:
        updateSql = "UPDATE " + COVID_19_EPIC_UPLOAD_TABLE + " set QUANTSTUDIO_SPECIMEN_ID = %s, 2019nCoV_N1 = %s, 2019nCoV_N2 = %s, 2019nCoV_N3 = %s, RP = %s, RESULT = %s where CONTAINER_ID = %s and SOURCE_EXCEL_FILE = %s;" 
        cur.execute(updateSql, (sample.completeSampleName, sample.nCoV_N1, sample.nCoV_N2, sample.nCoV_N3, sample.RP, sample.result, sample.name, excelFileName))
        SQL_CONNECTION.commit()
    else:
        insertSql = "INSERT INTO "+ COVID_19_EPIC_UPLOAD_TABLE +" VALUES(%s, %s,%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);"
        cur.execute(insertSql, (sample.completeSampleName, sample.name, PLMO, MRN, ptName, ptSex, ptAge, ordDept, excelFileName, "", sample.nCoV_N1, sample.nCoV_N2, sample.nCoV_N3, sample.RP, sample.result, "", ""))
        SQL_CONNECTION.commit()
    cur.close()

# method to write the entire database table to excel
def writeDataToExcel(excelName):
    xldf = pd.read_sql_query('select * from '+ COVID_19_EPIC_UPLOAD_TABLE +' where SOURCE_EXCEL_FILE = "'+ excelName +'" ;', SQL_CONNECTION)
    xldf = xldf.drop("SOURCE_EXCEL_FILE", 1)
    xldf = xldf.drop("ORDERING_DEPARTMENT", 1)
    xldf = xldf.drop("2019nCoV_N3", 1)
    xldf = xldf.drop("RLU_SCORE", 1)    
    xldf = xldf.drop("RLU_FLAG", 1)

    cols = xldf.columns.tolist()
    upload_col = cols.pop(cols.index("EPIC_UPLOAD_TIMESTAMP"))
    cols.insert(len(cols), upload_col)
    xldf = xldf.reindex(columns= cols)
    RESULT_LOG = COVID_19_TEST_SAMPLE_LOG + "/" + \
        excelName.split("/")[-1].replace("_SAMPLE_RESULTS_UPDATED_ID","_SAMPLE_EPIC_UPLOAD_LOG")
    SAMPLE_MAP_FILE = excelName.replace("_SAMPLE_RESULTS_UPDATED_ID", "_SAMPLE_MAP")

    sampleMapDf = pd.read_excel(SAMPLE_MAP_FILE)
    sampleMapDf["UPLOADED_TO_EPIC"] = "No"
    for index, row in sampleMapDf.iterrows():
        if xldf['QUANTSTUDIO_SPECIMEN_ID'].str.contains(str(row["Internal_Sample_ID"])).any() :
            sampleMapDf.at[index, "UPLOADED_TO_EPIC"] = "Yes"
    
    try:
        # xldf.to_excel(RESULT_LOG , index=False)
        with pd.ExcelWriter(RESULT_LOG) as writer:
            xldf.to_excel(writer, index=False, sheet_name='Sheet 1')
            #sampleMapDf.to_excel(writer, index=False, sheet_name='Sheet 2')
        print("done writeToExcel method and writing done to -->", RESULT_LOG )
    except:
        print("unable to save status excel, please close it")

def checkIncomingHl7(sampleDict, commentsDict, excelFile):
    UPLOAD_PATH = MIRTH_GATORSEQ + '/RESULTS'
    ORDERS_ARCHIVE_DIR = MIRTH_GATORSEQ + '/ORDERS_ARCHIVE/'
    ORDERS_DIR = MIRTH_GATORSEQ + '/ORDERS/'
    allhl7filenames = []
    for (dirpath, dirnames, filenames) in os.walk(ORDERS_DIR):
        allhl7filenames.extend(filenames)
        break
    for hl7_file_name in allhl7filenames:
        try:
            hl7file = open(ORDERS_DIR + hl7_file_name, mode="r").read()
        except:
            continue
        arr = hl7file.split("\n\n") #split by blank lines
        #Iterate all HL7 messages in file
        c=0
        for hl7msg in arr: 
            if hl7msg: #check if message not empty
                msg_unix_fmt = hl7msg.replace("\n","\r")
                h = hl7.parse(msg_unix_fmt)
                newHl7 = hl7update(h)
                #Read message id from HL7 message
                try:
                    messageId = str(h['OBR'][0][3]).replace("^", "")
                except:
                    continue
                if (not messageId):
                    continue
                plm = None
                try:
                    plm = h['ORC'][0][2]
                except:
                    print("---------could not find PLMO in hl7 file: ", hl7_file_name)
                    continue
                
                mrn = ""
                try:
                    mrn = h['PID'][0][3][0][0]
                except:
                    print("---------could not find MRN in hl7 file: ", hl7_file_name)
        
                ptName = ""
                try:
                    ptName = h['PID'][0][5][0]
                except:
                    print("---------could not find PATIENT_NAME in hl7 file: ", hl7_file_name) 

                ptSex = ""
                try:
                    ptSex = h['PID'][0][8][0]
                except:
                    print("---------could not find PATIENT_SEX in hl7 file: ", hl7_file_name) 

                ptAge = -1
                try:
                    ptAge = 2020 - int(h['PID'][0][7][0][:4]) 
                except:
                    print("---------could not find PATIENT_AGE in hl7 file: ", hl7_file_name) 

                ordDept = ""
                try:
                    ordDept = h['OBR'][0][15][0]
                except:
                    print("---------could not find Ordering_DEPT in hl7 file: ", hl7_file_name) 


                # search for messageId in the sampleDict
                #if messageId == "100047187": #100047166  100047187
                if messageId in sampleDict or plm[0] in sampleDict:
                   # print("--------found----------")
                    methodology_code = ""
                    if sampleDict.get(messageId) is not None:
                        givenSample =  sampleDict.get(messageId)
                        methodology_code = commentsDict[messageId]
                    else:
                        givenSample = sampleDict.get(plm[0])
                        methodology_code = commentsDict[plm[0]]
                    #print("processing hl7 input file: ", hl7_file_name)                   
                    newHl7.update_msh_segment()
                    newHl7.update_orc_segment()
                    newHl7.update_obr_segment()
                    newHl7.update_obx_segment()
                    if methodology_code:
                        #print("comments for {} = {}".format(givenSample.name, methodology_code.upper()))
                        newHl7.update_comments(open(comment_files[methodology_code.upper()], mode="r",  encoding='utf-8').read())
                    else:
                        print("Methodology is empty to assign comments for: ",givenSample.name)
                    h = newHl7.update_obx_seg_containing_gene( givenSample )
                    
                    out_file_path = UPLOAD_PATH + '/hl7-COVID_19-{}-output.txt'.format(messageId)
                    if h:
                        with open(out_file_path, 'w' ,  encoding='utf-8') as f:
                            f.write(str(h))
                        print("Out file available at :",out_file_path)
                        move(ORDERS_DIR + hl7_file_name, ORDERS_ARCHIVE_DIR + 'COVID_19_processed_' + get_current_formatted_date() + "-" + hl7_file_name) 
                        if plm:
                            updateRowInDatabase(givenSample, plm, str(mrn), str(ptName), str(ptSex), str(ptAge), str(ordDept), excelFile )
                    

class Sample:
    def __init__(self, sample_name, completeSampleName):
        self.name = str(sample_name)
        self.nCoV_N1 = None
        self.nCoV_N2 = None
        self.nCoV_N3 = None
        self.RP = None
        self.result = None
        self.completeSampleName = completeSampleName

    def __str__(self):
        return str( self.completeSampleName) + ": " + str(self.nCoV_N1) + " & " + str(self.nCoV_N2) + " & "  + str(self.nCoV_N3) + " & " + str(self.RP) + " & " + str(self.result)

    def __repr__(self):
        return str( self.completeSampleName) + ": " + str(self.nCoV_N1) + " & " + str(self.nCoV_N2) + " & "  + str(self.nCoV_N3) + " & " + str(self.RP) + " & " + str(self.result) + " | "

def isFloatValue(value, maxThreshold):
    try:
        val = float(value)
        #ToDo: check RP condition, if RP > 35 then is it considered to be Undetermined
        if maxThreshold and val > maxThreshold:
            return False
        return True

    except:
        if value == "Undetermined":
            return False
        else:
            print("-----------ERROR: unable to identify the value-------------->", value)
            return False
    
def addSampleDictToDatabase(sampleDict, excelName):
    for key in sampleDict.keys():
        sample = sampleDict[key]
        addRowInDatabase(sample, "", "", "", "", "", "", excelName)
    
    print("-> all samples added to database <-")

def heatMapToExcel(mask, writer, sheetname, run, filename):
    pivot = pd.pivot_table(mask, 'Val', index = ['Row', 'ID'], columns=['Column'], aggfunc='first')
    # Export to excel sheet with two worksheets for each run
    pivot.to_excel(writer, sheet_name=sheetname)
    # Creating pivot table with extracted letters as rows and columns
    tables = pd.pivot_table(run, 'CT',  index= ['Row'], columns=['Column'], aggfunc='first')
    # Replace 'Undetermined' values with zeroes and round the values to 1 digit.
    tables = tables.replace('Undetermined', 0)
    tables = tables.round(1)
    hmap = sns.heatmap(tables, annot=True, fmt='g', square=True, annot_kws={"size": 6}, cbar_kws={"orientation": "horizontal"})
    bottom, top = hmap.get_ylim()
    x = hmap.set_ylim(bottom + 0.5, top - 0.5)
    hmap.set_ylabel('')    
    hmap.set_xlabel('')
    plt.title(sheetname)
    fig = plt.figure()
    fig = hmap.get_figure()
    fig.savefig(filename + "_QC_" + sheetname +'.png')
    worksheet = writer.sheets['{}'.format(sheetname)]
    worksheet.insert_image('D20', filename + "_QC_" + sheetname +'.png')   
    return (pivot, tables) 

def createHeatMapDiagram(output, filename):
    writer = pd.ExcelWriter(filename + "_QC.xlsx", engine='xlsxwriter')
    runs = {}
    tables = {}
    mask = {}
    d1 = {}
    d2 = {}
    pivot = {}
    runs_list = output['Plate ID'].unique()

    for i in range(len(runs_list)):
        if not pd.isna(runs_list[i]):
            # Extracting dataframes for two runs
            runs[i] = output.loc[output['Plate ID']==runs_list[i]]
            # Masking data to duplicate rows for better view in the Output file
            mask[i] = runs[i]
            d1[i] = runs[i].assign(Val = runs[i]['CT'],ID = 'CT' ) 
            d2[i] = runs[i].assign(Val = runs[i]['CONTAINER_ID'], ID = 'CONTAINER_ID')
            mask[i] = pd.concat([d1[i], d2[i]]).sort_index().reset_index(drop = True)
            (pivot[i], tables[i]) = heatMapToExcel(mask[i], writer, runs_list[i], runs[i], filename)

    if len(mask) > 0:
        allsheets = pd.DataFrame(columns = mask[0].columns)
        all_maps = {"A":["A","B"],"B":["C","D"],"C":["E","F"],"D":["G","H"],"E":["I","J"],"F":["K","L"],"G":["M","N"],"H":["O","P"]}
        for c in range(1,13):
            for r in all_maps.keys():
                actualrow = mask[0].loc[(mask[0]["Row"] == r) & (mask[0]["Column"] == c)]
                if not actualrow.empty:
                    r1 = actualrow.copy(deep=False)
                    r1["Row"] = all_maps[r][0]
                    r1["Column"] = c
                    allsheets = allsheets.append(r1)
                if len(mask) > 2:
                    actualrow = mask[2].loc[(mask[2]["Row"] == r) & (mask[2]["Column"] == c)]
                    if not actualrow.empty:
                        r1 = actualrow.copy(deep=False)
                        r1["Row"] = all_maps[r][1]
                        r1["Column"] = c
                        allsheets = allsheets.append(r1)
                if len(mask) > 1:
                    actualrow = mask[1].loc[(mask[1]["Row"] == r) & (mask[1]["Column"] == c)]
                    if not actualrow.empty:
                        r1 = actualrow.copy(deep=False)
                        r1["Row"] = all_maps[r][0]
                        r1["Column"] = c+12
                        allsheets = allsheets.append(r1)
                if len(mask) > 3:        
                    actualrow = mask[3].loc[(mask[3]["Row"] == r) & (mask[3]["Column"] == c)]
                    if not actualrow.empty:
                        r1 = actualrow.copy(deep=False)
                        r1["Row"] = all_maps[r][1]
                        r1["Column"] = c+12
                        allsheets = allsheets.append(r1)

        heatMapToExcel(allsheets, writer, "FINAL", allsheets, filename)

    writer.save()    

def createHeatMapTable(df1, df2, filename):
    # Rename 'Internal_Sample-ID' to 'Sample Name' and merge two sheets based on 'Sample Name'
    df1 = df1.rename(columns = {'Internal_Sample_ID' : 'Sample Name'})
    output = pd.merge(df1,df2, on="Sample Name")

    # Cleaning data by considering only 2019nCoV_N1 samples
    output.drop(output.loc[output['Target Name']!='2019nCoV_N1'].index, inplace=True)
    
    # Creating two new columns by splitting Sample Name and deleting the old column
    output['Well_Name'] = output['Well_Name'].str[0:1] + ' ' + output['Well_Name'].str[1:]
    new = output["Well_Name"].str.split(" ", n = 1, expand = True) 
    output["Row"]= new[0] 
    output["Column"]= new[1] 
    output.drop(columns =["Well_Name"], inplace = True) 

    # Convert Column to numeric for sorted columns in the output
    output["Column"] = pd.to_numeric(output["Column"], errors='coerce')

    # Creating columns for better representation in the output
    output["ID"] = ' '
    output["Val"] = ' '

    createHeatMapDiagram(output, filename)    

if __name__ == "__main__":
    os.chdir(COVID_19_TEST_INPUT_FOLDER)
    _files = filter(os.path.isfile, os.listdir(COVID_19_TEST_INPUT_FOLDER))
    excel_files = [os.path.join(COVID_19_TEST_INPUT_FOLDER, f) for f in _files if "$" not in f] # add path to each file

    fileNames = []
    toProcess = []
    toUploadSamples = {}
    toUpdateComments = {}
    for f in excel_files:
        if "_SAMPLE_MAP" in f:
            sampleGroupName = f[:f.index("_SAMPLE_MAP")]
            if sampleGroupName + "_SAMPLE_RESULTS_UPDATED_ID.xlsx" not in excel_files:
                fileNames.append(sampleGroupName)
   
    for f in fileNames:
        sampleResult = f + "_SAMPLE_RESULTS.xlsx"
        sampleMap = f + "_SAMPLE_MAP.xlsx"
        sampleToContainer = {}
        hscBlankCheck = {}
        samplesUpload = {}
        updateComments = {}
        df = pd.read_excel(sampleMap)
        for i, row in df.iterrows():
            if "hsc" in str(row["Container_ID"]).lower() or "blank" in str(row["Container_ID"]).lower():
                hscBlankCheck[row["Internal_Sample_ID"]] = str(row["Container_ID"])

            sampleToContainer[row["Internal_Sample_ID"]] = row["Container_ID"]
            samplesUpload[str(row["Container_ID"])] = str(row["Upload"])
            updateComments[str(row["Container_ID"])] = str(row["Methodology"])
        
        results_df = pd.read_excel(sampleResult, skiprows=range(0,41))
        results_df["CONTAINER_ID"] = None
        
        error = False
        for index, row in results_df.iterrows():
            if sampleToContainer.get(row["Sample Name"]):
                results_df.at[index, "CONTAINER_ID"] = sampleToContainer[row["Sample Name"]]
                if row["Sample Name"] in hscBlankCheck.keys():
                    if ("hsc" in hscBlankCheck[row["Sample Name"]].lower() and ((row["Target Name"] == "RP" and not isFloatValue(row["CT"], 35)) or (row["Target Name"] != "RP" and not isFloatValue(row["CT"], 40)))):
                        error = True
                    elif "blank" in hscBlankCheck[row["Sample Name"]].lower() and row["CT"] != "Undetermined":
                        error = True   
            else:
                results_df.at[index, "CONTAINER_ID"] = str(row["Sample Name"]) + " <No containerId>"
        
        if error:
            error_data = {'Error' : ['There is probably an error in the input files. Check the HSC and BLANK samples.']}
            error_df = pd.DataFrame(error_data, columns = ['Error'])
            error_df.to_excel(f + "_ERROR.xlsx", index = False)
            continue

        toProcess.append(f + "_SAMPLE_RESULTS_UPDATED_ID.xlsx")
        results_df.to_excel(f + "_SAMPLE_RESULTS_UPDATED_ID.xlsx", index=False)
        #####logic to generate heatmap and table#####
        createHeatMapTable(df, results_df, f)
        toUploadSamples[f + "_SAMPLE_RESULTS_UPDATED_ID.xlsx"] = samplesUpload
        toUpdateComments[f + "_SAMPLE_RESULTS_UPDATED_ID.xlsx"] = updateComments

    for f in toProcess:
        sampleDict = {}
        eachExcel = os.path.join(COVID_19_TEST_INPUT_FOLDER, f)
        xldf = pd.read_excel(eachExcel)
        samplesUpload = toUploadSamples[f]
        updateComments = toUpdateComments[f]
        #xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        # generate sample dictionary in below format:
        # {<sample name>: <Class Sample>}
        #print(xldf.head())  
        for index, row in xldf.iterrows():
            sampleName = str(row["CONTAINER_ID"])
            if sampleName in samplesUpload.keys() and samplesUpload[sampleName].lower() == "yes":
                if 'PLMO' in sampleName:
                    plm = re.findall(r"PLMO\d+-\d+" , sampleName)
                    if len(plm):
                        sampleName = plm[0]

                sampleName = sampleName.replace("\\", "")            

                targetName = row["Target Name"]
                value = row["CT"]
                
                #ToDo: sampleName = <PLMO of a sample> or <id number of a sample>
                if sampleDict.get(sampleName) is None:
                    sampleDict[sampleName] = Sample(sampleName, str(row["Sample Name"]))
                if targetName == "RP":# and not math.isnan(value):
                    setattr(sampleDict[sampleName], "%s" % targetName, value)
                else:# not math.isnan(value):
                    #ToDO: is there any better of deriving 'nCoV_N1' from '2019nCoV_N1'? (python does not allow variable names to start with number)
                    setattr(sampleDict[sampleName], "%s" % targetName[4:], value)

        for sampleName in sampleDict.keys():
            sample = sampleDict[sampleName]
            '''
            #if both n1 and n2 are < 40 -> detected
            #elif if any of them is (>40 or undertermined) AND not all of them are (>40 or undertermined) -> Inderterminate
            #if RP is (>40 or undertermined) OR both n1 and n2 are float and > 40 -> invalid
            #if both n1 and n2 are "undetermined" -> not detected
            if all([isFloatValue(sample.nCoV_N1, 40), isFloatValue(sample.nCoV_N2, 40)]):#, isFloatValue(sample.nCoV_N3, None)]):
                sample.result = "Detected"
                continue
            elif any([not isFloatValue(sample.nCoV_N1, 40), not isFloatValue(sample.nCoV_N2, 40)]) and not all([not isFloatValue(sample.nCoV_N1, 40), not isFloatValue(sample.nCoV_N2, 40)]):
                sample.result = "Indeterminate"
                continue
            
            if (sample.RP and not isFloatValue(sample.RP, 40.0)) or ( type(sample.nCoV_N1) == "float" and sample.nCoV_N1 > 40 and type(sample.nCoV_N2) == "float" and sample.nCoV_N2 > 40  ):
                #INVALID Case is to be handled by pathologists separately and hence just continuing
                sample.result = "Invalid"
                continue
            if all([not isFloatValue(sample.nCoV_N1, None), not isFloatValue(sample.nCoV_N2, None)]):
                sample.result = "Not Detected"
            '''            
            if all([isFloatValue(sample.nCoV_N1, 40)]): #, isFloatValue(sample.nCoV_N3, None)]):
                sample.result = "Detected"
                continue
            elif any([not isFloatValue(sample.nCoV_N1, 40)]) and not all([not isFloatValue(sample.nCoV_N1, 40)]):
                sample.result = "Indeterminate"
                continue
            
            if (sample.RP and not isFloatValue(sample.RP, 40.0)) or ( type(sample.nCoV_N1) == "float" and sample.nCoV_N1 > 40):
                #INVALID Case is to be handled by pathologists separately and hence just continuing
                sample.result = "Invalid"
                continue
            if all([not isFloatValue(sample.nCoV_N1, None)]):
                sample.result = "Not Detected"           
            
            if sample.result is None:
                print("------ Last resort reached -----", sample)
                #del sampleDict[sampleName]
                sample.result = "Invalid"
        #print("below is the dictionary of all samples:")
        #print(sampleDict["PLMO20-000129"])
        #print(updateComments)
        if sampleDict:
            addSampleDictToDatabase(sampleDict, f)
            checkIncomingHl7(sampleDict, updateComments, f)
            writeDataToExcel(f)  
            #addSampleDictToExcel(sampleDict, f, True) 
        
    #writeDataToExcel("/ext/path/DRL/Molecular/COVID19/COVID_19_QuantStudio/ProdEnv/Results/2020-03-20 203810_QuantStudio_export_UPDATED_CONTAINER_ID.xlsx")
    SQL_CONNECTION.close()
