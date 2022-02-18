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
print("Run start time: ", str(datetime.datetime.now()) + "\n")

script_path = os.path.dirname(os.path.abspath( __file__ ))
parent_path = os.path.abspath(os.path.join(script_path, '..'))
comments_file_path = parent_path + "/COVID_19/POOL_Methodology_Comment_SingleLine.txt"
CONFIG_FILE = parent_path +"/linux_gatorseq.config.yaml"
RLU_NOMINAL_SCORE = 325
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
COVID_19_TEST_INPUT_FOLDER = replace_env(config_dict['POOL_COVID_19_TEST_INPUT_FOLDER']) #'G:\DRL\Molecular\Assays\PGX\PGX_Beaker_Interface' 
CONFIG_TOKENS_FILE = parent_path  + "/" + config_dict['CONFIG_TOKENS_FILE']
MIRTH_GATORSEQ = config_dict['MIRTH_GATORSEQ']
if CODE_ENV=='ProdEnv':
    MIRTH_GATORSEQ += '/PROD'
    COVID_19_EPIC_UPLOAD_TABLE = replace_env(config_dict['COVID_19_EPIC_UPLOAD_TABLE_PROD'])
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
   
    def update_obx_seg_containing_gene(self, result, runid):
        updates = 0
        temp_obx = self.h[:]
        l = len(self.h)
        for i in range(l):
            del temp_obx[l-i-1]
        new_obx_index = 1
        for obxSegment in self.h['OBX']:
            if obxSegment[3][0][1][0] == "SARS-COV-2, NAA":
                obxSegment[5][0] = result
                obxSegment[1] = new_obx_index
                machineId = ""
                if runid.startswith('003109'):
                    machineId = "UFHPL Panther1"
                elif runid.startswith('003253'):
                    machineId = "UFHPL Panther2"
                if len(machineId) > 0:
                    obxSegment[18][0] = machineId
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


def updateRowInDatabase(containerId, PLMO, MRN, ptName, ptSex, ptAge, ordDept, excelFileName):
    cur = SQL_CONNECTION.cursor()

    updateSql = "UPDATE "+ COVID_19_EPIC_UPLOAD_TABLE +" set PLMO_Number= %s, MRN = %s, PATIENT_NAME = %s, PATIENT_SEX = %s, PATIENT_AGE = %s, ORDERING_DEPARTMENT = %s, EPIC_UPLOAD_TIMESTAMP = %s WHERE CONTAINER_ID = %s and SOURCE_EXCEL_FILE = %s;"
    #print(sample.completeSampleName, type(sample.name), "!!!!!!",(sample.name, PLMO, get_current_formatted_date(), sample.nCoV_N1, sample.nCoV_N2, sample.nCoV_N3, sample.RP, sample.result))
    cur.execute(updateSql, (PLMO[0], MRN, ptName, ptSex, ptAge, ordDept, str(datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S")), containerId, excelFileName))
    SQL_CONNECTION.commit()
    cur.close()


# arguments are <class Sample> sample, str PLMO
def addRowInDatabase(containerId, result, RLUScore, RLUFlag, PLMO, MRN, ptName, ptSex, ptAge, ordDept, excelFileName):
    cur = SQL_CONNECTION.cursor()

    findsql = "SELECT * from " + COVID_19_EPIC_UPLOAD_TABLE + " where CONTAINER_ID = %s and SOURCE_EXCEL_FILE = %s;"
    cur.execute(findsql, (containerId, excelFileName))
    rows = cur.fetchall()
    # print(rows)
    # print("-----------------")
    if len(rows) > 0:
        updateSql = "UPDATE " + COVID_19_EPIC_UPLOAD_TABLE + " set QUANTSTUDIO_SPECIMEN_ID = %s, RESULT = %s, RLU_SCORE = %s, RLU_FLAG = %s where CONTAINER_ID = %s and SOURCE_EXCEL_FILE = %s;" 
        cur.execute(updateSql, ("pooled", result, RLUScore, RLUFlag, containerId, excelFileName))
        SQL_CONNECTION.commit()
    else:
        insertSql = "INSERT INTO "+ COVID_19_EPIC_UPLOAD_TABLE +" VALUES(%s, %s,%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);"
        cur.execute(insertSql, ("pooled", containerId, PLMO, MRN, ptName, ptSex, ptAge, ordDept, excelFileName, "", "", "", "", "", result, RLUScore, RLUFlag))
        SQL_CONNECTION.commit()
    cur.close()

# method to write the entire database table to excel
def writeDataToExcel(excelName, sampleToResult, sampleToPool):
    xldf = pd.read_sql_query('select CONTAINER_ID, EPIC_UPLOAD_TIMESTAMP, PLMO_NUMBER, PATIENT_NAME, PATIENT_SEX, PATIENT_AGE, RLU_SCORE, RLU_FLAG from '+ COVID_19_EPIC_UPLOAD_TABLE +' where SOURCE_EXCEL_FILE = "'+ excelName +'" ;', SQL_CONNECTION)
    dbdict = {}
    for i, row in xldf.iterrows():
        dbdict[row['CONTAINER_ID']] = [row['EPIC_UPLOAD_TIMESTAMP'], row['PLMO_NUMBER'], row['PATIENT_NAME'], row['PATIENT_SEX'], row['PATIENT_AGE'], row['RLU_SCORE'], row['RLU_FLAG']]
    
    #dbdict = xldf.set_index('CONTAINER_ID')['EPIC_UPLOAD_TIMESTAMP'].to_dict()
    writeToList = []
    #print("dbdict === {}".format(dbdict))
    for key, value in sampleToResult.items():    
        tempdict = {}
        tempdict['CONTAINER_ID'] = key
        tempdict['POOL_ID'] = sampleToPool[key]      
        tempdict['VALID_FLAG'] = sampleToResult[key][1]
        tempdict['RESULT'] = sampleToResult[key][2]
        tempdict['RLU_SCORE'] = sampleToResult[key][3]
        tempdict['RLU_FLAG'] = sampleToResult[key][4]         
        if key not in dbdict.keys() or dbdict[key][0] is None or dbdict[key][0] == "":
            tempdict['EPIC_UPLOADED'] = "NO"
            tempdict['TIMESTAMP'] = ""
        else:
            tempdict['EPIC_UPLOADED'] = "YES"
            tempdict['TIMESTAMP'] = dbdict[key][0]
            tempdict['PLMO_NUMBER'] = dbdict[key][1]
            tempdict['PATIENT_NAME'] = dbdict[key][2]
            tempdict['PATIENT_SEX'] = dbdict[key][3]
            tempdict['PATIENT_AGE'] = dbdict[key][4]
       
        writeToList.append(tempdict)

    xldf = pd.DataFrame(writeToList)
    columnsTitles = ['VALID_FLAG','POOL_ID','RLU_FLAG','RLU_SCORE','RESULT','CONTAINER_ID','EPIC_UPLOADED','TIMESTAMP','PLMO_NUMBER','PATIENT_NAME','PATIENT_SEX','PATIENT_AGE']
    xldf = xldf.reindex(columns=columnsTitles)
    xldf.sort_values(['RESULT','POOL_ID'], ascending=[False,True], inplace=True)
    RESULT_LOG = excelName + "_FINAL.xlsx"
    
    try:
        with pd.ExcelWriter(RESULT_LOG) as writer:
            xldf.to_excel(writer, index=False, sheet_name='Sheet 1')
        print("done writeToExcel method and writing done to -->", RESULT_LOG )
    except:
        print("unable to save status excel, please close it")

def checkIncomingHl7(sampleDict, sampleResultDict, excelFile):
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
        #c=0
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
                    if sampleDict.get(messageId) is not None:
                        givenSampleResult =  sampleDict.get(messageId)
                    else:
                        givenSampleResult = sampleDict.get(plm[0])                  
                    newHl7.update_msh_segment()
                    newHl7.update_orc_segment()
                    newHl7.update_obr_segment()
                    newHl7.update_obx_segment()
                    newHl7.update_comments(open(comments_file_path, mode="r",  encoding='utf-8').read())
                    h = newHl7.update_obx_seg_containing_gene( givenSampleResult, sampleResultDict[messageId][0] )
                    
                    out_file_path = UPLOAD_PATH + '/hl7-pooled-COVID_19-{}-output.txt'.format(messageId)
                    if h:
                        with open(out_file_path, 'w' ,  encoding='utf-8') as f:
                            f.write(str(h))
                        print("Out file available at :",out_file_path)
                        move(ORDERS_DIR + hl7_file_name, ORDERS_ARCHIVE_DIR + 'POOLED_COVID_19_processed_' + get_current_formatted_date() + "-" + hl7_file_name) 
                        if plm:
                            updateRowInDatabase(messageId, plm, str(mrn), str(ptName), str(ptSex), str(ptAge), str(ordDept), excelFile )
                    

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
    
def addSampleDictToDatabase(sampleResult, resultData, excelName):
    for containerId in sampleResult.keys():
        result = sampleResult[containerId]
        
        addRowInDatabase(containerId, result, resultData[containerId][3], resultData[containerId][4], "", "", "", "", "", "", excelName)
    
    print("-> all samples added to database <-")


if __name__ == "__main__":
    os.chdir(COVID_19_TEST_INPUT_FOLDER)
    _files = filter(os.path.isfile, os.listdir(COVID_19_TEST_INPUT_FOLDER))
    excel_files = [os.path.join(COVID_19_TEST_INPUT_FOLDER, f) for f in _files if "$" not in f] # add path to each file

    fileNames = []
    for f in excel_files:
        if "_SAMPLE_POOL_MAP" in f:
            sampleGroupName = f[:f.index("_SAMPLE_POOL_MAP")]
            if sampleGroupName + "_FINAL.xlsx" not in excel_files:
                fileNames.append(sampleGroupName)
    
    for f in fileNames:
        sampleResult = f + "_SAMPLE_POOL_RESULTS.lis"
        sampleMap = f + "_SAMPLE_POOL_MAP.xlsx"
        sampleToPool = {}
        df = pd.read_excel(sampleMap)
        for i, row in df.iterrows():
            if not pd.isna(row["Source Sample Barcode"]):
                sampleToPool[str(row["Source Sample Barcode"]).split(".")[0]] = str(row["Pooled Sample Barcode"]).split(".")[0]
        
        #print("sample to pool = {}".format(sampleToPool))
        results_df = pd.read_csv(sampleResult, sep='\t')
        pool_results = {}
        #print("results_df = {}".format(results_df))
        for i, row in results_df.iterrows():
            if row["Specimen Barcode"] in sampleToPool.values():
                if int(row["Interpretation 1"]) < RLU_NOMINAL_SCORE:
                    flag = "NORMAL"
                else:
                    flag = "HIGH"
                pool_results[row["Specimen Barcode"]] = [row["Run ID"], row["Interpretation 2"], row["Interpretation 3"], row["Interpretation 1"], flag]
        
        #print("poolresults = {}".format(pool_results))
        #contains results corresponding to each containerId
        sampleToResult = {}
        for i, row in df.iterrows():
            if str(row["Pooled Sample Barcode"]).split(".")[0] in pool_results.keys():
                sampleToResult[str(row["Source Sample Barcode"]).split(".")[0]] = pool_results[str(row["Pooled Sample Barcode"]).split(".")[0]]

        containerToResult = {}
        for containerId in sampleToResult:
            if sampleToResult[containerId][1].lower() == "valid" and sampleToResult[containerId][2].lower() == "negative" and int(sampleToResult[containerId][3]) < RLU_NOMINAL_SCORE:
                containerToResult[containerId] = "Not Detected"
        
        #print("sampleToResult=== {}".format(sampleToResult))        
        #print("container to result === {}".format(containerToResult))
        #add all samples to database
        addSampleDictToDatabase(containerToResult, sampleToResult, f)

        #logic to add the corresponding hl7 file to RESULTS folder
        checkIncomingHl7(containerToResult, sampleToResult, f)

        print(f)
        writeDataToExcel(f, sampleToResult, sampleToPool)

    SQL_CONNECTION.close()
