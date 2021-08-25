#!/usr/bin/env python

import numpy as np
import csv
import os
def LoadChimeCands(CandidateCSV='chime.csv'):
    """Summary or Description of the Function

    Parameters:
    ra (str): ra position (hh:mm:ss)
    dec(str): dec position (dd:mm:ss)
    dm (str): dm (m/cm^3)
    CandidateCSV (str): filepath to CandidateCSV file

    Returns:
    results(str): String containing the if the pulsar was found or not
    """
    chime_cand = {}
    with open (CandidateCSV,'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        for row in csv_reader:
            chime_ra = float(row[1])
            chime_dec = float(row[2])
            chime_dm =float(row[3])
            chime_cand.update({row[0]:[chime_ra,chime_dec,chime_dm]})
    return chime_cand

def load_header_localised(path):
    if not os.path.exists(path):
        return None,None,None,None
    files = os.listdir(path)
    if len(files)==0:
        return None,None,None,None
    else:
        ra_array=[]
        ra_errors=[]
        dec_array=[]
        dec_errors=[]
        for file in files:
            if '.pkl' in file:

                try:
                    my_localisations=np.load(os.path.join(path,file),allow_pickle=1)['results'][99]
                except:
                    continue
                ra_array.append(my_localisations[0])
                ra_errors.append(my_localisations[1])
                dec_array.append(my_localisations[2])
                dec_errors.append(my_localisations[3])
        ra_spread=max(ra_array)-min(ra_array)
        dec_spread = max(dec_array)-min(dec_array)
        ra_error = max(ra_errors[np.squeeze(np.argwhere(ra_array==min(ra_array)))],ra_errors[np.squeeze(np.argwhere(ra_array==max(ra_array)))])
        dec_error = max(dec_errors[np.squeeze(np.argwhere(dec_array==min(dec_array)))],dec_errors[np.squeeze(np.argwhere(dec_array==max(dec_array)))])

        if ((ra_spread+ra_error)>5) | ((dec_spread+dec_error)>5):
            print('look into this folder, spread too big! '+path+' number of events '+str(len(ra_array))+' ra spread '+str(ra_spread)+' err '+str(ra_error)+' dec spread '+str(dec_spread))
            return None,None,None,None
        else:
            ra = np.mean(ra_array)
            dec= np.mean(dec_array)
            ra_error =ra_error+ra_spread
            dec_error = dec_error+dec_spread
            return ra,ra_error,dec,dec_error

def load_new_sources(filename,header_localised_folder=None):
    '''This function will plot all the new confirmed sources'''
    if 'csv' in filename:
        with open(filename) as csv_file:
            my_new_sources={}
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                cluster_id = row[0]
                # if header_localised_folder:
                #     path=os.path.join(header_localised_folder,'cluster'+cluster_id+'/')
                #     ra,ra_tol,dec,dec_tol = load_header_localised(path)
                #     if ra:
                #         ra_hhmmss,dec_ddmmss = convert_to_hoursminsec(float(ra),float(dec))
                #         #do a check for the clustering ra and header-localise ra
                #         cluster_ra,cluster_dec = convert_to_deg(row[1],row[2])
                #         if (abs(cluster_ra-ra)<5) & (abs(cluster_dec-dec)<5):
                #             new_item = {row[0]:np.array([ra_hhmmss,ra_tol,dec_ddmmss,dec_tol,row[3]])}
                #         else:
                #             print('significant differences in '+'cluster'+cluster_id)

                #     else:
                #         new_item = {row[0]:np.array([row[1],ra_dec_tol,row[2],ra_dec_tol,row[3],row[-1]])}
                # else:
                #     new_item = {row[0]:np.array([row[1],ra_dec_tol,row[2],ra_dec_tol,row[3],row[-1]])}
                new_item = {row[0]:np.array([row[4],row[5],row[3]])}
                my_new_sources.update(new_item)
        return my_new_sources
    else:
        print('please ensure the extension is correct, this python script only takes .csv')
