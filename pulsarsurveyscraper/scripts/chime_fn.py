#!/usr/bin/env python

import numpy as np
import csv
import os
def load_clusters(candidates):
    #loads the candidates, should be a .npy file, this just does regular CHIME clusters without header localisation
    try:
        sp = np.load(candidates,allow_pickle=1).tolist()
    except:
        sp = np.load(candidates,allow_pickle=1)['data'].tolist()
    ul = set(sp.dbscan_labels)
    return_arr = []
    for l in ul:
        indexes = l==sp.dbscan_labels
        ra = sp.pos_ra_deg[indexes]
        dec = sp.pos_dec_deg[indexes]
        dm = sp.dm[indexes]
        times = sp.event_time[indexes]
        ara = np.mean(ra)
        adec = np.mean(dec)
        adm = np.mean(dm)
        era = 2.2/np.cos(np.radians(adec))
        edec = 0.5
        edm = max(dm)-min(dm)
        return_arr.append({'cc':l,'ara':ara,'adec':adec,'adm':adm,'era':era,'edec':edec,'edm':edm,'times':times})
    return np.array(return_arr)

def load_fine_tune(candidates):
    #loads the candidates, should be a .npy file
    sp = np.load(candidates,allow_pickle=1).tolist()
    fine_tune_dict = []
    for cluster in sp.fine_tune:
        #loop through all the fune tuned clusters
        ul = set(cluster.labels)
        for l in ul:
            if l>-2:
                mask = (l==cluster.labels)
                if len(cluster.labels) >1:
                    #use error weighted average
                    ftc_ra = cluster.localised_pos_ra_deg[mask]
                    ftc_dec = cluster.localised_pos_dec_deg[mask]
                    ftc_dm = cluster.dm[mask]
                    ftc_ra_err = cluster.ra_error[mask]
                    ftc_dec_err = cluster.dec_error[mask]
                    ftc_dm_err = cluster.dm_error[mask]

                    #inflate the errors by a factor of 2
                    ftc_ra_err = ftc_ra_err*2
                    ftc_dec_err = ftc_dec_err*2

                    #do an error weighted average
                    ara = np.average(ftc_ra,weights=1/np.square(ftc_ra_err))
                    adec = np.average(ftc_dec,weights=1/np.square(ftc_dec_err))
                    adm = np.average(ftc_dm,weights=1/np.square(ftc_dm_err))
                    #report the inflated errors, do this because unsure of systematics
                    era = np.average(ftc_ra_err)
                    edec = np.average(ftc_dec_err)
                    edm = np.average(ftc_dm_err)

                    fine_tune_dict.append({'c':cluster.c_num,'cc':l,'ara':ara,'adec':adec,'adm':adm,'era':era,'edec':edec,'edm':edm})
    return np.array(fine_tune_dict)
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
                new_item = {row[0]:np.array([row[4],row[5],row[3]])}
                my_new_sources.update(new_item)
        return my_new_sources
    else:
        print('please ensure the extension is correct, this python script only takes .csv')
