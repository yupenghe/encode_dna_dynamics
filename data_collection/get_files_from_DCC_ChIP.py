#!/usr/bin/env python2.7
import requests, json
import subprocess
import os

fullname2short = {"embryonic facial prominence":"CF",
		  "forebrain":"FB",
		  "heart":"HT",
		  "hindbrain":"HB",
		  "intestine":"IT",
		  "kidney":"KD",
		  "limb":"LM",
		  "liver":"LV",
		  "lung":"LG",
		  "midbrain":"MB",
		  "neural tube":"NT",
		  "stomach":"ST"}


# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}


all_datasets = {"E10_5":['ENCSR349UOB',#CF
                         'ENCSR714JXE',#HT
                         'ENCSR228WLP',#FB
                         'ENCSR050BJO',#HB
                         'ENCSR893YWL',#MB
                         'ENCSR386FAC'],#LM
                "E11_5":['ENCSR800JXR',#CF
                         'ENCSR843IAS',#MB
                         'ENCSR501OPC',#HB
                         'ENCSR415TUB',#FB
                         'ENCSR215ZYV',#NT
                         'ENCSR283NCE',#LM
                         'ENCSR016LTR',#HT
                         'ENCSR231EPI'],#LV
                "E12_5":['ENCSR273IZZ',#NT
                         'ENCSR605BEG',#HB
                         'ENCSR749IWM',#CF
                         'ENCSR192PSV',#LV
                         'ENCSR589KOM',#LM
                         'ENCSR296YDZ',#FB
                         'ENCSR524JAQ',#HT
                         'ENCSR553AAA'#MB
                ],
                "E13_5":['ENCSR629UDB',#HB
                         'ENCSR832HYG',#NT
                         'ENCSR911FIM',#FB
                         'ENCSR424HQH',#LM
                         'ENCSR266TXT',#LV
                         'ENCSR792SPH',#CF
                         'ENCSR107EDN',#HT
                         'ENCSR909IQE'],#MB
                "E14_5":['ENCSR026SJT',#IT
                         'ENCSR548TUQ',#HT ##Extra attension
                         'ENCSR250AOK',#ST
                         'ENCSR840QVS',#KD
                         'ENCSR944NEQ',#LG
                         'ENCSR006XBO',#CF
                         'ENCSR295OVP',#NT
                         'ENCSR705YBY',#LM
                         'ENCSR649KUM',#MB
                         'ENCSR514BWH',#HB
                         'ENCSR894URD',#FB
                         'ENCSR248SXA'],#LV
                "E15_5":['ENCSR189YPN',#KD
                         'ENCSR874LCF',#NT
                         'ENCSR438IXZ',#FB
                         'ENCSR871WIE',#LV
                         'ENCSR693OUQ',#IT
                         'ENCSR486MYL',#LM
                         'ENCSR881OKK',#LG
                         'ENCSR379NJH',#ST
                         'ENCSR356VGB',#HB
                         'ENCSR899FPM',#CF
                         'ENCSR520JSM',#HT
                         'ENCSR453LDG'],#MB
                "E16_5":['ENCSR483SIP',#MB
                         'ENCSR774QPF',#HB
                         'ENCSR502XBO',#FB
                         'ENCSR709SFR',#KD
                         'ENCSR287WEJ',#LG
                         'ENCSR967ZGY',#HT
                         'ENCSR124NGS',#IT
                         'ENCSR826FLO',#LV
                         'ENCSR827DHS'],#ST
                "P0":['ENCSR767BSH',#KD
                      'ENCSR571BOJ',#ST
                      'ENCSR420BKD',#IT
                      'ENCSR203RUV',#LG
                      'ENCSR687SNT',#LV
                      'ENCSR049ZYT',#HT
                      'ENCSR219XIH',#FB
                      'ENCSR803VSL',#HB
                      'ENCSR022KAT']#MB
}

os.system("mkdir -p raw_ChIP/")

for age in all_datasets.iterkeys():
    for accession in all_datasets[age]:
        # This URL locates the ENCODE biosample with accession number ENCSR800JXR
        URL = "https://www.encodeproject.org/reference-epigenomes/"+accession+"/?frame=embedded"

        # GET the object
        response = requests.get(URL, headers=HEADERS)

        # Extract the JSON response as a python dict
        response_json_dict = response.json()

        for ind_related_dataset in range(len(response_json_dict[u'related_datasets'])):
            related_dataset = response_json_dict["related_datasets"][ind_related_dataset]
            # Get information
            if related_dataset['assay_term_name'] != "ChIP-seq":
                continue
            short_name = fullname2short[related_dataset['biosample_term_name']]

            ## Filter
            #if not (short_name == "HB"):
            #    continue
            
            mark = related_dataset['target']['label']

            #if not (mark == "H3K27ac"):
            if not (mark == "Control"):
                continue            

            if mark == "H3K9me3" or mark == "H3K36me3":
                continue;
            
            #print("_".join([age,short_name,mark]))
            # Download the fastq files and label them with biological replication number
            for ind_file in range(len(related_dataset['files'])):
                #if related_dataset['files'][ind_file]['file_format'] != 'fastq':
                if related_dataset['files'][ind_file]['file_format'] != 'bam':
                    continue
                if related_dataset['files'][ind_file]['output_type'] != 'alignments':
                    continue
                if related_dataset['files'][ind_file]['assembly'] != "mm10-minimal":
                    continue
                bio_rep = related_dataset['files'][ind_file]['biological_replicates']
                output_prefix = "_".join([age,short_name,str(bio_rep[0]),mark])
                #if len(bio_rep)!=1:
                #print output_prefix
                # Download
                url = 'https://www.encodeproject.org'+related_dataset['files'][ind_file]['href']
                output_filename = "raw_ChIP/" + output_prefix + "." + str(ind_file) + ".bam"
                subprocess.check_call(['curl',
                                       '-RL',
                                       url,
                                       "-o",
                                       output_filename
                ])
