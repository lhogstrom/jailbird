import pandas as pd

class label_loader:
    """load in drug labels"""
    def load_TTD(self):
        # ## TTD pcl class list - updated grpSetept/2013
        drugFile = '/xchip/cogs/hogstrom/notes/TTD_annotations/TTD_targetID_pert_id_mapping_20130927.txt'
        drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
        # replace ugly characters
        drugLabels['Name'] = drugLabels['Name'].str.replace("/","-")
        drugLabels['Name'] = drugLabels['Name'].str.replace(" ","_")
        drugLabels['Name'] = drugLabels['Name'].str.replace("&","_")
        drugLabels['Name'] = drugLabels['Name'].str.replace("?","_")
        drugLabels['Name'] = drugLabels['Name'].str.replace("(","_")
        drugLabels['Name'] = drugLabels['Name'].str.replace(")","_")
        drugLabels['Name'] = drugLabels['Name'].str.replace("'","_")
        # specify directionality
        brdLstLst = drugLabels.pert_ids_merged.str.split('|').tolist()
        drugLabels['pert_ids_merged'] = brdLstLst
        grouped = drugLabels.groupby(['TTD_target_ID','Category'])
        ttd_cp_dict = {}
        for ig, group in enumerate(grouped):
            # if ig > 20: # shorten diction 
            #     continue
            groupName = group[1]['Name'].values[0]
            groupCat = group[1]['Category'].values[0]
            cpLstLst = group[1]['pert_ids_merged'].values
            cpLst = [item for sublist in cpLstLst for item in sublist] 
            cpLst = list(set(cpLst)) # make sure list is unique
            ttd_cp_dict[groupName+'-'+groupCat] = cpLst
        #set inames
        inameDict = {}
        for x in drugLabels.iterrows():
            brds = x[1]['pert_ids_merged']
            iname = x[1]['Drug_name']
            if len(brds) > 1:
                for brd in brds:
                    inameDict[brd] = iname
            else:
                inameDict[brds[0]] = iname
        return ttd_cp_dict
