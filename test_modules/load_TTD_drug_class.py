import pandas as pd
import cmap.io.gmt as gmt
'parse the TTD file for drug annotations'

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

    def load_drugbank_by_gene(self,group_by_action=True):
        drugFile = '/xchip/cogs/hogstrom/notes/TTD_annotations/drugbanktarget_pert_id.csv'
        drugLabels = pd.io.parsers.read_csv(drugFile)
        drugLabels['action_calc'] = drugLabels['action_calc'].str.replace("other/unknown","other")
        # specify directionality
        brdLstLst = drugLabels.pert_ids.str.split('|').tolist()
        drugLabels['pert_ids'] = brdLstLst
        # group by gene and action
        ttd_cp_dict = {}
        if group_by_action:
            grouped = drugLabels.groupby(['gene','action_calc'])
            for group in grouped:
                groupName = group[1]['gene'].values[0]
                groupCat = group[1]['action_calc'].values[0]
                cpLstLst = group[1]['pert_ids'].values
                cpLst = [item for sublist in cpLstLst for item in sublist] 
                cpLst = list(set(cpLst)) # make sure list is unique
                if groupCat == '-666':
                    ttd_cp_dict[groupName] = cpLst
                else:
                    ttd_cp_dict[groupName+'-'+groupCat] = cpLst
        else:
            grouped = drugLabels.groupby(['gene'])
            for group in grouped:
                groupName = group[1]['gene'].values[0]
                # groupCat = group[1]['action_calc'].values[0]
                cpLstLst = group[1]['pert_ids'].values
                cpLst = [item for sublist in cpLstLst for item in sublist] 
                cpLst = list(set(cpLst)) # make sure list is unique
                ttd_cp_dict[groupName] = cpLst
        # group only be gene
        grouped = drugLabels.groupby(['gene'])
        gene_cp_dict = {}
        for ig, group in enumerate(grouped):
            # if ig > 40:
            #     continue
            if group[0] == '-666':
                continue    
            groupName = group[1]['gene'].values[0]
            cpLstLst = group[1]['pert_ids'].values
            cpLst = [item for sublist in cpLstLst for item in sublist] 
            cpLst = list(set(cpLst)) # make sure list is unique
            gene_cp_dict[groupName] = cpLst
        return ttd_cp_dict

    def load_clique_set_n69(self):
        ''' 
        load drug label set currated by Rajiv

        '''
        #load in clique annotations and matrix
        cFile = '/xchip/cogs/sig_tools/sig_cliquescore_tool/sample/cp_clique_n69/clique.gmt'
        cliqueGMT = gmt.read(cFile)
        cliqFrm = pd.DataFrame(cliqueGMT)
        # set grouping structures 
        pclDict = {}
        for x in cliqFrm.iterrows():
            pclDict[x[1]['id']] = set(x[1]['sig'])
        return pclDict