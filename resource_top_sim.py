import cobra
import cobra.test
import numpy as np
import os, signal, json
import networkx as nx
import pandas as pd

# end products Ec_biomass_iJO1366_WT_53p95M Ec_biomass_iJO1366_core_53p95M

path = os.path.expanduser("~/GitHub/EcoEvoMet")
model = cobra.test.create_test_model("ecoli")
S_df = cobra.util.array.create_stoichiometric_matrix(model, array_type = 'DataFrame')

S = S_df.values
metabs = S_df.columns.values.tolist()
to_keep = [x for x in metabs if 'EX_' in x]
EX_Ec_biomass_iJO1366_WT_53p95M = []


def get_A():
    rxns_to_remove = ['Ec_biomass_iJO1366_WT_53p95M', 'Ec_biomass_iJO1366_core_53p95M']
    model = cobra.test.create_test_model("ecoli")
    S_df = cobra.util.array.create_stoichiometric_matrix(model, array_type = 'DataFrame')
    S_df.drop(rxns_to_remove,axis=1)
    S_df[S_df != 0] = 1
    S_np = S_df.values
    A_np = np.transpose(S_np).dot(S_np)
    A_df = pd.DataFrame(A_np, index=S_df.columns.values.tolist(), columns=S_df.columns.values.tolist())
    A_df[A_df != 0] = 1
    return A_df


def get_dist_A():
    A_df = get_A()
    A_nx = nx.from_numpy_matrix(A_df.values)
    d_A_nx = nx.all_pairs_shortest_path_length(A_nx)
    d_A_df = pd.DataFrame.from_dict(dict(d_A_nx))#, orient='index', columns=metabs)
    d_A_df = pd.DataFrame(data=d_A_df.values, index=metabs, columns=metabs)
    return d_A_df


def get_gene_rxm_df():
    gene_dict = {}
    rxn_dict = {}
    rxn_props = {}
    # Reading the json as a dict
    with open(path + '/iJO1366.json') as json_data:
        data = json.load(json_data)
        for gene in data['genes']:
            gene_dict[gene['id']] = gene['name']

    for reaction in model.reactions:
        reaction_split = str(reaction).split(':')[0]
        rxn_genes = list(model.reactions.get_by_id(reaction_split).genes)
        rxn_genes_names = [model.genes.get_by_id(str(x)).name for x in rxn_genes ]
        if len(rxn_genes_names) == 0:
            continue
        for rxn_genes_name in rxn_genes_names:
            rxn_dict[rxn_genes_name] = reaction_split

    S_df = cobra.util.array.create_stoichiometric_matrix(model, array_type = 'DataFrame')
    #S = S_df.values
    for column in S_df:
        column_np = S_df[column].values
        pos = column_np[column_np>0]
        neg = column_np[column_np<0]
        rxn_props[column] = len(column_np[column_np!=0])

    df_genes = pd.DataFrame(gene_dict.items(), columns=['bigg_id', 'gene_name'])
    df_rxn = pd.DataFrame(rxn_dict.items(), columns=['bigg_id', 'reaction'])
    df_merge = pd.merge(df_genes, df_rxn, on='bigg_id', how='outer')
    df_out = path + '/gene_rxn_table.txt'
    df_merge.to_csv(df_out, sep = '\t', index = True)


def get_m5_mutations():
    clade_hmm_states = {'A':0,'E':1,'FB':2,'FM':3, 'Fm':4,'PB':5,'PM':6,'Pm':7,'PB*':8}
    M_m5 = {}
    m_m5 = {}
    with open(path + '/gene_convergence_matrix.txt', "r") as file:
        for i, line in enumerate(file):
            if i == 0:
                continue
            line = line.strip().split(',')
            m5 = line[2].strip()
            if len(m5) == 0:
                continue
            muts = m5.split(';')
            gene = line[0]
            M_m5[gene] = {}
            m_m5[gene] = {}
            for mut in muts:
                mut_split = mut.split(':')
                status = int(mut_split[2])
                time = int(mut_split[0])
                M_m5[gene][time] = 0
                m_m5[gene][time] = 0
                if (status == 3) or (status == 6):
                    M_m5[gene][time] += 1

                elif (status == 4) or (status == 7):
                    m_m5[gene][time] += 1

                else:
                    continue
    print(M_m5)



# just modify your LTDE code...
class good_et_al:

    def __init__(self):
        self.populations = ['m6']

    def parse_convergence_matrix(self, filename):
        convergence_matrix = {}
        convergence_matrix_file = open(filename,"r")
        # Header line
        line = convergence_matrix_file.readline()
        populations = [item.strip() for item in line.split(",")[2:]]
        for line in convergence_matrix_file:
            items = line.split(",")
            gene_name = items[0].strip()
            length = float(items[1])
            convergence_matrix[gene_name] = {'length':length, 'mutations': {population: [] for population in populations}}
            for population, item in zip(populations,items[2:]):
                if item.strip()=="":
                    continue
                subitems = item.split(";")
                for subitem in subitems:
                    subsubitems = subitem.split(":")
                    mutation = (float(subsubitems[0]), float(subsubitems[1]), float(subsubitems[2]), float(subsubitems[3]))
                    convergence_matrix[gene_name]['mutations'][population].append(mutation)

        return convergence_matrix


    def reformat_convergence_matrix(self):#, mut_type = 'F'):
        conv_dict = self.parse_convergence_matrix(path + "/gene_convergence_matrix.txt")
        time_points = []
        new_dict = {}
        for gene_name, gene_data in conv_dict.items():
            for pop_name, mutations in gene_data['mutations'].items():
                for mutation in mutations:
                    time = int(mutation[0])
                    time_points.append(time)
        time_points = sorted(list(set(time_points)))
        for gene_name, gene_data in conv_dict.items():
            if gene_name not in new_dict:
                new_dict[gene_name] = {}
            for pop_name, mutations in gene_data['mutations'].items():
                if len(mutations) == 0:
                    continue
                if pop_name != "m6":
                    continue
                mutations.sort(key=lambda tup: tup[0])
                # keep only fixed mutations
                #clade_hmm_states = {'A':0,'E':1,'FB':2,'FM':3, 'Fm':4,'PB':5,'PM':6,'Pm':7,'PB*':8}
                #mutations = [list(x) for x in mutations]
                M_muts_F = [x for x in mutations if (int(x[2]) == 3)]
                M_muts_P = [x for x in mutations if (int(x[2]) == 6)]
                m_muts_F = [x for x in mutations if (int(x[2]) == 4)]
                m_muts_P = [x for x in mutations if (int(x[2]) == 7)]
                muts_lists = [M_muts_F, M_muts_P, m_muts_F, m_muts_P]
                if all(len(x) == 0 for x in muts_lists) == True:
                    continue
                mut_dict = {}
                if len(M_muts_F) > 0:
                    mut_dict['M_muts_F'] = list(M_muts_F[0])
                if len(M_muts_P) > 0:
                    mut_dict['M_muts_P'] = list(M_muts_P[0])
                if len(m_muts_F) > 0:
                    mut_dict['m_muts_F'] = list(m_muts_F[0])
                if len(m_muts_P) > 0:
                    mut_dict['m_muts_P'] = list(m_muts_P[0])
                for key, value in mut_dict.items():
                    value = list(value)
                    time = value[0]
                    if int(time) == 4250:
                        continue
                    pop_name = "m6_" + key[0]
                    if "F" in key:
                        remaining_time_points = time_points[time_points.index(time):]
                        for time_point in remaining_time_points:
                            pop_time = pop_name +'_' + str(int(time_point))
                            if pop_time not in new_dict[gene_name]:
                                new_dict[gene_name][pop_time] = 1
                            else:
                                new_dict[gene_name][pop_time] += 1
                    else:
                        pop_time = pop_name +'_' + str(int(mutation[0]))
                        if pop_time not in new_dict[gene_name]:
                            new_dict[gene_name][pop_time] = 1
                        else:
                            new_dict[gene_name][pop_time] += 1
        df = pd.DataFrame.from_dict(new_dict)
        df = df.fillna(0)
        df = df.loc[:, (df != 0).any(axis=0)]
        df_out = path + '/m6_gene_by_pop.txt'
        df.to_csv(df_out, sep = '\t', index = True)


def get_clade_dist():
    d = get_dist_A()
    print(d)







#good_et_al().reformat_convergence_matrix()
#get_m5_mutations()
#get_gene_rxm_df()

get_clade_dist()
