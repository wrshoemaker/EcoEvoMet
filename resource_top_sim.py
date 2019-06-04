from __future__ import division
import cobra
import cobra.test
import numpy as np
import os, json, itertools
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


def get_rxns_A():
    rxns_to_remove = ['Ec_biomass_iJO1366_WT_53p95M', 'Ec_biomass_iJO1366_core_53p95M']
    model = cobra.test.create_test_model("ecoli")
    S_df = cobra.util.array.create_stoichiometric_matrix(model, array_type = 'DataFrame')
    S_df.drop(rxns_to_remove, axis=1)
    S_df[S_df != 0] = 1
    S_np = S_df.values
    A_np = np.transpose(S_np).dot(S_np)
    A_df = pd.DataFrame(A_np, index=S_df.columns.values.tolist(), columns=S_df.columns.values.tolist())
    A_df[A_df != 0] = 1
    return A_df


def get_dist_A():
    A_df = get_rxns_A()
    A_nx = nx.from_numpy_matrix(A_df.values)
    d_A_nx = nx.all_pairs_shortest_path_length(A_nx)
    d_A_df = pd.DataFrame.from_dict(dict(d_A_nx)) #, orient='index', columns=metabs)
    d_A_df = pd.DataFrame(data=d_A_df.values, index=metabs, columns=metabs)
    df_out = path + '/d_rxns_iJO1366.txt'
    d_A_df.to_csv(df_out, sep = '\t', index = True)



def get_directed_rxn_A():
    #rxns_to_remove = ['Ec_biomass_iJO1366_WT_53p95M', 'Ec_biomass_iJO1366_core_53p95M']
    model = cobra.test.create_test_model("ecoli")
    S_df = cobra.util.array.create_stoichiometric_matrix(model, array_type = 'DataFrame')
    rxns = S_df.columns.tolist()
    df_direct_A = pd.DataFrame(index=S_df.columns.values, columns=S_df.columns.values)
    df_direct_A = df_direct_A.fillna(0)
    rxn_pairs = list(itertools.combinations(rxns, 2))
    num_pairs = len(rxn_pairs)
    for rxn_pair in rxn_pairs:
        if num_pairs % 1000 == 0:
            print(str(num_pairs) + " pairs to go!")
        num_pairs -= 1
        rxn_pair = list(rxn_pair)
        if rxn_pair[0] == rxn_pair[1]:
            continue
        S_df_rxn_pair = S_df[rxn_pair]
        S_df_rxn_pair = S_df_rxn_pair[(S_df_rxn_pair != 0).all(1)]
        if S_df_rxn_pair.shape[0] == 0:
            continue
        S_df_rxn_pair[S_df_rxn_pair > 0] = 1
        S_df_rxn_pair[S_df_rxn_pair < 0] = -1
        # rxn1 => rxn2
        rxn1_rxn2 = S_df_rxn_pair.loc[(S_df_rxn_pair[rxn_pair[0]] == -1) & (S_df_rxn_pair[rxn_pair[1]] == 1)]
        # rxn2 => rxn1
        rxn2_rxn1 = S_df_rxn_pair.loc[(S_df_rxn_pair[rxn_pair[1]] == -1) & (S_df_rxn_pair[rxn_pair[0]] == 1)]
        if rxn1_rxn2.shape[0] > 0:
            df_direct_A.at[rxn_pair[0], rxn_pair[1]] = 1
        if rxn2_rxn1.shape[0] > 0:
            df_direct_A.at[rxn_pair[1], rxn_pair[0]] = 1

    df_out = path + '/directed_rxn_A.txt'
    df_direct_A.to_csv(df_out, sep = '\t', index = True)


def get_directed_rxn_d():
    A_df = pd.read_csv(path + '/directed_rxn_A.txt', sep = '\t', header = 'infer', index_col = 0)
    A_nx = nx.from_numpy_matrix(A_df.values)
    d_A_nx = nx.all_pairs_shortest_path_length(A_nx)
    d_A_df = pd.DataFrame.from_dict(dict(d_A_nx)) #, orient='index', columns=metabs)
    d_A_df = pd.DataFrame(data=d_A_df.values, index=metabs, columns=metabs)
    df_out = path + '/d_rxns_iJO1366.txt'
    d_A_df.to_csv(df_out, sep = '\t', index = True)


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



#def get_pair_dist():

def get_clade_dist():
    d_df = pd.read_csv(path + '/d_iJO1366.txt', sep = '\t', header = 'infer', index_col = 0)
    gene_df = pd.read_csv(path + '/m6_gene_by_pop.txt', sep = '\t', header = 'infer', index_col = 0)
    rxn_df = pd.read_csv(path + '/gene_rxn_table.txt', sep = '\t', header = 'infer', index_col = 0)
    rxn_genes = rxn_df.gene_name.tolist()
    rxn_genes = [x for x in rxn_genes if x in gene_df.columns.values]
    rxns = []
    rxn_series = []
    for rxn_gene in rxn_genes:
        rxn = rxn_df.loc[rxn_df['gene_name'] == rxn_gene, 'reaction'].iloc[0]
        rxns.append(rxn)
        gene_column = gene_df[rxn_gene ]
        rxn_column = gene_column.rename(rxn)
        rxn_series.append(rxn_column)
    # no duplicate rxns
    rxn_time_df = pd.concat(rxn_series, axis=1, keys=[r.name for r in rxn_series])
    major_clade = rxn_time_df.loc['m6_M_62750']
    minor_clade = rxn_time_df.loc['m6_m_62750']
    major_clade = major_clade[major_clade != 0]
    minor_clade = minor_clade[minor_clade != 0]
    major_clade_rxns = major_clade.index.values
    minor_clade_rxns = minor_clade.index.values

    def get_mean_dist(list1, list2):
        dists = []
        for i, rxn_i in enumerate(list1):
            for j, rxn_j in enumerate(list2):
                if (i < j) or (rxn_i != rxn_j):
                    dists.append(d_df.loc[rxn_i][rxn_j])

        return np.mean(dists)

    time_points = list(set([int(x.split('_')[2]) for x in rxn_time_df.index.tolist()]))
    time_points.remove(11750)
    time_points.remove(11250)
    time_points.remove(6250)
    time_points.sort()
    dist_diff = []
    dist_diff_m = []
    dist_diff_M = []
    for time_point in time_points:
        rxn_timepoint_m = rxn_time_df.loc['m6_m_' + str(time_point)]
        rxn_timepoint_M = rxn_time_df.loc['m6_M_' + str(time_point)]
        rxn_timepoint_m = rxn_timepoint_m[rxn_timepoint_m !=0]
        rxn_timepoint_M = rxn_timepoint_M[rxn_timepoint_M !=0]
        rxn_timepoint_m_rxns = rxn_timepoint_m.index.values
        rxn_timepoint_M_rxns = rxn_timepoint_M.index.values

        m_within = get_mean_dist(rxn_timepoint_m_rxns, rxn_timepoint_m_rxns)
        M_within = get_mean_dist(rxn_timepoint_M_rxns, rxn_timepoint_M_rxns)
        m_M_between = get_mean_dist(rxn_timepoint_m_rxns, rxn_timepoint_M_rxns)
        dist_diff.append(m_M_between / np.mean([m_within, M_within]))
        dist_diff_m.append(m_M_between / m_within)
        dist_diff_M.append(m_M_between / M_within)

    df_dists = pd.DataFrame(zip(time_points, dist_diff, dist_diff_m, dist_diff_M), columns=['Time', 'Distance', 'Distance_m', 'Distance_M'])
    df_out = path + '/time_dist.txt'
    df_dists.to_csv(df_out, sep = '\t', index = True)







#good_et_al().reformat_convergence_matrix()
#get_m5_mutations()
#get_gene_rxm_df()
#get_clade_dist()

get_directed_rxn_A()
