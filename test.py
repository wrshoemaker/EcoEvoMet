from __future__ import print_function
import cobra
import cobra.test
import numpy as np

N = np.asarray([[1, -1, -1, 0, 0, 0, 0],
                [0, -2, 2, 0, 4, 1, 1],
                [0, 2, 0, -1, -2, 0, 0],
                [-2, 4, 0, 0, -4, -1, 0],
                [0, 2, -2, -1, 3, 0, 0]])


model = cobra.test.create_test_model("ecoli")
#model = cobra.test.create_test_model("textbook")

N_df = cobra.util.array.create_stoichiometric_matrix(model, array_type = 'DataFrame')
N = N_df.values
metabs = N_df.columns.values.tolist()
to_keep = [x for x in metabs if 'EX_' in x]
e_coli_carbon = []
print(to_keep)

#test = cobra.flux_analysis.deletion.single_reaction_deletion(model, model.reactions[:5])
#print(model.objective.value)

def sim_alpha(model, to_keep):
    s_0 = cobra.sampling.sample(model, 100)
    s_0_in = s_0[to_keep]
    s_0_in_mean = s_0_in.mean(axis=0)

    for reaction in model.reactions[:5]:
        with model as model:
            reaction.knock_out()
            model.optimize()
            growth_rate = model.objective.value
            print(growth_rate)
            s = cobra.sampling.sample(model, 20)
            s_in = s[to_keep]
            s_in_mean = s_in.mean(axis=0)
            delta_alpha = s_in_mean - s_0_in_mean
            print(np.mean(delta_alpha) , np.linalg.norm(delta_alpha))
