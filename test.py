import cobra
import cobra.test
import numpy as np
import os, signal, random
from timeit import default_timer as timer
from datetime import timedelta


path = os.path.expanduser("~/GitHub/EcoEvoMet")

model = cobra.test.create_test_model("ecoli")
#model = cobra.test.create_test_model("textbook")

N_df = cobra.util.array.create_stoichiometric_matrix(model, array_type = 'DataFrame')
N = N_df.values
metabs = N_df.columns.values.tolist()
to_keep = [x for x in metabs if 'EX_' in x]
e_coli_carbon = []

#test = cobra.flux_analysis.deletion.single_reaction_deletion(model, model.reactions[:5])
def handler(signum, frame):
    print("Optimization took > 4 min")
    raise Exception("end of time")

to_skip = ['NTD11pp']

def sim_alpha(model, to_keep, iter=100):
    df_out = open(path + '/alpha_data_2.txt', 'w')
    header = ['reaction', 'alpha_mean', 'alpha_mag' ,'growth_rate']
    df_out.write('\t'.join(header) + '\n')
    s_0 = cobra.sampling.sample(model, iter)
    s_0_in = s_0[to_keep]
    s_0_in_mean = s_0_in.mean(axis=0)
    count = len([str(x) for x in model.reactions]) - len(to_keep)
    signal.signal(signal.SIGALRM, handler)
    #reactions = list(model.reactions)
    #random.shuffle(reactions)
    # 1000 +
    #for i, reaction in enumerate(model.reactions[676: 2258-666]):
    for i, reaction in enumerate(model.reactions):
        # 4 min * 60 sec/min = 240 sec
        signal.alarm(240)
        try:
            start = timer()
            reaction_split = str(reaction).split(':')[0]
            if reaction_split in to_skip:
                continue
            #NTD11
            if reaction_split in to_keep:
                continue
            count -= 1
            with model as model:
                reaction.knock_out()
                model.optimize()
                growth_rate = model.objective.value
                try:
                    s = cobra.sampling.sample(model, iter)
                except:
                    continue
                else:
                    s_in = s[to_keep]
                    s_in_mean = s_in.mean(axis=0)
                    delta_alpha = s_in_mean - s_0_in_mean
                    out_line = [reaction_split, str(np.mean(delta_alpha)), str(np.linalg.norm(delta_alpha)), str(growth_rate)]
                    df_out.write('\t'.join(out_line) + '\n')

            end = timer()
            print(reaction_split + " complete, " + str(count) + " reactions to go!, " +  str(timedelta(minutes=end-start)) + " min." )

        except exc:
            count -= 1
            print(exc)

    df_out.close()



sim_alpha(model, to_keep)
