import pickle, os

path = os.path.expanduser("~/GitHub/EcoEvoMet")



with open(path + '/iJO1366.pickle', 'rb') as pickle_file:
    d = pickle.load(pickle_file)

    print(dir(d))

    print(d.description)
