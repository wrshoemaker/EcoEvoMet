from __future__ import division
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
from scipy import stats

path = os.path.expanduser("~/GitHub/EcoEvoMet")


def CV_KDE(oneD_array):
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple


def resource_fig():
    df = pd.read_csv(path + '/alpha_data_final.txt', sep = '\t', header = 'infer', index_col = 0)

    KDE = CV_KDE(np.log10(df.alpha_mag.values))

    print("Bandwidth = " + str(KDE[2]))
    print("CV = " + str(np.std(np.log10(df.alpha_mag.values)) / np.mean(np.log10(df.alpha_mag.values)) ))

    fig = plt.figure()
    plt.plot(KDE[0], KDE[1], linewidth=3, alpha=0.5)#, label='bw=%.2f' % KDE[2])
    plt.xlabel('Magnitude of change in \nresource usage strategy ' + r'$\| \Delta\vec{\alpha}   \|, \; \mathrm{log}_{10}$', fontsize = 16)
    plt.ylabel('Frequency' , fontsize = 16)


    fig_name = path + '/resource_kde.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def pca_dist_fig():
    gene_by_pop_df = pd.read_csv(path + '/m6_gene_by_pop.txt', sep = '\t', header = 'infer', index_col = 0)
    X = gene_by_pop_df.values - np.mean(gene_by_pop_df.values, axis=0)
    pca = PCA(n_components=2)
    X_out = pca.fit_transform(X)
    df_pca = pd.DataFrame(data=X_out, index=gene_by_pop_df.index.values)
    gens =  [int(x.split('_')[2]) for x in df_pca.index.values.tolist() ]
    df_pca['Generations'] = gens
    m_df_pca = df_pca[df_pca.index.str.contains('_m_')]
    M_df_pca = df_pca[df_pca.index.str.contains('_M_')]

    df = pd.read_csv(path + '/time_dist.txt', sep = '\t', header = 'infer', index_col = 0)
    df = df[np.isfinite(df['Distance'])]

    fig = plt.figure()
    plt.subplot(211)
    plt.scatter(M_df_pca.Generations.values, M_df_pca[0].values, label="Ara-6 major lineage", color='#ff7f0e')
    plt.scatter(m_df_pca.Generations.values, m_df_pca[0].values, label="Ara-6 minor lineage", color='#1f77b4')
    plt.ylabel('PCA 1 (' + str(round(pca.explained_variance_ratio_[0],3)*100) + '%)' , fontsize = 12)
    plt.xlim(9500, 65000)
    plt.legend(loc='center right')
    plt.text(10000, 4, r"$\mathbf{a)}$", fontsize=12)


    plt.subplot(212)
    plt.plot(df.Time.values, df.Distance_m.values, color='#1f77b4')
    plt.plot(df.Time.values, df.Distance_M.values, color='#ff7f0e')
    plt.xlabel('Generation', fontsize = 16)
    #plt.ylabel('Between-clade \nreaction network ' , fontsize = 12)
    plt.ylabel( r'$\frac{\mathrm{between-lineage \; graph\;distance}}{\mathrm{within-lineage \; graph\;distance}}$', fontsize = 13)
    plt.hlines(y=1, xmin=9800, xmax=65000, color='k', linestyle=':', alpha = 0.8, zorder=1)
    plt.xlim(9500, 65000)
    plt.ylim(0.6, 2)
    plt.text(10000, 1.85, r"$\mathbf{b)}$", fontsize=12)

    fig_name = path + '/pca_dist_time.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def alpha_network():
    df_alpha = pd.read_csv(path + '/alpha_data_final.txt', sep = '\t', header = 'infer', index_col = 0)
    df_k = pd.read_csv(path + '/directed_rxn_k.txt', sep = '\t', header = 'infer', index_col = 0)
    df_c = pd.read_csv(path + '/directed_rxn_c.txt', sep = '\t', header = 'infer', index_col = 0)
    df_d = pd.read_csv(path + '/directed_rxn_d_ex_rxns.txt', sep = '\t', header = 'infer', index_col = 0)

    df_k_merged = df_alpha.merge(df_k, how='outer', left_index=True, right_index=True)
    df_c_merged = df_alpha.merge(df_c, how='outer', left_index=True, right_index=True)
    df_d_merged = df_alpha.merge(df_d, how='outer', left_index=True, right_index=True)

    df_d_merged = df_d_merged.dropna(subset=['dist_mean', 'alpha_mag'])


    fig = plt.figure()
    x = df_d_merged.dist_mean.values
    y = df_d_merged.alpha_mean.values
    plt.scatter(x, y, color='#ff7f0e')

    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print("slope is " + str(slope))
    print("r2-value is " + str(r_value **2))
    print("p-value is " + str(p_value))

    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')


    fig_name = path + '/test_alpha_network.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


#alpha_network()
#resource_fig()
pca_dist_fig()
