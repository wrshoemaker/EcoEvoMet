from __future__ import division
import os
import numpy as np

import  matplotlib.pyplot as plt


path = os.path.expanduser("~/GitHub/EcoEvoMet")


def fig1():


    fig = plt.figure()
    plt.subplot(211, aspect='equal')
    #plt.plot([-1, 0, 1], [0, 0.5, 1])
    #plt.plot([-1, 0, 1], [-1, 0, 1])
    plt.axvline(-0.6, color='grey', linestyle='--')
    plt.axvline(0.6, color='grey', linestyle='--')

    plt.axvspan(-1.2, -0.6, alpha=0.5, color='red',zorder=1)
    plt.axvspan(0.6, 1.2, alpha=0.5, color='red',zorder=1)
    plt.axvspan(-0.6, 0.6, alpha=0.5, color='royalblue',zorder=1)

    plt.axvline(-0.6, color='k', linestyle='--',zorder=2)
    plt.axvline(0.6, color='k', linestyle='--',zorder=2)

    plt.plot([0.6, 1.2], [1,1], c='k',zorder=3)
    plt.plot([-1.2, -0.6], [-1,-1], c='k',zorder=3)
    plt.plot([-0.6, 0.6], [-1, 1], c='k',zorder=3)

    plt.plot([0, 0], [-1.2, 0], c='k',linestyle='--',zorder=3)
    plt.plot([-1.2, 0], [0, 0], c='k',linestyle='--',zorder=3)
    plt.scatter([0],[0],c='k',alpha=1)

    plt.yticks([-1,-0.5,0,0.5,1], [0,None, None, None,1], rotation='horizontal')

    plt.text(0.5, -1.39, r'$X_{max}$')
    plt.text(-0.67, -1.39, r'$X_{min}$')
    plt.text(-0.93, 0.1, r'$f_{0}^{*}$')

    plt.text(-1.15, 1.01, "a)")

    plt.xlim([-1.2, 1.2])
    plt.ylim([-1.2, 1.2])

    plt.ylabel('Equilib. frequency, ' + r'$f^{*}(\Delta X)$')


    # plot 2!
    plt.subplot(212, aspect='equal')
    plt.xlabel('Fitness difference, ' + r'$\Delta X$')

    plt.xlim([-1.2, 1.2])
    plt.ylim([-1.2, 1.2])

    plt.yticks([-0.8, 0.8], [r'$-$', r'$+$'], rotation='horizontal')

    plt.plot([-1, 1], [-0.8, 0.8], c='#87CEEB')
    plt.plot([-1, 1], [0.8, -0.8], c='#FF6347')

    plt.scatter([-1, 1],[-0.8, 0.8],c='#87CEEB',alpha=1, label="Glucose muts.", zorder=100)
    plt.scatter([-1, 1],[0.8, -0.8],c='#FF6347',alpha=1, label="Acetate muts.", zorder=101)

    plt.text(-1.15, 1.01, "b)")

    plt.legend(loc='upper center',prop={'size': 6.5})

    plt.ylabel('Selection coefficient, ' + r'$s$')



    fig_name = path + '/fitness.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



fig1()
