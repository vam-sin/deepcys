import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import numpy as np
import string
import math
import sys
import re


input_file = sys.argv[1]
A = []
for index, line in enumerate(open(input_file)):
    A.append(line.split()[1])
myset = set(A)
chains = list(myset)
no_of_chains = len(chains)

if no_of_chains < 5:
    xDict = {}
    yDict = {}
    for item in chains:
        listx = []
        listy = []
        with open(input_file) as f:
            for line in f:
                tok = line.split()
                if tok[1] == item:
                    listx.append(re.sub('[a-zA-Z]+','',tok[2]))
                    listy.append(tok[9])
        f.close()
        xDict[item] = np.asarray(listx).astype(np.float)
        yDict[item] = np.asarray(listy).astype(np.float)    

    number_of_subplots = no_of_chains

    m = cm.ScalarMappable(cmap=cm.jet)

    for v,item in zip(range(number_of_subplots),chains):
        v = v+1
        m.set_array(yDict[item])
        ax1 = plt.subplot(number_of_subplots,1,v)
        #ax1.set_ylim(-0.5,1.5)
        l=ax1.scatter(xDict[item],yDict[item],c=yDict[item], cmap=plt.cm.get_cmap("gist_rainbow"), edgecolors='None',alpha=0.75 ,label='Colour based on rHpy values')
        #l = ax1.plot(xDict[item],yDict[item])
        ax1.set_title('chain: %s'%(item))
        ax1.set_xlabel('Residue Numbers')
        ax1.set_ylabel('rHpy', fontsize=15)
        
    plt.tight_layout(rect = [0,0.01,1,0.95])

    protein = input_file[:-9]

    plt.suptitle('Graphical representation of Microenvironment of %s protein'%(protein.split('/')[-1:][0]))

    # fig = matplotlib.pyplot.gcf()
    # fig.set_size_inches(19.5,10.5)

    #plt.show()
    plt.savefig("%s.png"%(protein),bbox_inches='tight',dpi=200)

elif 4 < no_of_chains < 10: 

    subplots_adjust(hspace=0.75)
    subplots_adjust(top=0.9)
    subplots_adjust(bottom=-1.2)

    xDict = {}
    yDict = {}
    for item in chains:
    	listx = []
    	listy = []
    	with open(input_file) as f:
    		for line in f:
    			tok = line.split()
    			if tok[1] == item:
    				listx.append(re.sub('[a-zA-Z]+','',tok[2]))
    				listy.append(tok[9])
    	f.close()
    	xDict[item] = np.asarray(listx).astype(np.float)
    	yDict[item] = np.asarray(listy).astype(np.float)	

    number_of_subplots = no_of_chains

    m = cm.ScalarMappable(cmap=cm.jet)

    for v,item in zip(range(number_of_subplots),chains):
        v = v+1
        m.set_array(yDict[item])
        ax1 = plt.subplot(number_of_subplots,3,v)
        #ax1.set_ylim(-0.5,1.5)
        l=ax1.scatter(xDict[item],yDict[item],c=yDict[item], cmap=plt.cm.get_cmap("gist_rainbow"), edgecolors='None',alpha=0.75 ,label='Colour based on rHpy values')
        #l = ax1.plot(xDict[item],yDict[item])
        ax1.set_title('chain: %s'%(item))
        ax1.set_xlabel('Residue Numbers')
        ax1.set_ylabel('rHpy', fontsize=15)
        #plt.colorbar(l)
        
    #plt.tight_layout()
    #rect = [0,0.01,1,0.95]


    protein = input_file[:-9]

    plt.suptitle('Graphical representation of Microenvironment of %s protein'%(protein.split('/')[-1:][0]))

    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(19.5,10.5)

    #plt.show()
    plt.savefig("%s.png"%(protein),bbox_inches='tight')

elif 9 < no_of_chains < 13:
    subplots_adjust(hspace=0.75)
    subplots_adjust(top=0.9)
    subplots_adjust(bottom=-1.5)

    xDict = {}
    yDict = {}
    for item in chains:
        listx = []
        listy = []
        with open(input_file) as f:
            for line in f:
                tok = line.split()
                if tok[1] == item:
                    listx.append(re.sub('[a-zA-Z]+','',tok[2]))
                    listy.append(tok[9])
        f.close()
        xDict[item] = np.asarray(listx).astype(np.float)
        yDict[item] = np.asarray(listy).astype(np.float)    

    number_of_subplots = no_of_chains

    m = cm.ScalarMappable(cmap=cm.jet)

    for v,item in zip(range(number_of_subplots),chains):
        v = v+1
        m.set_array(yDict[item])
        ax1 = plt.subplot(number_of_subplots,3,v)
        #ax1.set_ylim(-0.5,1.5)
        l=ax1.scatter(xDict[item],yDict[item],c=yDict[item], cmap=plt.cm.get_cmap("gist_rainbow"), edgecolors='None',alpha=0.75 ,label='Colour based on rHpy values')
        #l = ax1.plot(xDict[item],yDict[item])
        ax1.set_title('chain: %s'%(item))
        ax1.set_xlabel('Residue Numbers')
        ax1.set_ylabel('rHpy', fontsize=15)
        #plt.colorbar(l)
        
    #plt.tight_layout()
    #rect = [0,0.01,1,0.95]


    protein = input_file[:-9]

    plt.suptitle('Graphical representation of Microenvironment of %s protein'%(protein.split('/')[-1:][0]))

    # manager = plt.get_current_fig_manager()
    # manager.resize(*manager.window.maxsize())
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(19.5,10.5)
    #plt.show()
    plt.savefig("%s.png"%(protein),bbox_inches='tight')

else :
    subplots_adjust(hspace=0.75)
    subplots_adjust(top=0.9)
    subplots_adjust(bottom=-1.2)

    xDict = {}
    yDict = {}
    for item in chains:
        listx = []
        listy = []
        with open(input_file) as f:
            for line in f:
                tok = line.split()
                if tok[1] == item:
                    listx.append(re.sub('[a-zA-Z]+','',tok[2]))
                    listy.append(tok[9])
        f.close()
        xDict[item] = np.asarray(listx).astype(np.float)
        yDict[item] = np.asarray(listy).astype(np.float)    

    number_of_subplots = no_of_chains

    m = cm.ScalarMappable(cmap=cm.jet)

    for v,item in zip(range(number_of_subplots),chains):
        v = v+1
        m.set_array(yDict[item])
        ax1 = plt.subplot(number_of_subplots,5,v)
        #ax1.set_ylim(-0.5,1.5)
        l=ax1.scatter(xDict[item],yDict[item],c=yDict[item], cmap=plt.cm.get_cmap("gist_rainbow"), edgecolors='None',alpha=0.75 ,label='Colour based on rHpy values')
        #l = ax1.plot(xDict[item],yDict[item])
        ax1.set_title('chain: %s'%(item))
        ax1.set_xlabel('Residue Numbers')
        ax1.set_ylabel('rHpy', fontsize=15)
        #plt.colorbar(l)
        
    #plt.tight_layout()
    #rect = [0,0.01,1,0.95]


    protein = input_file[:-9]

    plt.suptitle('Graphical representation of Microenvironment of %s protein'%(protein.split('/')[-1:][0]))

    # manager = plt.get_current_fig_manager()
    # manager.resize(*manager.window.maxsize())
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(19.5,10.5)
    #plt.show()
    plt.savefig("%s.png"%(protein),bbox_inches='tight')
