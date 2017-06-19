import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
from xml.etree import ElementTree as ET
from matplotlib import numpy as np

def plotLIE(trainSet):
    
    x=np.array([ i['Gexp'] for i in trainSet ])
    y=np.array([ i['Gcalc'] for i in trainSet ])

    minAxi=min([x.min(), y.min()])
    maxAxi=max([x.max(), y.max()])
    minAxi= 5*((minAxi//5))
    maxAxi= 5*((maxAxi//5)+1)
    
    fig = plt.figure()       
    ax = fig.add_subplot(111, aspect='equal')
    plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=20)
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=20)

    plt.xlabel('$\Delta$G$_{obs}$', fontsize=20, color='black')
    plt.ylabel(r'$\Delta$G$_{calc}$', fontsize=20, color='black')
    plt.axis([minAxi,maxAxi, minAxi, maxAxi])
    plt.plot(x,y, 'ro')
    plt.plot([minAxi,maxAxi],[minAxi,maxAxi],color='black', linestyle='-')
    plt.plot([minAxi,maxAxi-5],[minAxi+5,maxAxi],color='black', linestyle='-.')
    plt.plot([minAxi+5,maxAxi],[minAxi,maxAxi-5],color='black', linestyle='-.')

    for i,(x,y) in enumerate(zip(x,y)):
        ax.annotate(i+1,xy=(x,y),xytext=(5,0),textcoords='offset points')
    
    plt.tight_layout()
    imgdata = BytesIO()
   
    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)  # rewind the data

    return imgdata


def fixSVG(imgdata,h=200):
    #remove xml header:
    figdata_svg = '<svg' + imgdata.getvalue().split('<svg')[1]
    
    #return unicode(figdata_svg,'utf-8')  
    ET.register_namespace("","http://www.w3.org/2000/svg")
    ET.register_namespace("xlink","http://www.w3.org/1999/xlink")

    # fix size    
    root=ET.fromstring(figdata_svg)
    root.attrib['height']='%d'%h
    root.attrib['preserveAspectRatio']="xMinYMin meet"
    
    #export to unicode string
    figout=BytesIO()
    tree=ET.ElementTree(root)

    tree.write(figout, encoding='UTF-8')

    return unicode(figout.getvalue(),'utf-8') #svg_dta