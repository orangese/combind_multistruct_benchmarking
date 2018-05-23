
from matplotlib import pyplot as plt
import numpy as np

def stats_hist(ev, f_name, raw=True, smoothed=True, conditional=True):
    green_ind = 1
    blue_ind = 0
    if not conditional:
        blue_ind = -1

    native_prop = ev.raw_data[f_name][green_ind]
    decoy_prop = ev.raw_data[f_name][blue_ind]

    u = ev.mean[f_name][-1]#np.mean(data)
    std = ev.std[f_name][-1]#(np.mean([(val-u)**2 for val in data]))**0.5
    if u == None: return

    plt.rcParams['figure.figsize'] = (10,7.5)
    plt.rcParams['figure.dpi'] = 400
    plt.rcParams['font.size'] = 20
    
    xmin, xmax = -0.5, 1.5# -0.18,1.18#min(data) - 0.5*std, max(data) + 0.5*std
    binwidth = 0.05#0.4*std
    bins = np.arange(xmin, xmax + binwidth, binwidth)
    
    if raw:
        xmin, xmax = 0, 1
        binwidth = 0.1#0.4*std
        bins = np.arange(xmin, xmax + binwidth, binwidth)

        plt.hist(decoy_prop, bins=bins, color='b', normed=True,
                 histtype='step', linewidth=2)#alpha=0.5)
        plt.hist(decoy_prop, bins=bins, color='b', normed=True,
                 histtype='stepfilled', linewidth=0, alpha=0.2)
        plt.hist(native_prop, bins=bins, color='g', normed=True,
                 histtype='step', linewidth=2)#,alpha=0.8)
        plt.hist(native_prop, bins=bins, color='g', normed=True,
                 histtype='stepfilled', linewidth=0, alpha=0.2)

        ax = plt.subplot(111)  
        ax.spines["top"].set_visible(False)  
        ax.spines["right"].set_visible(False)
        ax.get_xaxis().tick_bottom()  
        ax.get_yaxis().tick_left()

        #plt.xlim(u-4*std,u+4*std)
        plt.xlim(xmin, xmax)
        #xmin, xmax = plt.xlim()

        plt.show()
    
    if smoothed:
        xvals = np.arange(xmin - std,xmax + std,std*0.01)#0.02)

        norm = 1
        #if not conditional:
        #    norm = len(native_prop)/float(len(decoy_prop))

        y_native = [norm*ev.evaluate(f_name, x, green_ind) for x in xvals]
        y_decoy = [ev.evaluate(f_name, x, blue_ind) for x in xvals]

        plt.plot(xvals, y_decoy, linewidth=2, color='b')
        plt.plot(xvals, y_native, linewidth=2, color='g')

        #if not conditional:
        #    plt.fill_between(xvals, 0, y_native, facecolor='b',alpha=0.02)
        #    plt.fill_between(xvals, y_native, y_decoy, facecolor='b',alpha=0.2)
        #else:
        plt.fill_between(xvals, 0, y_decoy, facecolor='b',alpha=0.2)
        plt.fill_between(xvals, 0, y_native, facecolor='g',alpha=0.2)

        plt.xlim(xmin, xmax)

        ax = plt.subplot(111)  
        ax.spines["top"].set_visible(False)  
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.get_xaxis().tick_bottom()  
        ax.get_yaxis().set_visible(False)#.tick_left()

        #plt.gca().axes.yaxis.set_ticklabels([])
        plt.xticks([0,1])
        plt.tick_params(axis='x', direction='inout', length=13, width=3, pad=15)
        ax.axhline(linewidth=4, color='k')#, color="g")

        plt.show()

    return
    lines = [':','--']
    for i,p in enumerate(prot):
        #native_p = [(x.props[prop]-mean)/std for x in all_pairs[p] if x.correct()]
        #decoy_p = [(x.props[prop]-mean)/std for x in all_pairs[p] if not x.correct()]
        #y11 = yvals(native_p, xvals)
        #y21 = yvals(decoy_p, xvals)
        plt.plot(xvals, y11, 
                 lines[i%len(lines)],color='g',label=p,linewidth=3)
        plt.plot(xvals, y21, 
                 lines[i%len(lines)],color='b',label=p,linewidth=3)
        plt.fill_between(xvals, 0, y11, facecolor='g',alpha=0.1)#, linewidth=2, color='b')
        plt.fill_between(xvals, 0, y21, facecolor='b',alpha=0.1)#, linewidth=2, color='g')
        
    ax = plt.subplot(111)  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left()
        
    #plt.legend()
    plt.show()
    #break






