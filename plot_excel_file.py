"""
A script for plotting all spectra in an excel file as mass spectra and saving them
This allows for spectra without the same dimension to be compared

new:
    functional. creates simple figures with little customization
    ---1.0---
"""

def barplotthese(filename,startmz=50,endmz=700,colour='bl'):
    """
    plots all sheets in an excel file as bar plots
    """
    def msbarplot(x,y):
        """plot the spectrum"""
        fs = 16 # fontsize
        lw = 1.5 # line width
        pl.clf()
        pl.close()
        fig = pl.figure(figsize = (7.87,4.87))
        ax = fig.add_subplot(111)
        font = {'fontname':'Arial'} #font parameters for axis/text labels
        tickfont = pl.matplotlib.font_manager.FontProperties(family='Arial',size=fs) # font parameters for axis ticks
        ax.spines["bottom"].set_position(('axes',-0.01))
        pl.tick_params(axis='x', length=lw*3, width=lw, direction='out',top='off')
        pl.tick_params(axis='y', length=lw*3, width=lw, direction='out',right='off')
        for label in ax.get_yticklabels():
            label.set_fontproperties(tickfont)
        for label in ax.get_xticklabels():
            label.set_fontproperties(tickfont)
        ax.spines["right"].set_visible(False) #hide right axis
        ax.spines["top"].set_visible(False) # hide top axis
        ax.spines["bottom"].set_visible(False) # hide bottom axis (comment to toggle)
        ax.set_xlabel('m/z', style='italic', fontsize = fs,**font)
        ax.set_ylabel('normalized intensity', fontsize = fs,**font)
        for axis in ["top","bottom","left","right"]:
            ax.spines[axis].set_linewidth(lw)
        ax.set_xlim(startmz,endmz)
        ax.set_ylim(0.,100.)
        
        ax.bar(x,y,1.0,linewidth=0,color=Colour(colour).mpl) # plot bar
        
        pl.subplots_adjust(left = 0.10, right = 0.96, bottom = 0.12, top = 0.98) # manually adjusted for 16 fontsize
        pl.savefig(filename+' - '+sheet+' ('+`startmz`+'-'+`endmz`+')'+'.png',dpi=200,transparent=True)
    
    import os,sys
    sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/_classes')
    from _XLSX import XLSX
    from _Colour import Colour
    import pylab as pl
    from tome_v02 import normalize
    
    xlfile = XLSX(filename)
    
    for sheet in xlfile.wb.get_sheet_names(): # for each sheet
        sys.stdout.write('Plotting spectrum "%s"' %sheet)
        sys.stdout.flush()
        spec = xlfile.pullspectrum(sheet,8)[0] # pull spectrum
        spec[1] = normalize(spec[1],100.) # normalize to 100
        msbarplot(spec[0],spec[1]) # plot and save bar plot
        sys.stdout.write(' DONE\n')
    sys.stdout.write('fin.\n')
        
if __name__ == '__main__':
    startmz = 50
    endmz = 700
    colour = 'bl'
    filename = 'Copy of Book1'
    barplotthese(filename,startmz,endmz,colour) 