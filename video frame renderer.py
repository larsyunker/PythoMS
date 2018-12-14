"""
Extracts data from an mzML file and uses the matplotlib animation class to create an animation of the mass spectrum
and the abundance traces for a highlighted region.
The left subplot is a mass spectrum at the current time
The right subplot is a chart of intensity traces leading up to the current time

Some modification of the script behaviour is supported with the variable assignments at the top of the script.
The video is coded to save using ffmpeg (https://www.ffmpeg.org/). Installation instructions for Windows can be found
here: https://www.wikihow.com/Install-FFmpeg-on-Windows . Simplified modification of environment modification can be
accomplished using Rapid Environment Editor (https://www.rapidee.com/).


"""
import matplotlib.animation as animation
from PyRSIR import pyrsir
from pythoms.tome import bindata, binnspectra
from pythoms.colour import Colour
from bisect import bisect_left, bisect_right
import pylab as pl
import sys
import os

# input *.raw filename
filename = os.path.join(os.getcwd(), 'validation_files', 'LY-2015-09-15 06.mzML.gz')

# species to track and plot
# sub and superscripts can be denoted by TeX formatting
# e.g. an Ar group with a positive charge is denoted by 'Ar$^+$', and CH2 would be denoted CH$_2$
# (see http://matplotlib.org/users/mathtext.html#subscripts-and-superscripts for more info)
# m/z bounds, the desired colour, and affinity must be specified for every species
sp = {
    'Ar$^+$I': {'bounds': [478.865, 481.622], 'colour': '#2078b4', 'affin': '+'},
    'Ar$^+$Ar': {'bounds': [443.016, 446.711], 'colour': '#34a048', 'affin': '+'},
}

# set number of scans to sum
n = 2

# define scan range
# (scan range that you wish to render)
scr = [1, 2147]
# scr = [261,427]

# define mz range for spectrum image
# (these will be the bounds of the x-axis)
mz = [440., 489.]

# left-right scalar for scan info placement
# 0 is fuly right, 1 is fully left
infop = 0.1

# reaction start point (eg. catalyst injection)
# provide time in minutes
inj = 4.612

# provide timepoints for injections/additions
# give the form 'injection name':time in min, for each timepoint (use actual time, not shifted time)
timepoints = {
    'catalyst injection': 4.612
}

# axis line width
axwidth = 1.5

# trace line width
lw = 1.5

# save an image every n scans 
save = 2

# scans per second
sps = 30

# font size
fs = 16


#####################################################################
def msfignorm(x, y):
    """
    Normalizes the height of a mass spectrum
    The height will be the sum of the heights of the base peaks in the window

    The function will normalize the y-values (assumes intensity) and return them
    """
    height = 0  # starting point
    for key in sp:
        # index location of selected peak in spectrum
        left = bisect_right(x, sp[key]['bounds'][0])
        right = bisect_left(x, sp[key]['bounds'][1])
        try:
            height += max(y[left:right])  # add maximum in selected region to height
        except ValueError:  # if no intensity in region, add some small number
            height += 0.01

    for ind, val in enumerate(y):  # normalizes all y values
        y[ind] = val / height

    return x, y


def timelimits(index):
    """
    finds the appropriate time limits for the traces
    """
    mintime = 10000
    maxtime = -10000
    for mode in mskeys:
        sumkey = f'{n}sum{mode}'
        if sumkey in rtime.keys():
            if rtime[sumkey][0] < mintime:
                mintime = rtime[sumkey][0]
            if index == 0:
                index += 1
            if rtime[sumkey][index] > maxtime:
                maxtime = rtime[sumkey][index]
    return mintime, maxtime


def animate(i):
    specplot.set_data(*spec[i])  # update mass spectrum
    for key in sp:  # update plots for each species
        curr = sp[key]
        curr['traceplot'].set_data(
            rtime[str(n) + 'sum' + curr['affin']][:i + 1],
            curr[str(n) + 'norm'][:i + 1]
        )
        # determine mass spectrum slice indicies
        left = bisect_left(
            spec[i][0],
            sp[key]['bounds'][0]
        )
        right = bisect_right(
            spec[i][0],
            sp[key]['bounds'][1]
        )
        curr['msplot'].set_data(  # colour mass spectrum appropriately
            spec[i][0][left:right],
            spec[i][1][left:right]
        )
    mintime, maxtime = timelimits(i)
    traceax.set_xlim((mintime, maxtime))  # update right trace x limits
    textx = maxtime - (maxtime - mintime) * infop  # calculate location for scan number and time text
    scantext.set_text('scan %d' % (i * n + 1))  # text for scan number
    scantext.set_position((textx, 0.96))
    timetext.set_text('%.1f min' % maxtime)  # text for time
    timetext.set_position((textx, 0.92))
    for ax in [specax, traceax]:
        for label in ax.get_yticklabels():
            label.set_fontproperties(tickfont)
        for label in ax.get_xticklabels():
            label.set_fontproperties(tickfont)
    return (specax, traceax)


if save < n:  # if the script is told to save more often than it sums
    save = n
mskeys = ['+', '-']
pyrsirkw = {
    'plot': False,  # plot the data for a quick look
    # 'verbose': True,  # chatty
    'bounds confidence': 0.99,  # confidence interval for automatically generated bounds
    'sumspec': False,  # whether or not to output a summed spectrum
    'return': True,  # whether to return data (if the data from the function is required by another function)
}
mzml, sp, rtime = pyrsir(filename, sp, 1, **pyrsirkw)[:3]  # run pyrsir

sstart = mzml.scan_index(scr[0])  # index of start scan
send = mzml.scan_index(scr[1])  # index of last scan
for key in sp:
    sp[key]['raw'] = sp[key]['raw'][sstart:send + 1]  # trim to scan range
for key in rtime:
    rtime[key] = rtime[key][sstart:send + 1]  # trim to scan range
spec = mzml.retrieve_scans(scr[0], scr[1], mz[0], mz[1], outside=True)  # pull all spectra within scan range
sys.stdout.write('%s summing and normalizing species traces' % str(n))
sumkey = str(n) + 'sum'
normkey = str(n) + 'norm'
sumsp = []
for key in sp:
    sp[key][sumkey] = bindata(n, sp[key]['raw'])  # bin each species
    sp[key]['colour'] = Colour(sp[key]['colour']).mpl  # convert colour into matplotlib format
    for ind, val in enumerate(sp[key][sumkey]):  # for normalization
        try:
            sumsp[ind] += val
        except IndexError:
            sumsp.append(val)

for mode in mskeys:
    sumkey = str(n) + 'sum' + mode
    modekey = 'raw' + mode
    if modekey in rtime.keys():  # if there is data for that mode
        rtime[sumkey] = bindata(n, rtime[modekey], n)
        for ind, val in enumerate(rtime[sumkey]):
            rtime[sumkey][ind] = val - inj  # shift time data to zero at injection point
        for key in sp:  # for each species
            if sp[key]['affin'] in mskeys:  # if species has affinity
                spkey = str(n) + 'sum'
                sp[key][normkey] = []
                for ind, val in enumerate(sp[key][spkey]):
                    sp[key][normkey].append(val / (sumsp[ind] + 0.01))  # +0.01 to avoid div/0 errors
sys.stdout.write(' DONE\n')
sys.stdout.flush()

# bin and normalize mass spectra
spec = [msfignorm(*spectrum) for spectrum in binnspectra(spec, n, start=mz[0], end=mz[1])]

# initial figure setup
fig, [specax, traceax] = pl.subplots(
    1,
    2,
    figsize=(19.2, 10.8),
    sharey=True
)

specplot, = specax.plot(  # plot for mass spectrum
    [],
    [],
    'k-',
    lw=1.5
)
specax.set_xlabel(
    'm/z',
    style='italic',
    fontname='Arial',
    fontsize=fs,
)
specax.set_ylabel(
    'Relative Intensity',
    fontname='Arial',
    fontsize=fs
)
specax.set_xlim(*mz)
specax.spines['bottom'].set_visible(False)

traceax.set_xlabel(
    'time (min)',
    fontname='Arial',
    fontsize=fs,
)

for ax in [specax, traceax]:
    for side in ['right', 'top']:  # set spines as invisible
        ax.spines[side].set_visible(False)
    for side in ['left', 'bottom']:
        ax.spines[side].set_linewidth(axwidth)
    ax.tick_params(
        labelsize=fs,
        length=axwidth * 3,
        width=axwidth,
        direction='out',
        right='off',
    )
    ax.set_ylim([
        -0.001,
        1.
    ])

for key in sp:  # create plot instances and labels
    sp[key]['traceplot'], = traceax.plot(
        [],
        [],
        linewidth=lw,
        label=key,
        color=sp[key]['colour'],
    )
    sp[key]['msplot'], = specax.plot(
        [],
        [],
        linewidth=1.5,
        color=sp[key]['colour'],
    )
    specax.text(
        sp[key]['bounds'][0],
        1.01,
        key,
        fontname='Arial',
        fontsize=fs,
        color=sp[key]['colour'],
    )

for key in timepoints:
    traceax.axvline(
        x=(timepoints[key] - inj),
        ymin=0,
        ymax=1,
        linewidth=lw,
        color='b',
        linestyle=':',
    )
    traceax.text(
        timepoints[key] - inj,
        0.5,
        key,
        fontsize=fs,
        color='b',
        backgroundcolor='w',
        rotation='vertical',
        horizontalalignment='center',
        verticalalignment='center',
        alpha=0.75,
        fontname='Arial',
    )

tickfont = pl.matplotlib.font_manager.FontProperties(family='Arial', size=fs)

scantext = pl.text(  # text for scan number
    0.5,
    0.96,
    '',
    fontname='Arial',
    fontsize=fs,
)
timetext = pl.text(
    0.5,
    0.92,
    '',
    fontname='Arial',
    fontsize=fs,
)

pl.subplots_adjust(  # hard coded subplot tightening
    left=0.07,
    right=0.99,
    bottom=0.095,
    top=0.96,
    wspace=0.06,
    hspace=0.05
)

ani = animation.FuncAnimation(
    fig,
    animate,
    frames=len(spec),
    interval=1 / sps * 1000,
    repeat=False,
)

print('Writing animation to file...')
ani.save(
    filename + '.mp4',
    writer='ffmpeg',
    bitrate=1000,
    codec='libx264',
)

