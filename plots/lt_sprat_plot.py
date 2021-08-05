import itertools
import os

from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as patches
import matplotlib.ticker as plticker
import matplotlib as mpl
import numpy as np
from scipy.signal import find_peaks

from rascal import util

# Load the LT SPRAT data
if '__file__' in locals():
    base_dir = os.path.dirname(__file__)
else:
    base_dir = os.getcwd()

colors = list(plt.cm.tab10(np.arange(10))) + ["black", "firebrick"]
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)


def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])


max_tries = [100, 150, 200, 250, 500, 1000, 2000, 5000]

# Number of repetance
N = 1000

best_p_mt = np.load(os.path.join(base_dir, 'sprat_best_p_auto_mt.npy'))
rms_mt = np.load(os.path.join(base_dir, 'sprat_rms_auto_mt.npy'))
residual_mt = np.load(os.path.join(base_dir, 'sprat_residual_auto_mt.npy'),
                      allow_pickle=True)
peak_utilisation_mt = np.load(os.path.join(
    base_dir, 'sprat_peak_utilisation_auto_mt.npy'),
                              allow_pickle=True)
atlas_utilisation_mt = np.load(os.path.join(
    base_dir, 'sprat_atlas_utilisation_auto_mt.npy'),
                               allow_pickle=True)

best_p_manual_mt = np.load(os.path.join(base_dir,
                                        'sprat_best_p_manual_mt.npy'))
rms_manual_mt = np.load(os.path.join(base_dir, 'sprat_rms_manual_mt.npy'))
residual_manual_mt = np.load(os.path.join(base_dir,
                                          'sprat_residual_manual_mt.npy'),
                             allow_pickle=True)
peak_utilisation_manual_mt = np.load(os.path.join(
    base_dir, 'sprat_peak_utilisation_manual_mt.npy'),
                                     allow_pickle=True)
atlas_utilisation_manual_mt = np.load(os.path.join(
    base_dir, 'sprat_atlas_utilisation_manual_mt.npy'),
                                      allow_pickle=True)
'''
# Figure 2 - polynomial coefficients
p0_range = np.nanpercentile(best_p_mt[-1][:, 0], [1., 99.])
p1_range = np.nanpercentile(best_p_mt[-1][:, 1], [1., 99.])
p2_range = np.nanpercentile(best_p_mt[-1][:, 2], [1., 99.])
p3_range = np.nanpercentile(best_p_mt[-1][:, 3], [1., 99.])
p4_range = np.nanpercentile(best_p_mt[-1][:, 4], [1., 99.])
p0_range_manual = np.nanpercentile(best_p_manual_mt[-1][:, 0], [1., 99.])
p1_range_manual = np.nanpercentile(best_p_manual_mt[-1][:, 1], [1., 99.])
p2_range_manual = np.nanpercentile(best_p_manual_mt[-1][:, 2], [1., 99.])
p3_range_manual = np.nanpercentile(best_p_manual_mt[-1][:, 3], [1., 99.])
p4_range_manual = np.nanpercentile(best_p_manual_mt[-1][:, 4], [1., 99.])

p0_kurtosis = np.zeros(len(max_tries))
p1_kurtosis = np.zeros(len(max_tries))
p2_kurtosis = np.zeros(len(max_tries))
p3_kurtosis = np.zeros(len(max_tries))
p4_kurtosis = np.zeros(len(max_tries))
p0_manual_kurtosis = np.zeros(len(max_tries))
p1_manual_kurtosis = np.zeros(len(max_tries))
p2_manual_kurtosis = np.zeros(len(max_tries))
p3_manual_kurtosis = np.zeros(len(max_tries))
p4_manual_kurtosis = np.zeros(len(max_tries))

fig2, ax2 = plt.subplots(2, 5, sharey=True)
fig2.set_figheight(10)
fig2.set_figwidth(10)
for i, mt in enumerate(max_tries):
    # First row - auto lines
    ax2[0, 0].hist(best_p_mt[i][:, 0],
                   bins=50,
                   range=p0_range,
                   histtype='step',
                   label=str(mt))
    p0_kurtosis[i] = kurtosis(best_p_mt[i][:, 0], bias=False)
    ax2[0, 1].hist(best_p_mt[i][:, 1],
                   bins=50,
                   range=p1_range,
                   histtype='step')
    p1_kurtosis[i] = kurtosis(best_p_mt[i][:, 1], bias=False)
    ax2[0, 2].hist(best_p_mt[i][:, 2],
                   bins=50,
                   range=p2_range,
                   histtype='step')
    p2_kurtosis[i] = kurtosis(best_p_mt[i][:, 2], bias=False)
    ax2[0, 3].hist(best_p_mt[i][:, 3],
                   bins=50,
                   range=p3_range,
                   histtype='step')
    p3_kurtosis[i] = kurtosis(best_p_mt[i][:, 3], bias=False)
    ax2[0, 4].hist(best_p_mt[i][:, 4],
                   bins=50,
                   range=p4_range,
                   histtype='step')
    p4_kurtosis[i] = kurtosis(best_p_mt[i][:, 4], bias=False)
    # Second row - manual lines
    ax2[1, 0].hist(best_p_manual_mt[i][:, 0],
                   bins=50,
                   range=p0_range_manual,
                   histtype='step')
    p0_manual_kurtosis[i] = kurtosis(best_p_manual_mt[i][:, 0])
    ax2[1, 1].hist(best_p_manual_mt[i][:, 1],
                   bins=50,
                   range=p1_range_manual,
                   histtype='step')
    p1_manual_kurtosis[i] = kurtosis(best_p_manual_mt[i][:, 1])
    ax2[1, 2].hist(best_p_manual_mt[i][:, 2],
                   bins=50,
                   range=p2_range_manual,
                   histtype='step')
    p2_manual_kurtosis[i] = kurtosis(best_p_manual_mt[i][:, 2])
    ax2[1, 3].hist(best_p_manual_mt[i][:, 3],
                   bins=50,
                   range=p3_range_manual,
                   histtype='step')
    p3_manual_kurtosis[i] = kurtosis(best_p_manual_mt[i][:, 3])
    ax2[1, 4].hist(best_p_manual_mt[i][:, 4],
                   bins=50,
                   range=p4_range_manual,
                   histtype='step')
    p4_manual_kurtosis[i] = kurtosis(best_p_manual_mt[i][:, 4])

ax2[0, 0].grid()
ax2[0, 1].grid()
ax2[0, 2].grid()
ax2[0, 3].grid()
ax2[0, 4].grid()

ax2[1, 0].grid()
ax2[1, 1].grid()
ax2[1, 2].grid()
ax2[1, 3].grid()
ax2[1, 4].grid()

ax2[0, 0].set_title('coeff 0')
ax2[0, 1].set_title('coeff 1')
ax2[0, 2].set_title('coeff 2')
ax2[0, 3].set_title('coeff 3')
ax2[0, 4].set_title('coeff 4')

ax2[0, 0].set_ylabel('Auto lines')
ax2[1, 0].set_ylabel('Manual lines')

handles, labels = ax2[0, 0].get_legend_handles_labels()

fig2.legend(flip(handles, 6),
            flip(labels, 6),
            loc='lower center',
            mode='expand',
            ncol=6)
fig2.tight_layout(w_pad=0.01)
fig2.subplots_adjust(bottom=0.1)

fig2.savefig('figure_2_polynomial_coefficients.png')
'''

# Figure 3 - RMS
fig3, ax3 = plt.subplots(2, 1, sharex=True, sharey=True)
fig3.set_figheight(6)
fig3.set_figwidth(6)

ax3[0].violinplot(rms_mt.T, showmedians=True)
ax3[1].violinplot(rms_manual_mt.T, showmedians=True)

ax3[0].set_xticks(range(1, len(max_tries) + 1))
ax3[0].set_xticklabels(max_tries)
ax3[1].set_xticks(range(1, len(max_tries) + 1))
ax3[1].set_xticklabels(max_tries)

ax3[0].grid()
ax3[1].grid()

ax3[0].set_ylabel(r'RMS / $\AA$')
ax3[1].set_xlabel('max_tries')
ax3[1].set_ylabel(r'RMS / $\AA$')

fig3.tight_layout()
fig3.savefig('figure_2_rms.png')

# Figure 4 - Peak Utilisation
fig4, ax4 = plt.subplots(2, 1, sharex=True, sharey=True)
fig4.set_figheight(6)
fig4.set_figwidth(6)

labels = []


def add_label(violin, label):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))


add_label(ax4[0].violinplot(peak_utilisation_mt.T * 100., showmedians=True),
          label='Peak Utilisation')
add_label(ax4[0].violinplot(atlas_utilisation_mt.T * 100., showmedians=True),
          label='Atlas Utilisation')
ax4[1].violinplot(peak_utilisation_manual_mt.T * 100., showmedians=True)
ax4[1].violinplot(atlas_utilisation_manual_mt.T * 100., showmedians=True)

ax4[0].set_xticks(range(1, len(max_tries) + 1))
ax4[0].set_xticklabels(max_tries)
ax4[1].set_xticks(range(1, len(max_tries) + 1))
ax4[1].set_xticklabels(max_tries)

ax4[0].grid()
ax4[1].grid()

ax4[1].set_xlabel('max_tries')
ax4[0].set_ylabel('Percentage')
ax4[1].set_ylabel('Percentage')

ax4[0].legend(*zip(*labels))

fig4.tight_layout()
fig4.savefig('figure_3_peak_atlas_utilisation.png')

# Figure 5 - wavelength at chosen pixels
# Note: start counting from ZERO
pix = [150., 350., 550., 750., 950.]
wave = np.array([np.array([np.zeros(len(pix))] * N)] * len(max_tries))
wave_manual = np.array([np.array([np.zeros(len(pix))] * N)] * len(max_tries))
for i, mt in enumerate(max_tries):
    for j in range(N):
        best_p = best_p_mt[i][j]
        best_p_manual = best_p_manual_mt[i][j]
        wave[i][j] = np.polynomial.polynomial.polyval(pix, best_p)
        wave_manual[i][j] = np.polynomial.polynomial.polyval(
            pix, best_p_manual)

resolution1 = np.nanmean(wave_manual[-1][:, 0] / 300.)
resolution2 = np.nanmean(wave_manual[-1][:, 1] / 300.)
resolution3 = np.nanmean(wave_manual[-1][:, 2] / 300.)
resolution4 = np.nanmean(wave_manual[-1][:, 3] / 300.)
resolution5 = np.nanmean(wave_manual[-1][:, 4] / 300.)

pix0_range = (np.nanmean(wave_manual[-1][:, 0]) - resolution1 * 50,
              np.nanmean(wave_manual[-1][:, 0]) + resolution1 * 50)
pix1_range = (np.nanmean(wave_manual[-1][:, 1]) - resolution1 * 50,
              np.nanmean(wave_manual[-1][:, 1]) + resolution1 * 50)
pix2_range = (np.nanmean(wave_manual[-1][:, 2]) - resolution1 * 50,
              np.nanmean(wave_manual[-1][:, 2]) + resolution1 * 50)
pix3_range = (np.nanmean(wave_manual[-1][:, 3]) - resolution1 * 50,
              np.nanmean(wave_manual[-1][:, 3]) + resolution1 * 50)
pix4_range = (np.nanmean(wave_manual[-1][:, 4]) - resolution1 * 50,
              np.nanmean(wave_manual[-1][:, 4]) + resolution1 * 50)

pix0_range_manual = (np.nanmean(wave_manual[-1][:, 0]) - resolution1 * 5,
                     np.nanmean(wave_manual[-1][:, 0]) + resolution1 * 5)
pix1_range_manual = (np.nanmean(wave_manual[-1][:, 1]) - resolution1 * 5,
                     np.nanmean(wave_manual[-1][:, 1]) + resolution1 * 5)
pix2_range_manual = (np.nanmean(wave_manual[-1][:, 2]) - resolution1 * 5,
                     np.nanmean(wave_manual[-1][:, 2]) + resolution1 * 5)
pix3_range_manual = (np.nanmean(wave_manual[-1][:, 3]) - resolution1 * 5,
                     np.nanmean(wave_manual[-1][:, 3]) + resolution1 * 5)
pix4_range_manual = (np.nanmean(wave_manual[-1][:, 4]) - resolution1 * 5,
                     np.nanmean(wave_manual[-1][:, 4]) + resolution1 * 5)

fig5, ax5 = plt.subplots(2, 5, sharey=True)
fig5.set_figheight(10)
fig5.set_figwidth(10)
for i, mt in enumerate(max_tries):
    # First row - auto lines
    ax5[0, 0].hist(wave[i][:, 0],
                   bins=20,
                   range=pix0_range,
                   histtype='step',
                   label=str(mt))
    ax5[0, 1].hist(wave[i][:, 1], bins=20, range=pix1_range, histtype='step')
    ax5[0, 2].hist(wave[i][:, 2], bins=20, range=pix2_range, histtype='step')
    ax5[0, 3].hist(wave[i][:, 3], bins=20, range=pix3_range, histtype='step')
    ax5[0, 4].hist(wave[i][:, 4], bins=20, range=pix4_range, histtype='step')
    # Second row - manual lines
    ax5[1, 0].hist(wave_manual[i][:, 0],
                   bins=20,
                   range=pix0_range_manual,
                   lw=2,
                   histtype='step')
    ax5[1, 1].hist(wave_manual[i][:, 1],
                   bins=20,
                   range=pix1_range_manual,
                   lw=2,
                   histtype='step')
    ax5[1, 2].hist(wave_manual[i][:, 2],
                   bins=20,
                   range=pix2_range_manual,
                   lw=2,
                   histtype='step')
    ax5[1, 3].hist(wave_manual[i][:, 3],
                   bins=20,
                   range=pix3_range_manual,
                   lw=2,
                   histtype='step')
    ax5[1, 4].hist(wave_manual[i][:, 4],
                   bins=20,
                   range=pix4_range_manual,
                   lw=2,
                   histtype='step')

ax5[0, 0].grid()
ax5[0, 1].grid()
ax5[0, 2].grid()
ax5[0, 3].grid()
ax5[0, 4].grid()

ax5[1, 0].grid()
ax5[1, 1].grid()
ax5[1, 2].grid()
ax5[1, 3].grid()
ax5[1, 4].grid()

ax5[0, 0].set_title('Pix 150')
ax5[0, 1].set_title('Pix 350')
ax5[0, 2].set_title('Pix 550')
ax5[0, 3].set_title('Pix 750')
ax5[0, 4].set_title('Pix 950')

ax5[0, 0].set_ylim(0, 1000)
ax5[0, 1].set_ylim(0, 1000)
ax5[0, 2].set_ylim(0, 1000)
ax5[0, 3].set_ylim(0, 1000)
ax5[0, 4].set_ylim(0, 1000)

ax5[1, 0].set_ylim(0, 1000)
ax5[1, 1].set_ylim(0, 1000)
ax5[1, 2].set_ylim(0, 1000)
ax5[1, 3].set_ylim(0, 1000)
ax5[1, 4].set_ylim(0, 1000)

ax5[0, 0].set_xlim(
    np.nanmean(wave_manual[-1][:, 0]) - resolution1 * 50,
    np.nanmean(wave_manual[-1][:, 0]) + resolution1 * 50)
ax5[0, 1].set_xlim(
    np.nanmean(wave_manual[-1][:, 1]) - resolution2 * 50,
    np.nanmean(wave_manual[-1][:, 1]) + resolution2 * 50)
ax5[0, 2].set_xlim(
    np.nanmean(wave_manual[-1][:, 2]) - resolution3 * 50,
    np.nanmean(wave_manual[-1][:, 2]) + resolution3 * 50)
ax5[0, 3].set_xlim(
    np.nanmean(wave_manual[-1][:, 3]) - resolution4 * 50,
    np.nanmean(wave_manual[-1][:, 3]) + resolution4 * 50)
ax5[0, 4].set_xlim(
    np.nanmean(wave_manual[-1][:, 4]) - resolution5 * 50,
    np.nanmean(wave_manual[-1][:, 4]) + resolution5 * 50)

loc1 = plticker.MultipleLocator(base=resolution1 * 10)
loc2 = plticker.MultipleLocator(base=resolution2 * 10)
loc3 = plticker.MultipleLocator(base=resolution3 * 10)
loc4 = plticker.MultipleLocator(base=resolution4 * 10)
loc5 = plticker.MultipleLocator(base=resolution5 * 10)

ax5[1, 0].set_xlim(
    np.nanmean(wave_manual[-1][:, 0]) - resolution1 * 5,
    np.nanmean(wave_manual[-1][:, 0]) + resolution1 * 5)
ax5[1, 1].set_xlim(
    np.nanmean(wave_manual[-1][:, 1]) - resolution2 * 5,
    np.nanmean(wave_manual[-1][:, 1]) + resolution2 * 5)
ax5[1, 2].set_xlim(
    np.nanmean(wave_manual[-1][:, 2]) - resolution3 * 5,
    np.nanmean(wave_manual[-1][:, 2]) + resolution3 * 5)
ax5[1, 3].set_xlim(
    np.nanmean(wave_manual[-1][:, 3]) - resolution4 * 5,
    np.nanmean(wave_manual[-1][:, 3]) + resolution4 * 5)
ax5[1, 4].set_xlim(
    np.nanmean(wave_manual[-1][:, 4]) - resolution5 * 5,
    np.nanmean(wave_manual[-1][:, 4]) + resolution5 * 5)

loc1_manual = plticker.MultipleLocator(base=resolution1)
loc2_manual = plticker.MultipleLocator(base=resolution2)
loc3_manual = plticker.MultipleLocator(base=resolution3)
loc4_manual = plticker.MultipleLocator(base=resolution4)
loc5_manual = plticker.MultipleLocator(base=resolution5)

ax5[0, 0].xaxis.set_major_locator(loc1)
ax5[0, 1].xaxis.set_major_locator(loc2)
ax5[0, 2].xaxis.set_major_locator(loc3)
ax5[0, 3].xaxis.set_major_locator(loc4)
ax5[0, 4].xaxis.set_major_locator(loc5)

ax5[1, 0].xaxis.set_major_locator(loc1_manual)
ax5[1, 1].xaxis.set_major_locator(loc2_manual)
ax5[1, 2].xaxis.set_major_locator(loc3_manual)
ax5[1, 3].xaxis.set_major_locator(loc4_manual)
ax5[1, 4].xaxis.set_major_locator(loc5_manual)

ax5[1, 2].set_xlabel('Wavelength / A')
ax5[0, 0].set_ylabel('Frequency')
ax5[1, 0].set_ylabel('Frequency')

handles, labels = ax5[0, 0].get_legend_handles_labels()

fig5.legend(flip(handles, 4),
            flip(labels, 4),
            loc='lower center',
            mode='expand',
            ncol=4)

fig5.tight_layout(w_pad=0.01)
fig5.subplots_adjust(bottom=0.14, wspace=0.15, hspace=0.15)

label1 = ax5[0, 0].get_xticklabels()
label2 = ax5[0, 1].get_xticklabels()
label3 = ax5[0, 2].get_xticklabels()
label4 = ax5[0, 3].get_xticklabels()
label5 = ax5[0, 4].get_xticklabels()

label1_manual = ax5[1, 0].get_xticklabels()
label2_manual = ax5[1, 1].get_xticklabels()
label3_manual = ax5[1, 2].get_xticklabels()
label4_manual = ax5[1, 3].get_xticklabels()
label5_manual = ax5[1, 4].get_xticklabels()

ax5[0, 0].set_xticklabels(label1, rotation=270)
ax5[0, 1].set_xticklabels(label2, rotation=270)
ax5[0, 2].set_xticklabels(label3, rotation=270)
ax5[0, 3].set_xticklabels(label4, rotation=270)
ax5[0, 4].set_xticklabels(label5, rotation=270)

ax5[1, 0].set_xticklabels(label1_manual, rotation=270)
ax5[1, 1].set_xticklabels(label2_manual, rotation=270)
ax5[1, 2].set_xticklabels(label3_manual, rotation=270)
ax5[1, 3].set_xticklabels(label4_manual, rotation=270)
ax5[1, 4].set_xticklabels(label5_manual, rotation=270)

rect1 = patches.Rectangle(
    (np.nanmean(wave_manual[-1][:, 0]) - resolution1 * 10, -10),
    resolution1 * 20,
    1100,
    lw=2,
    edgecolor='grey',
    facecolor='lightcyan')
rect2 = patches.Rectangle(
    (np.nanmean(wave_manual[-1][:, 1]) - resolution2 * 10, -10),
    resolution2 * 20,
    1100,
    lw=2,
    edgecolor='grey',
    facecolor='lightcyan')
rect3 = patches.Rectangle(
    (np.nanmean(wave_manual[-1][:, 2]) - resolution3 * 10, -10),
    resolution3 * 20,
    1100,
    lw=2,
    edgecolor='grey',
    facecolor='lightcyan')
rect4 = patches.Rectangle(
    (np.nanmean(wave_manual[-1][:, 3]) - resolution4 * 10, -10),
    resolution4 * 20,
    1100,
    lw=2,
    edgecolor='grey',
    facecolor='lightcyan')
rect5 = patches.Rectangle(
    (np.nanmean(wave_manual[-1][:, 4]) - resolution5 * 10, -10),
    resolution5 * 20,
    1100,
    lw=2,
    edgecolor='grey',
    facecolor='lightcyan')

ax5[0, 0].set_facecolor('cornsilk')
ax5[0, 1].set_facecolor('cornsilk')
ax5[0, 2].set_facecolor('cornsilk')
ax5[0, 3].set_facecolor('cornsilk')
ax5[0, 4].set_facecolor('cornsilk')

ax5[0, 0].add_patch(rect1)
ax5[0, 1].add_patch(rect2)
ax5[0, 2].add_patch(rect3)
ax5[0, 3].add_patch(rect4)
ax5[0, 4].add_patch(rect5)

ax5[1, 0].set_facecolor('lightcyan')
ax5[1, 1].set_facecolor('lightcyan')
ax5[1, 2].set_facecolor('lightcyan')
ax5[1, 3].set_facecolor('lightcyan')
ax5[1, 4].set_facecolor('lightcyan')

fig5.savefig('figure_4_wavelengths.png')

# Figure 6 - 2D heatmap of the solution


fig6, ax6 = plt.subplots(2, 1, sharex=True, sharey=False)
fig6.set_figheight(6)
fig6.set_figwidth(6)

pix = np.arange(1024) + 1
wave = np.zeros((N, len(pix)))
wave_manual = np.zeros((N, len(pix)))

i = -1
for j in range(N):
    best_p = best_p_mt[i][j]
    best_p_manual = best_p_manual_mt[i][j]
    wave[j] = np.polynomial.polynomial.polyval(pix, best_p)
    wave_manual[j] = np.polynomial.polynomial.polyval(pix, best_p_manual)

delta_wave = wave - np.nanmedian(wave, axis=0)
delta_wave_manual = wave_manual - np.nanmedian(wave_manual, axis=0)

dw_min, dw_max = -500., 500.
dw_manual_min, dw_manual_max = -50., 50.

delta_wave_heatmap = []
delta_wave_manual_heatmap = []

for i in range(N):
    delta_wave_heatmap.append(
        np.histogram(delta_wave[:, i], bins=100, range=(dw_min, dw_max))[0])
    delta_wave_manual_heatmap.append(
        np.histogram(delta_wave_manual[:, i],
                     bins=100,
                     range=(dw_manual_min, dw_manual_max))[0])

dw_yedges = np.histogram(delta_wave, bins=100, range=(dw_min, dw_max))[1]
dw_manual_yedges = np.histogram(delta_wave_manual,
                                bins=100,
                                range=(dw_manual_min, dw_manual_max))[1]

ax6[0].imshow(np.array(delta_wave_heatmap).T, origin='lower', aspect='auto')
ax6[0].set_yticks(np.linspace(0, 100, len(dw_yedges[::10])) - 0.5)
ax6[0].set_yticklabels(dw_yedges[::10].astype('int'))
ax6[0].set_ylabel(r'$\Delta\lambda\ /\ \AA$')

ax6[1].imshow(np.array(delta_wave_manual_heatmap).T,
              origin='lower',
              aspect='auto')
ax6[1].set_yticks(np.linspace(0, 100, len(dw_manual_yedges[::10])) - 0.5)
ax6[1].set_yticklabels(dw_manual_yedges[::10].astype('int'))

ax6[1].set_ylabel(r'$\Delta\lambda\ /\ \AA$')
ax6[1].set_xticks(pix[::200] - 1)
ax6[1].set_xticklabels((pix[::200] - 1).astype('int'))
ax6[1].set_xlabel('Calibration Run #')


ax6[0].vlines(
    [173,
     1000],
    ymin=0,
    ymax=99,
    color='grey',
    ls=':',
    lw=1)
ax6[1].vlines(
    [173,
     1000],
    ymin=0,
    ymax=99,
    color='grey',
    ls=':',
    lw=1)
res = wave[-1] / 300.

# auto calibration
ax6[0].plot(50 + res, color='black', lw=1, ls='dashed')
ax6[0].plot(50 + res * 2, color='black', lw=1, ls='dashed')
ax6[0].plot(50 + res * 3, color='black', lw=1, ls='dashed')
ax6[0].plot(50 + res * 4, color='black', lw=1, ls='dashed')
ax6[0].plot(50 + res * 5, color='black', lw=1, ls='dashed')
ax6[0].plot([0, 1000], [50, 50], color='black', lw=1, ls='dashed')
ax6[0].plot(50 - res, color='black', lw=1, ls='dashed')
ax6[0].plot(50 - res * 2, color='black', lw=1, ls='dashed')
ax6[0].plot(50 - res * 3, color='black', lw=1, ls='dashed')
ax6[0].plot(50 - res * 4, color='black', lw=1, ls='dashed')
ax6[0].plot(50 - res * 5, color='black', lw=1, ls='dashed')
ax6[0].plot([0, 1000], [55, 55], color='white', lw=1, ls=':')
ax6[0].plot([0, 1000], [45, 45], color='white', lw=1, ls=':')
ax6[0].set_xlim(0, 1000)
ax6[0].set_ylim(0, 99.5)

# manual calibration
ax6[1].plot(50 + res, color='black', lw=1, ls=':')
ax6[1].plot(50 + res * 2, color='black', lw=1, ls=':')
ax6[1].plot(50 + res * 3, color='black', lw=1, ls=':')
ax6[1].plot(50 + res * 4, color='black', lw=1, ls=':')
ax6[1].plot(50 + res * 5, color='black', lw=1, ls=':')
ax6[1].plot([0, 1000], [50, 50], color='black', lw=1, ls=':')
ax6[1].plot(50 - res, color='black', lw=1, ls=':')
ax6[1].plot(50 - res * 2, color='black', lw=1, ls=':')
ax6[1].plot(50 - res * 3, color='black', lw=1, ls=':')
ax6[1].plot(50 - res * 4, color='black', lw=1, ls=':')
ax6[1].plot(50 - res * 5, color='black', lw=1, ls=':')
ax6[1].set_xlim(0, 1000)
ax6[1].set_ylim(-0.5, 100)

fig6.tight_layout()
fig6.subplots_adjust(hspace=0)
fig6.savefig('figure_5_heatmap.png')

print('Auto calibration:')
for i, mt in enumerate(max_tries):
    print('max_tries = {}.'.format(i))
    print('C0 = {} +/- {}.'.format(np.nanmean(best_p_mt[i][:, 0]),
                                   np.nanstd(best_p_mt[i][:, 0])))
    print('C1 = {} +/- {}.'.format(np.nanmean(best_p_mt[i][:, 1]),
                                   np.nanstd(best_p_mt[i][:, 1])))
    print('C2 = {} +/- {}.'.format(np.nanmean(best_p_mt[i][:, 2]),
                                   np.nanstd(best_p_mt[i][:, 2])))
    print('C3 = {} +/- {}.'.format(np.nanmean(best_p_mt[i][:, 3]),
                                   np.nanstd(best_p_mt[i][:, 3])))
    print('C4 = {} +/- {}.'.format(np.nanmean(best_p_mt[i][:, 4]),
                                   np.nanstd(best_p_mt[i][:, 4])))
    print('RMS = {} +/- {}.'.format(np.nanmean(rms_mt[i]),
                                    np.nanstd(rms_mt[i])))
    print('Peak Utilisation = {} +/- {}.'.format(
        np.nanmean(peak_utilisation_mt[i]), np.nanstd(peak_utilisation_mt[i])))
    print('Atlas Utilisation = {} +/- {}.'.format(
        np.nanmean(atlas_utilisation_mt[i]),
        np.nanstd(atlas_utilisation_mt[i])))
    wave150 = []
    wave350 = []
    wave550 = []
    wave750 = []
    wave950 = []
    for j in range(N):
        wave150.append(np.polynomial.polynomial.polyval(150, best_p_mt[i][j]))
        wave350.append(np.polynomial.polynomial.polyval(350, best_p_mt[i][j]))
        wave550.append(np.polynomial.polynomial.polyval(550, best_p_mt[i][j]))
        wave750.append(np.polynomial.polynomial.polyval(750, best_p_mt[i][j]))
        wave950.append(np.polynomial.polynomial.polyval(950, best_p_mt[i][j]))

    print('Pix 150 = {} +/- {}.'.format(np.nanmean(wave150),
                                        np.nanstd(wave150)))
    print('Pix 350 = {} +/- {}.'.format(np.nanmean(wave350),
                                        np.nanstd(wave350)))
    print('Pix 550 = {} +/- {}.'.format(np.nanmean(wave550),
                                        np.nanstd(wave550)))
    print('Pix 750 = {} +/- {}.'.format(np.nanmean(wave750),
                                        np.nanstd(wave750)))
    print('Pix 950 = {} +/- {}.'.format(np.nanmean(wave950),
                                        np.nanstd(wave950)))

print('Manual calibration:')
for i, mt in enumerate(max_tries):
    print('max_tries = {}.'.format(i))
    print('C0 = {} +/- {}.'.format(np.nanmean(best_p_manual_mt[i][:, 0]),
                                   np.nanstd(best_p_manual_mt[i][:, 0])))
    print('C1 = {} +/- {}.'.format(np.nanmean(best_p_manual_mt[i][:, 1]),
                                   np.nanstd(best_p_manual_mt[i][:, 1])))
    print('C2 = {} +/- {}.'.format(np.nanmean(best_p_manual_mt[i][:, 2]),
                                   np.nanstd(best_p_manual_mt[i][:, 2])))
    print('C3 = {} +/- {}.'.format(np.nanmean(best_p_manual_mt[i][:, 3]),
                                   np.nanstd(best_p_manual_mt[i][:, 3])))
    print('C4 = {} +/- {}.'.format(np.nanmean(best_p_manual_mt[i][:, 4]),
                                   np.nanstd(best_p_manual_mt[i][:, 4])))
    print('RMS = {} +/- {}.'.format(np.nanmean(rms_manual_mt[i]),
                                    np.nanstd(rms_manual_mt[i])))
    print('Peak Utilisation = {} +/- {}.'.format(
        np.nanmean(peak_utilisation_manual_mt[i]),
        np.nanstd(peak_utilisation_manual_mt[i])))
    print('Atlas Utilisation = {} +/- {}.'.format(
        np.nanmean(atlas_utilisation_manual_mt[i]),
        np.nanstd(atlas_utilisation_manual_mt[i])))
    wave150 = []
    wave350 = []
    wave550 = []
    wave750 = []
    wave950 = []
    for j in range(N):
        wave150.append(
            np.polynomial.polynomial.polyval(150, best_p_manual_mt[i][j]))
        wave350.append(
            np.polynomial.polynomial.polyval(350, best_p_manual_mt[i][j]))
        wave550.append(
            np.polynomial.polynomial.polyval(550, best_p_manual_mt[i][j]))
        wave750.append(
            np.polynomial.polynomial.polyval(750, best_p_manual_mt[i][j]))
        wave950.append(
            np.polynomial.polynomial.polyval(950, best_p_manual_mt[i][j]))

    print('Pix 150 = {} +/- {}.'.format(np.nanmean(wave150),
                                        np.nanstd(wave150)))
    print('Pix 350 = {} +/- {}.'.format(np.nanmean(wave350),
                                        np.nanstd(wave350)))
    print('Pix 550 = {} +/- {}.'.format(np.nanmean(wave550),
                                        np.nanstd(wave550)))
    print('Pix 750 = {} +/- {}.'.format(np.nanmean(wave750),
                                        np.nanstd(wave750)))
    print('Pix 950 = {} +/- {}.'.format(np.nanmean(wave950),
                                        np.nanstd(wave950)))
