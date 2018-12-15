from matplotlib import rc, rcParams
from collections import OrderedDict

# ---------------------------------------------------
# Color sets
# ---------------------------------------------------
# Standard tableau 20 set
tableau = OrderedDict([
    ('blue',        '#1F77B4'),
    ('orange',      '#FF7F0E'),
    ('green',       '#2CA02C'),
    ('red',         '#D62728'),
    ('purple',      '#9467BD'),
    ('brown',       '#8C564B'),
    ('pink',        '#E377C2'),
    ('grey',        '#7F7F7F'),
    ('yellow',      '#BCBD22'),
    ('turquoise',   '#17BECF'),
    ('lightblue',   '#AEC7E8'),
    ('lightorange', '#FFBB78'),
    ('lightgreen',  '#98DF8A'),
    ('lightred',    '#FF9896'),
    ('lightpurple', '#C5B0D5'),
    ('lightbrown',  '#C49C94'),
    ('lightpink',   '#F7B6D2'),
    ('lightgrey',   '#C7C7C7'),
    ('lightyellow', '#DBDB8D'),
    ('lightturquoise', '#9EDAE5')
])

fontsize=18
nearly_black = '#161616'
light_grey = '#EEEEEE'
lighter_grey = '#F5F5F5'
white = '#ffffff'
grey = '#7F7F7F'

colours = { 'v_o' : tableau['blue'],
            'v_la' : tableau['yellow'],
            'v_zr' : tableau['green'],
            'zr_li_oct' : tableau['red'],
            'zr_li_tet' : tableau['grey'],
            'zr_la' : tableau['orange'],
            'la_zr' : tableau['turquoise'],
            'o_i' : tableau['brown'],
            'v_li' : tableau['purple'],
            'i_li' : tableau['pink'] }

formatting = { 'axes.formatter.limits': (-3,3),
               'xtick.major.pad': 7,
               'ytick.major.pad': 7,
               'ytick.color': nearly_black,
               'xtick.color': nearly_black,
               'axes.labelcolor': nearly_black,
               'legend.facecolor': light_grey,
               'pdf.fonttype': 42,
               'ps.fonttype': 42,
               'mathtext.fontset': 'custom',
               'font.size': fontsize,
               'text.usetex' : True,
               'font.family' : 'serif',
               'font.serif' : 'Palatino Linotype',
               'text.usetex': True,
               'savefig.bbox':'tight',
               'axes.facecolor': light_grey,
               'axes.labelpad': 10.0,
               'axes.labelsize': fontsize,
               'axes.titlepad': 30,
               'lines.markersize': 5.0,
               'lines.scale_dashes': False }

stab_region_formatting = {'axes.formatter.limits': (-3,3),
                          'xtick.major.pad': 7,
                          'ytick.major.pad': 7,
                          'ytick.color': nearly_black,
                          'xtick.color': nearly_black,
                          'axes.labelcolor': nearly_black,
                          #'legend.facecolor': light_grey,
                          'pdf.fonttype': 42,
                          'ps.fonttype': 42,
                          'mathtext.fontset': 'custom',
                          'font.size': fontsize,
                          'text.usetex' : True,
                          'font.family' : 'serif',
                          'font.serif' : 'Palatino Linotype',
                          'text.usetex': True,
                          'savefig.bbox':'tight',
                          'axes.facecolor': white,
                          'axes.labelpad': 10.0,
                          'axes.labelsize': fontsize,
                          'axes.titlepad': 30,
                          'lines.markersize': 5.0,
                          'lines.scale_dashes': False,
                          'xtick.labelsize': fontsize * 0.8,
                          'ytick.labelsize': fontsize * 0.8}

transitionfig_formatting = {'axes.formatter.limits': (-3,3),
                          'xtick.major.pad': 7,
                          'ytick.major.pad': 7,
                          'ytick.color': nearly_black,
                          'xtick.color': nearly_black,
                          'axes.labelcolor': nearly_black,
			  'axes.spines.bottom': True,
                          #'axes.spines.left': True,
                          #'legend.facecolor': light_grey,
                          'pdf.fonttype': 42,
                          'ps.fonttype': 42,
                          'mathtext.fontset': 'custom',
                          'font.size': fontsize,
                          'text.usetex' : True,
                          'font.family' : 'serif',
                          'font.serif' : 'Palatino Linotype',
                          'text.usetex': True,
                          'savefig.bbox':'tight',
                          'axes.facecolor': white,
                          'axes.labelpad': 10.0,
		          'axes.grid.axis': 'y',
			  'axes.grid': True,
                          'grid.color': light_grey,
                          'axes.labelsize': fontsize,
                          'axes.titlepad': 30,
                          'lines.markersize': 7.0,
                          'lines.scale_dashes': False,
                          'xtick.labelsize': fontsize * 0.8,
                          'ytick.labelsize': fontsize * 0.8}

compare_distinct_series = {'axes.formatter.limits': (-3,3),
                          'xtick.major.pad': 7,
                          'ytick.major.pad': 7,
                          'ytick.color': nearly_black,
                          'xtick.color': nearly_black,
                          'axes.labelcolor': nearly_black,
                          'axes.spines.bottom': False,
                          'axes.spines.left': False,
                          'axes.spines.right': False,
                          'axes.spines.top': False,
                          'axes.axisbelow': False,
                          'legend.frameon': False,
                          'pdf.fonttype': 42,
                          'ps.fonttype': 42,
                          'mathtext.fontset': 'custom',
                          'font.size': fontsize,
                          'text.usetex' : True,
                          'font.family' : 'serif',
                          'font.serif' : 'Palatino Linotype',
                          'text.usetex': True,
                          'savefig.bbox':'tight',
                          'axes.facecolor': light_grey,
                          'axes.labelpad': 10.0,
                          'axes.labelsize': fontsize,
                          'axes.titlepad': 30,
                          'axes.grid': True,
                          'axes.axisbelow': True,
                          'grid.color': white,
                          'lines.markersize': 7.0,
                          'lines.scale_dashes': False,
                          'xtick.labelsize': fontsize * 0.8,
                          'ytick.labelsize': fontsize * 0.8,
                          'legend.fontsize': fontsize * 0.8}

for k, v in formatting.items():
    rcParams[k] = v

color_cycle = tableau.values()
try:
    from matplotlib import cycler
    rcParams['axes.prop_cycle'] = cycler(color=color_cycle)
except Exception:
    raise
