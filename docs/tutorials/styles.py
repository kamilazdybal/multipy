# ------------------------------------------------------------------------------
#
# This file sets plotting parameters and styles.
#
# ------------------------------------------------------------------------------

from matplotlib import rcParams

csfont = {'fontname':'Charter', 'fontweight':'regular'}
hfont = {'fontname':'Charter', 'fontweight':'bold'}
ifont = {'fontname':'Charter', 'fontweight':'regular', 'style':'italic'}
rcParams["font.family"] = "serif"
rcParams["font.serif"] = "Charter"
rcParams["font.sans-serif"] = "Charter"
rcParams["font.cursive"] = "Charter"
rcParams["font.monospace"] = "Charter"
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Charter'
rcParams['mathtext.it'] = 'Charter:italic'
rcParams['mathtext.bf'] = 'Charter:bold'
rcParams['font.size'] = 12
rcParams["text.usetex"] = False

grid_opacity = 0.3

font_axes = 16
font_labels = 20
font_annotations = 20
font_title = 20
font_text = 16
font_legend = 16
font_colorbar = 24
font_colorbar_axes = 18

marker_size = 50
marker_scale_legend = 1
marker_scale_legend_clustering = 10

scatter_point_size = 2

line_width = 1

eigenvector_bar_width = 0.4

color_species_1 = '#C7254E'
color_species_2 = '#BBBBBB'
color_species_3 = '#008CBA'
color_species_4 = '#ff9e1b'

colors = [color_species_1, color_species_2, color_species_3, color_species_4]

color_greys = ['#e1e1e1', '#c9c9c9', '#929292', '#4b4b4b', '#1e1e1e']
