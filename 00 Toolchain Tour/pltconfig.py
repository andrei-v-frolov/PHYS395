# matplotlib configuration for publication-quality output

from matplotlib import rc, rcParams

# custom style for matplotlib plots
params = {
	'axes.labelsize': 16,
	'legend.fontsize': 16,
	'xtick.labelsize': 16,
	'ytick.major.pad': 11,
	'xtick.major.pad': 11,
	'ytick.labelsize': 16,
	'font.size': 16,
	'font.family': 'sans-serif',
	'font.sans-serif': 'FreeSans'
}

rcParams.update(params)

# output plot as PDF file
output = {
	'backend': 'pdf',
	'savefig.dpi': 300
}

rcParams.update(output)

# render labels using LaTeX
rcParams.update({'text.usetex': True})
rc('text.latex', preamble=r'\usepackage{sfmath}')
