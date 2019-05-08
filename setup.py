from distutils.core import setup

setup(
	name = 'fermimodel',
	packages = ['fermimodel',],
	version = '0.1',
	description = 'Tools for generating and editing models for use with the Fermitools.',
	keywords = 'fermitools gtobssim gtlike',
	url = 'https://github.com/jaulbric/fermimodel',
	author = 'J. F. Ulbricht',
	author_email = 'julbrich@ucsc.edu',
	install_requires = [
		'astropy',
		'scipy',
	],
	scripts = ['bin/make4FGLxml']
	)