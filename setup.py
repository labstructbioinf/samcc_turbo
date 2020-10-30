from distutils.core import setup

install_reqs = ['biopython>=1.76', 'networkx==2.2', 'seaborn>=0.10.0']

setup(
  name = 'samcc-turbo',
  packages = ['samcc'],
  version = '0.0.2',
  license='MIT',
  description = 'Software for automatic detection and measurement of coiled coils in PDB structures.',
  author = 'kszczepaniak',
  author_email = 'k.szczepaniak@cent.uw.edu.pl',
  url = 'https://github.com/labstructbioinf/samcc_turbo',
  download_url = 'https://github.com/labstructbioinf/samcc_turbo/archive/v0.0.2.tar.gz',
  keywords = ['bioinformatics', 'protein structure', 'coiled-coil'],
  install_requires=install_reqs,
  classifiers=[
    'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    # 'Intended Audience :: Developers',      # Define that your audience are developers
    # 'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    # 'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    # 'Programming Language :: Python :: 3.4',
    # 'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.7',
  ],
)
