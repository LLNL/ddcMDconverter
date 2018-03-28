from setuptools import setup, find_packages

setup(name='ddcmdconverter',
      description='A tool for converting other simulation formats to ddcMD'
                  'compatible formats.',
      version='0.0.1dev',
      author='Xiaohua Zhang',
      author_email='zhang30@llnl.gov',
      packages=find_packages(),
      entry_points={
        'console_scripts': [
            'martiniobj2pdb = ddcmdconverter.objmartini2pdb:main',
            'martini2obj = ddcmdconverter.martini2obj:main',
            'pdbmartini2obj = ddcmdconverter.pdbmartini2obj:main',
            'pdb2obj = ddcmdconverter.pdb2obj:main',
            'obj2pdb = ddcmdconverter.obj2pdb:main',
            'restraint = ddcmdconverter.restraint:main',
        ]
      },
      install_requires=[],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        ],
      )
