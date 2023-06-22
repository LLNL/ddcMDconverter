from setuptools import setup, find_packages

setup(name='ddcmdconverter',
      description='A tool for converting other simulation formats to ddcMD'
                  'compatible formats.',
      version='1.0.0dev',
      author='Xiaohua Zhang',
      author_email='zhang30@llnl.gov',
      packages=find_packages(),
      entry_points={
        'console_scripts': [
            'martiniobj2pdb = ddcmdconverter.martini.objmartini2pdb:main',
            'martini2obj = ddcmdconverter.martini.martini2obj:main',
            'pdbmartini2obj = ddcmdconverter.martini.pdbmartini2obj:main',
            'pdb2obj = ddcmdconverter.charmm.pdb2obj:main',
            'obj2pdb = ddcmdconverter.charmm.obj2pdb:main',
            'restraint = ddcmdconverter.martini.restraint:main',
            'ddcmdconverter = ddcmdconverter.gmx.ddcmdParameter:main',
        ]
      },
      install_requires=[],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        ],
      )
