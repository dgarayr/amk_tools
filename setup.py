from setuptools import setup,find_packages
setup(name='amk_tools',
      version='0.8',
      author="Diego Garay-Ruiz",
      description="RxN management for AutoMeKin",
      py_modules=['RXReader','RXVisualizer'],
      install_requires=['networkx','bokeh>=2.3.2','jsmol_bokeh_extension','matplotlib','numpy']
      )
