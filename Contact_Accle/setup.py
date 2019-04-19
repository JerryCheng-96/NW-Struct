from distutils.core import setup, Extension

module = Extension('ContactAccel', sources = ["Contact_Accle.c"])

setup (name = 'ContactAccel',
       version = '1.0', 
       description = '',
       ext_modules = [module]
       )
