from distutils.core import setup, Extension

module = Extension('TestMod', sources = ["test.c"])

setup (name = 'TestMod',
       version = '1.0', 
       description = '',
       ext_modules = [module]
       )