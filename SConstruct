import sys
from ACGenerateFile import * #Autoconf support! Important for SWIG

class colors:
	RED = "\033[91m"
	BLACK = "\033[0m"
	GREEN = "\033[92m"
	BLUE = "\033[94m"
	YELLOW = "\033[93m"

platforms = ['darwin','linux2','cygwin','os2']
platform = sys.platform

import os
copy_env = ['PATH','INCLUDE','LIB','TMP']
ENV = {}
for key in copy_env:
	if key in os.environ.keys():
		ENV[key] = os.environ[key]


env = Environment(SHLIBPREFIX="",ENV=ENV)

#Add command line options to scons script
build_types = ['debug','debug_optimized','release']
AddOption('--build',dest='build_type',type='string',nargs=1,action='store',default='debug')
AddOption('--interface',dest='interface_type',type='string',nargs=1,action='store',default='none')
build_interface = GetOption('interface_type')

build_type = GetOption('build_type')
if not (build_type in build_types):
	if platform in platforms:
		print colors.RED + "ERROR: expected build options are 'debug', 'debug_optimized' and 'release'" + colors.BLACK
	else:
		print "ERROR: expected build options are 'debug', 'debug_optimized' and 'release'"
	print
	quit(-1)


#add src path to Environment
src_path = ['#/src','#/External/','#/External/Eigen/','#/External/maxflow','#/External/Cephes']
env.Append(CPPPATH = src_path)

#set compiler flags for different build types
compiler_flags = []
debug_flags = []
optimize_flags = ['-O3','-DRELEASE']
if platform in platforms:
	compiler_flags.append('-fPIC')
	optimize_flags.extend(['-msse','-msse2'])
	debug_flags.extend(['-Wextra','-g','-Wall'])

env.Append(CCFLAGS=compiler_flags)
if build_type == 'debug':
	env.Append(CCFLAGS=debug_flags)
elif build_type == 'debug_optimized':
	env.Append(CCFLAGS=debug_flags)
	env.Append(CCFLAGS=optimize_flags)
else:
	env.Append(CCFLAGS=optimize_flags)

#env['STATIC_AND_SHARED_OBJECTS_ARE_THE_SAME']=1 

conf = Configure(env)

env.SConsignFile()

##BUILD static libraries
Export('env','src_path','platform','conf')

#build easycore lib
easycore = SConscript('#/src/CEasyGWAS/SConscript',variant_dir='#/build/' + platform + '/lib/CEasyGWAS/',duplicate=0)
env.Append(LIBS=easycore)

#build maxflow
maxflow = SConscript('#/External/maxflow/SConscript',variant_dir='#/build/' + platform + '/lib/maxflow/',duplicate=0)
env.Append(LIBS=maxflow)

#build cephes
cephes = SConscript('#/External/Cephes/SConscript',variant_dir='#/build/' + platform + '/lib/cephes/',duplicate=0)
env.Append(LIBS=cephes)

Export('easycore','maxflow','cephes')

##Build python interface
if build_interface=="python":
    CEasyGWAS = SConscript('#src/interfaces/python/CEasyGWAS/SConscript',variant_dir='#build/' + platform + '/interfaces/python/CEasyGWAS/',duplicate=0)

##BUILD test files
SConscript('#/src/tools/SConscript',variant_dir='#/build/' + platform + '/tools/',duplicate=0)

scones = SConscript('#/src/testing/scones/SConscript',variant_dir='#/build/' + platform + '/testing/scones/',duplicate=0)
env.Install('#bin/' + platform + '/testing/',scones)

#fastanova = SConscript('#/src/testing/fastanova/SConscript',variant_dir='#/build/' + platform + '/testing/fastanova/',duplicate=0)
#env.Install('#bin/' + platform + '/testing/',fastanova)

test = SConscript('#/src/testing/general/SConscript',variant_dir='#/build/' + platform + '/testing/general/',duplicate=0)
env.Install('#bin/' + platform + '/testing/',test)
