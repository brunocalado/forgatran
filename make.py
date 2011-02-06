#!/usr/bin/python3
# -*- coding: UTF-8 -*-

#----------------------------------------------------------
# Compile and generate the executable for the binary gene
# sample
#----------------------------------------------------------




#----------------------------------------------------------
# Variables and Constants
#----------------------------------------------------------
import os
import sys
      



#----------------------------------------------------------
# Variables and Constants
#----------------------------------------------------------
#dir receives the actual path 
#dir = os.getcwd()

compilerOptionsFileName = 'intelCompilerOptions.cfg'
sourceFilesFileName = 'sourceFiles.cfg'

compiler = 'ifort'
extension = 'f90'

compilerOptions = []
sourceDirs = [ 'app/', 'genes/', 'src/' ]
sourceFiles = []
selectedSourceFiles = []
exeFilename = 'FORGATRAN'


out = ''
outSrc= ''
debug = 1

if __name__ == '__main__':
    if len(sys.argv) == 2: 
        if '-' in sys.argv[1]: 
            if sys.argv[1] == '-config': 
 
                fileHandle = open( sourceFilesFileName, 'w')
                for i in sourceDirs:
                    for j in os.listdir(i):
                        if extension in j:
                            fileHandle.writelines( '1 ' + i+j  + '\n')
                fileHandle.close()
                sys.exit(None)
            elif 'debug' in sys.argv[1]: 
                debug = 1
            else:
                exeFilename = sys.argv[1]  


        else:
            compilerOptionsFileName = sys.argv[1]
 


    fileHandle = open( sourceFilesFileName, 'r')
    sourceFiles = fileHandle.read().splitlines()
    fileHandle.close()

    if debug:
        print( '\n' + 50*'-' )
        print( 'Files to be compiled:' )  
        
    for i in range(0, len(sourceFiles)):
        if '1' in sourceFiles[i].split()[0]:
            if debug: print( '1 -> ' + sourceFiles[i].split()[1] )
            selectedSourceFiles.append(sourceFiles[i].split()[1])
        else:
            if debug: print( '0 -> ' + sourceFiles[i].split()[1] )

    if debug: print( 50*'-' + '\n' )



#----------------------------------------------------------
# Compiler Options
#----------------------------------------------------------
#Open files with the compiler flags
fileHandle = open( compilerOptionsFileName, 'r')
compilerOptions = (fileHandle.read()).splitlines()
fileHandle.close()

#Temporary variable
out = compiler


if debug:
    print( '\n' + 50*'-' )
    print( 'Selected options:' )  
for i in compilerOptions:
    if debug: print( i )
    out += ' ' + i
out += ' ' + exeFilename
if debug: print( 50*'-' + '\n' )





#----------------------------------------------------------
# Main
#----------------------------------------------------------
#outtmp = out
#for j in range(0,3):
#out = outtmp

for i in selectedSourceFiles:
    out += ' ' + i

os.system( out )
os.system( out )
os.system( out )
os.system( out )
os.system( out )
os.system( out )


if debug: print(out)
if debug: print( 'Success' )


os.system( 'rm -f *.mod' )
os.system( 'mv ' + exeFilename + ' tmp' )



