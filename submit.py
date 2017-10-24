import os

user    = 'varaujof'
version = '01'
athena  = '20.7.7.6'
rootversion = '6.08.06'
outputFile = 'analysis'

listfile = open('TILEMU.list.tmp')
datasets = listfile.readlines()
listfile.close()

for data in datasets:
    data = data.rstrip()
    command =  'prun --exec "TMDB %%IN %s" --bexec "make clean;make"' % outputFile
    command += ' --athenaTag=%s' % athena
    #ommand += ' --rootVer %s' % rootversion
    #command += ' --nGBPerJob=MAX'
    #command += ' --nFilesPerJob=1'
    #command += ' --forceStaged'
    command += ' --express'

    command += ' --inDS=group.det-muon.%s' % data
    command += ' --outDS=user.%s.%s.%s' % (user, data, version)

    command += ' --outputs "%s.root"' % outputFile
    command += '--mergeOutput'

    print command
    os.system(command)
