#!/usr/bin/env python

# Copyright (C) 2004 Rune Linding & Lars Juhl Jensen - EMBL
# The DisEMBL is licensed under the GPL license
# (http://www.opensource.org/licenses/gpl-license.php)
# DisEMBL pipeline

# Modified to work with current versions of Biopython (1.7+)
# by Shyam Saladi (saladi1@illinois.edu), Janauary 2013 
# Bio:SeqIO completely replaces Bio:Fasta

from string import *
from sys import argv
from Bio import File
from Bio import SeqIO
import fpformat
import sys
import tempfile
import os
from os import system,popen3

# change these to the correct paths
NN_bin = os.environ['NN_bin']
SG_bin = os.environ['SG_bin']

def JensenNet(sequence):
    outFile = tempfile.mktemp()
    inFile= tempfile.mktemp()
    open(inFile,'w').write(sequence+'\n')
    system(NN_bin + '< ' + inFile +' > ' + outFile)
    REM465 = []
    COILS = []
    HOTLOOPS = []
    resultsFile = open(outFile,'r')
    results = resultsFile.readlines()
    resultsFile.close()
    for result in results:
        coil = float(fpformat.fix(split(result)[0],6))
        COILS.append(coil)
        hotloop = float(fpformat.fix(split(result)[1],6))
        HOTLOOPS.append(hotloop)
        rem465 = float(fpformat.fix(split(result)[2],6))
        REM465.append(rem465)
    os.remove(inFile)
    os.remove(outFile)
    return COILS, HOTLOOPS, REM465


def SavitzkyGolay(window,derivative,datalist):
    if len(datalist) < 2*window:
        window = len(datalist)/2
    elif window == 0:
        window = 1
    stdin, stdout, stderr = popen3(SG_bin + ' -V0 -D' + str(derivative) + ' -n' + str(window)+','+str(window))
    for data in datalist:
        stdin.write(`data`+'\n')
    try:
        stdin.close()
    except:
        print stderr.readlines()
    results = stdout.readlines()
    stdout.close()
    SG_results = []
    for result in results:
        f = float(fpformat.fix(result,6))
        if f < 0:
            SG_results.append(0)
        else:
            SG_results.append(f)
    return SG_results

def getSlices(NNdata, fold, join_frame, peak_frame, expect_val):
    slices = []
    inSlice = 0
    for i in range(len(NNdata)):
        if inSlice:
            if NNdata[i] < expect_val:
                if maxSlice >= fold*expect_val:
                    slices.append([beginSlice, endSlice])
                inSlice = 0
            else:
                endSlice += 1
                if NNdata[i] > maxSlice:
                    maxSlice = NNdata[i]
        elif NNdata[i] >= expect_val:
            beginSlice = i
            endSlice = i
            inSlice = 1
            maxSlice = NNdata[i]
    if inSlice and maxSlice >= fold*expect_val:
        slices.append([beginSlice, endSlice])

    i = 0
    while i < len(slices):
        if i+1 < len(slices) and slices[i+1][0]-slices[i][1] <= join_frame:
            slices[i] = [ slices[i][0], slices[i+1][1] ]
            del slices[i+1]
        elif slices[i][1]-slices[i][0]+1 < peak_frame:
            del slices[i]
        else:
            i += 1
    return slices


def reportSlicesTXT(slices, sequence):
    if slices == []:
        s = lower(sequence)
    else:
        if slices[0][0] > 0:
            s = lower(sequence[0:slices[0][0]])
        else:
            s = ''
        for i in range(len(slices)):
            if i > 0:
                sys.stdout.write(', ')
            sys.stdout.write( str(slices[i][0]+1) + '-' + str(slices[i][1]+1) )
            s = s + upper(sequence[slices[i][0]:(slices[i][1]+1)])
            if i < len(slices)-1:
                s = s + lower(sequence[(slices[i][1]+1):(slices[i+1][0])])
            elif slices[i][1] < len(sequence)-1:
                s = s + lower(sequence[(slices[i][1]+1):(len(sequence))])
    print ''
    print s



def runDisEMBLpipeline():
    try:
        smooth_frame = int(sys.argv[1])
        peak_frame = int(sys.argv[2])
        join_frame = int(sys.argv[3])
        fold_coils = float(sys.argv[4])
        fold_hotloops = float(sys.argv[5])
        fold_rem465 = float(sys.argv[6])
        file = str(sys.argv[7])
        try:
            mode = sys.argv[8]
        except:
            mode = 'default'
    except:
        print '\nDisEMBL.py smooth_frame peak_frame join_frame fold_coils fold_hotloops fold_rem465 sequence_file [mode]\n'
        print 'A default run would be: ./DisEMBL.py 8 8 4 1.2 1.4 1.2  fasta_file'
        print 'Mode: "default"(nothing) or "scores" which will give scores per residue in TAB seperated format'
        raise SystemExit
    db = open(file,'r')
    print ' ____  _     _____ __  __ ____  _       _  _  _'
    print '|  _ \(_)___| ____|  \/  | __ )| |     / || || |'
    print '| | | | / __|  _| | |\/| |  _ \| |     | || || |_'
    print '| |_| | \__ \ |___| |  | | |_) | |___  | ||__   _|'
    print '|____/|_|___/_____|_|  |_|____/|_____| |_(_) |_|'
    print '# Copyright (C) 2004 - Rune Linding & Lars Juhl Jensen '
    print '# EMBL Biocomputing Unit - Heidelberg - Germany        '
    print '#'
    for cur_record in SeqIO.parse(db, "fasta"):
	    sequence = upper(str(cur_record.seq.tostring()))
            # Run NN
            COILS_raw, HOTLOOPS_raw, REM465_raw = JensenNet(sequence)
            # Run Savitzky-Golay
            REM465_smooth = SavitzkyGolay(smooth_frame,0,REM465_raw)
            COILS_smooth = SavitzkyGolay(smooth_frame,0,COILS_raw)
            HOTLOOPS_smooth = SavitzkyGolay(smooth_frame,0,HOTLOOPS_raw)
            if mode == 'default':
                sys.stdout.write('> '+cur_record.id+'_COILS ')
                reportSlicesTXT( getSlices(COILS_smooth, fold_coils, join_frame, peak_frame, 0.43), sequence )
                sys.stdout.write('> '+cur_record.id+'_REM465 ')
                reportSlicesTXT( getSlices(REM465_smooth, fold_rem465, join_frame, peak_frame, 0.50), sequence )
                sys.stdout.write('> '+cur_record.id+'_HOTLOOPS ')
                reportSlicesTXT( getSlices(HOTLOOPS_smooth, fold_hotloops, join_frame, peak_frame, 0.086), sequence )
                sys.stdout.write('\n')
            elif mode == 'scores':
                sys.stdout.write('# RESIDUE COILS REM465 HOTLOOPS\n')
                for i in range(len(REM465_smooth)):
                    sys.stdout.write(sequence[i]+'\t'+fpformat.fix(COILS_smooth[i],5)+'\t'+fpformat.fix(REM465_smooth[i],5)+'\t'+fpformat.fix(HOTLOOPS_smooth[i],5)+'\n')
            else:
                sys.stderr.write('Wrong mode given: '+mode+'\n')
                raise SystemExit
    db.close()
    return

runDisEMBLpipeline()
