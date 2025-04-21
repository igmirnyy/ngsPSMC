#!/usr/bin/env python3

#    Copyright (c) 2018 Vladimir Shchur (vlshchur@gmail.com)
#
#    This file is part of MiSTI.
#
#    MiSTI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MiSTI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MiSTI.  If not, see <https://www.gnu.org/licenses/>.


import sys
import matplotlib.pyplot as plt
from math import log
import argparse
import random

class Units:#This is a class of static variables
    mutRate = 2.5e-8
    binsize = 100
    N0 = 10000
    genTime = 1
    firstCall = True
    def __init__(self, **kwargs):
        if "mutRate" in kwargs:
            Units.mutRate = kwargs["mutRate"]
        if "binsize" in kwargs:
            Units.binsize = kwargs["binsize"]
        if "N0" in kwargs:
            Units.N0 = kwargs["N0"]
        if "genTime" in kwargs:
            Units.genTime = kwargs["genTime"]
        if "inpFile" in kwargs:
            Units.SetUnitsFromFile(kwargs["inpFile"])
#        if Units.firstCall or len(kwargs) > 0:
            #print("Units")
            #print("mutation rate =", Units.mutRate, "\tbinsize =", Units.binsize, "\tN0 =", Units.N0, "\tgeneration time =", Units.genTime)
#            self.PrintUnits()
#            Units.firstCall = False

    def PrintUnits(self):
        print("Units: mutation rate =", Units.mutRate, "\tbinsize =", Units.binsize, "\tN0 =", Units.N0, "\tgeneration time =", Units.genTime)

    def SetUnitsFromFile(self, fn):
        try:
            with open(fn) as f:
                for line in f:
                    line = line.split("=")
                    if line[0] == "mutRate" and len(line) == 2:
                        try:#
                            Units.mutRate = float(line[1])
                        except:
                            print("Cannot read mutation rate entry from file, using default or previous values")
                            pass
                    elif line[0] == "binsize" and len(line) == 2:
                        try:
                            Units.binsize = float(line[1])
                        except:
                            print("Cannot read bin size entry from file, using default or previous values")
                            pass
                    elif line[0] == "N0" and len(line) == 2:
                        try:
                            Units.N0 = float(line[1])
                        except:
                            print("Cannot read N0 entry from file, using default or previous values")
                            pass
                    elif line[0] == "genTime" and len(line) == 2:
                        try:
                            Units.genTime = float(line[1])
                        except:
                            print("Cannot read generation time entry from file, using default values")
                            pass
        except:
            print("Units input file not found, using default values.")

class psmc_data:
    def __init__(self, name, filen, date, hetloss, rd):
        self.name = name
        self.file = filen
        self.date = date
        self.hetloss = hetloss
        self.rd = rd
class ms_data:
    def __init__(self,name, filen, date):
        self.name = name
        self.file = filen
        self.date = float(date)
        
def PrintErr(*args, sep="", endl="\n"):
    message = sep.join(args)
    message += endl
    sys.stderr.write(message)

def ReadPSMCFile(fn, RD = -1):
    maxRD = -1
    Tk = []
    Lk = []
    th = 0
    rh = 0
    
    with open(fn) as f:
        for line in f:
            line = line.split()
            if line[0] == "RD":
                maxRD = int( line[1] )
        if maxRD == -1:
            print("Corrupted or empty input file")
            sys.exit(0)
        if RD == -1 or RD > maxRD:
            RD = maxRD

    with open(fn) as f:  
        for line in f:
            line = line.split()
            if line[0] != "RD" or int( line[1] ) != RD:
                continue
            while line[0] != "RS":
                if line[0] == "TR":
                    th = float(line[1])
                    rh = float(line[2])
                line = next(f)
                line = line.split()
            while line[0] != "PA":
                if line[0] != "RS":
                    print("Unexpected line.")
                    sys.exit(0)
                Tk.append( float(line[2]) )
                Lk.append( float(line[3]) )
                line = next(f)
                line = line.split()
            break
    data = [Tk, Lk, RD, th, rh]
    return( data )

def ReadPSMC(psmc_objs: list[psmc_data], maxY = None):
    u = Units()
    theta = 4.0*u.binsize*u.mutRate*u.N0
    scaleTime = 2*u.genTime*u.N0
    scaleEPS = 1
    for psmc in psmc_objs:
        data = ReadPSMCFile(psmc.file, psmc.rd)
        data[3] = data[3]/(1.0-psmc.hetloss)
        data[0] = [v*data[3]/theta for v in data[0]]
        data[1] = [v*data[3]/theta for v in data[1]]
        sdResc = psmc.date/2/u.N0/u.genTime
        sdRect = psmc.date
        data[0] = [v + sdResc for v in data[0]]
        
        x = [v*scaleTime for v in data[0]]
        y = [scaleEPS*v for v in data[1]]
        if maxY is not None:
            y = [min(1.0/v,maxY) for v in y]
        AddToPlot(x, y, lbl = psmc.name)
#    plt.show()


def ReadMS(ms_objs: list[ms_data]):
    u = Units()
    for ms in ms_objs:
        with open(ms.file) as f:
            command_line = f.readline()
        command = command_line.split()
        t = float(command[command.index('-t')+1])
        nsites = float(command[command.index('-r')+2])
        simulation_N0 = t / nsites /4 / u.mutRate
        time_shift = ms.date/2/u.N0/u.genTime
            
            
        size_changes = command_line.split(' -eN ')
        times = []
        lambdas = []
        for i in size_changes[1:]:
            t_k, lambda_k = map(float, i.split()[:2])
            times.append(t_k * simulation_N0 * 4 * u.genTime + time_shift)
            lambdas.append(lambda_k * simulation_N0 / u.N0)
        times = [times[0] / 2]  + times
        lambdas = [simulation_N0 / u.N0] + lambdas
        AddToPlot(times, lambdas, "--", lbl = ms.name)

def AddToPlot(times, lambdas, style="-", lbl = "", id=1):
    #plt.figure(id)
    plt.step(times+[2*times[-1]], [lambdas[0]]+lambdas, style, alpha=0.7, label=lbl)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Migration inference from PSMC.')
    parser.add_argument('-psmc', nargs='*', action = 'append',
                        help='psmc data: sample name, file name, sample date, sample hetloss, psmc round (optional).')#-mi [npop:1/2] [migStart] 
    parser.add_argument('-ms', nargs='*', action = 'append',
                        help='ms simulation data: sample name, file name, sample date')
    parser.add_argument('--funits', '-u', nargs=1, type=str, default="setunits.txt",
                        help='File name with units to be used to rescale times and EPS.')
    parser.add_argument('--fout', '-o', nargs=1, type=str, default="plot.pdf",
                        help='Output file.')
    parser.add_argument('--maxY', '-Y', nargs=1, type=float, default=None,
                        help='Range for Y axis (upper bound).')
    clargs = parser.parse_args()

    if isinstance(clargs.funits, list):
        clargs.funits = clargs.funits[0]
    if isinstance(clargs.fout, list):
        clargs.fout = clargs.fout[0]
    if isinstance(clargs.maxY, list):
        clargs.maxY = clargs.maxY[0]
        
    units = Units()
    units.SetUnitsFromFile(clargs.funits)
    units.PrintUnits()


    if clargs.psmc == None and clargs.ms == None:
        PrintErr("Provide at least one -psmc or -ms argument.")
        sys.exit(0)
    clargs.psmc = clargs.psmc or []
    clargs.ms = clargs.ms or []
    plt.semilogx()
    psmcs = []
    for el in clargs.psmc:
        if len(el) == 4:
            el.append(-1)
        if len(el) != 5:
            PrintErr("psmc Argument should have 4 or 5 parameters.")
            sys.exit(1)
        psmc_d = psmc_data(el[0], el[1], float(el[2]), float(el[3]), int(el[4]))
        psmcs.append(psmc_d)
    mss = []
    for el in clargs.ms:
        if len(el) != 3:
            PrintErr("ms Argument should have 3 parameters.")
            sys.exit(1)
        ms_d = ms_data(el[0], el[1], el[2])
        mss.append(ms_d)
    
    ReadPSMC(psmcs)
    ReadMS(mss)
    plt.ylim(bottom=0)
    if clargs.maxY != None:
        plt.ylim(top=clargs.maxY)
    plt.legend()
    plt.savefig(clargs.fout)