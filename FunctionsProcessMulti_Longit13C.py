import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from sklearn.metrics import auc
import os 
from tkinter import simpledialog
from operator import contains
from tkinter.filedialog import askdirectory
import tkinter as tk
from ipywidgets import interact, widgets, fixed, interactive, HBox, Layout
import pickle
import pybaselines
import matplotlib as mpl
import pandas as pd

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Funtion to select directory with the data and extract subexperiments folders automatically for the processing
def getDirsAll(intdr = 'C:/IBECPostDocDrive/2023_10_04_ProcessSlabsChip/LongitudinalExps'):
    # GUI to select the experiment directory desired
    root = tk.Tk()
    # root.withdraw()
    root.attributes('-topmost',True)
    root.iconify()

    test = askdirectory(parent=root, initialdir=intdr)

    root.destroy()
    
    
    rml = ["99999","AdjProtocols","AdjResult","AdjStatePerStudy","ResultState","ScanProgram.scanProgram","study.MR","subject","Plots", "Results"] # Fixed dirs to remove since they contain no experiments
    fldslst = os.listdir(test) # Get list of dirs from selected folder
    [fldslst.remove(rml[i]) for i in range(len(rml)) if rml[i] in fldslst] # Remove undesired dirs
    
    
    fldslst_tmp = [fldslst[i] for i in range(len(fldslst)) if fldslst[i][0:3] == 'raw']

    if not fldslst_tmp:
        machi = 'MRI' # Our Bruker MRI sytem
    else:
        machi = 'Pulsar' # Our benchtop spectrometer
        
    if machi == 'Pulsar':
        fldslst2 = [int(x[3::]) for x in fldslst_tmp] # Convert string list into numbers to sort them as experiments happen
        fldslst2.sort() # Sort dirs
        fldslst3 = ['raw'+str(x) for x in fldslst2] # Go back to dir as strings instead of numericals
        
        # Loop thorugh all the directories and extract the scan name to identify the singlepulses of interest
        dircns = dict()
        for i in range(len(fldslst3)):
            dircns[fldslst3[i]] = 'Pulsar_Scan'
            
        lngrps = 1   

        lgtdt = {}
        lgtdt["0"] = fldslst3
        stck = True
        
    elif machi == 'MRI':
    
    
        fldslst2 = [int(x) for x in fldslst] # Convert string list into numbers to sort them as experiments happen
        fldslst2.sort() # Sort dirs
        fldslst3 = [str(x) for x in fldslst2] # Go back to dir as strings instead of numericals
        
        # Loop thorugh all the directories and extract the scan name to identify the singlepulses of interest
        dircns = dict()
        for i in range(len(fldslst3)):
            with open(test+"/"+fldslst3[i]+"/acqp", 'r+') as f:
                for line in f:
                    if 'ACQ_scan_name' in line.strip():
                        dircns[fldslst3[i]] = next(f,'')[1:-2]
                        break
                    
        # Check first 13C Scan to see if it is a single scan or a stack
        insc = ["Singlepulse_13C" in dircns[fldslst3[i]] for i in range(len(fldslst3))].index(True)
        NR = []
        with open(test+'/'+fldslst3[insc]+'/acqp', 'r+') as f:
            for line in f:
                if '##$NR=' in line.strip():
                    # NR = float(next(f,''))
                    NR = int(line.strip()[6:])
                    break       
                    
        if NR == 1:   
            stck = False         
            lngrps = 0 # Number of longitudinal experiments in one experiment folder
            for i in range(len(dircns)-1):
                if ("Localizer" in dircns[fldslst3[i]] and "Singlepulse_13C" in dircns[fldslst3[i+1]]) or ("Localized" in dircns[fldslst3[i]] and "Singlepulse_13C" in dircns[fldslst3[i+1]]) or ("T1" in dircns[fldslst3[i]] and "Singlepulse_13C" in dircns[fldslst3[i+1]]) or ("T2" in dircns[fldslst3[i]] and "Singlepulse_13C" in dircns[fldslst3[i+1]]):
                    lngrps += 1
                
            # Put in a dictionary with as many entries as longitudinal experiments the directory names for each pertinent single pulse 13C scan
            lgtdt = dict()
            tmprm = 0
            for i in range(lngrps):
                lgtdt[str(i)] = []
                for j in range(tmprm, len(dircns)-1):
                    if "Singlepulse_13C" in dircns[fldslst3[j]] and "Singlepulse_13C" in dircns[fldslst3[j+1]]:
                        lgtdt[str(i)].append(fldslst3[j])
                        if j == len(dircns)-2:
                            lgtdt[str(i)].append(fldslst3[j+1])
                    elif "Singlepulse_13C" in dircns[fldslst3[j]] and not contains(dircns[fldslst3[j+1]], "Singlepulse_13C"):
                        lgtdt[str(i)].append(fldslst3[j])
                        tmprm = j+1
                        break
                    
            for i in range(lngrps):
                if not lgtdt[str(i)]:
                    lngrps -= 1
                    del lgtdt[str(i)]
        else:
            stck = True
            lngrps = sum(["Singlepulse_13C" in dircns[fldslst3[i]] for i in range(len(fldslst3))])
            lgtdt = dict()
            tmprm = 0
            for i in range(lngrps):
                lgtdt[str(i)] = []
                for j in range(tmprm, len(dircns)):
                    if "Singlepulse_13C" in dircns[fldslst3[j]]:
                        lgtdt[str(i)].append(fldslst3[j])
                        tmprm = j+1
                        break
            
    # This is to check if the user has processed this already and extract some info. 
    printInfo(test)
    
    return(test, fldslst3, dircns, lngrps, lgtdt, stck, machi)



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Function to ask the user the number of wells scaned in the chip to loop over the number later on
def getWellNum(stck):

    if stck == False:
        WELLS = -1
        while (WELLS > 16 or WELLS < 0):
            WELLS = simpledialog.askstring("Wells", "Enter Number Of Wells Used in Experiment:")
            try:
                WELLS = int(WELLS)
            except:
                WELLS = -1
    else:
        WELLS = 1
            
    return(WELLS)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Function to extract necessary metadata for each single pulse scan folder for each subexperiment and each well
def extrMDAll(lngrps, test, lgtdt, machi):
    metDat = {}
    indx = [i for i in range(len(test)) if test[i].find('/')==0]
    foldPath = test[0::]

    # Generate directory to save all plots
    if os.path.exists(foldPath+'/Plots') == False:
        os.mkdir(foldPath+'/Plots')
        
    if os.path.exists(foldPath+'/Results') == False:
        os.mkdir(foldPath+'/Results')
        
    if machi == 'MRI':
        for j in range(lngrps):
            metDat[str(j)] = {}
            for i in lgtdt[str(j)]:                
                rot, RECO_ft_mode, bw, bwc, ACQ_repetition_time, ACQ_Size, NR, bwHz, ti = getMetaDat(foldPath, str(i))
                metDat[str(j)][str(i)] = rot, RECO_ft_mode, bw, bwc, ACQ_repetition_time, ACQ_Size, NR, bwHz, ti
    elif machi == 'Pulsar':
        metDat[str(0)] = {}
        for i in range(len(lgtdt[str(0)])):
            metDat[str(0)][lgtdt['0'][i]] = getMetaDat_Pulsar(foldPath, lgtdt['0'][i])
            
            
    return(metDat, indx, foldPath)



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to extract and compute the time vector for each time well and for each sub-experiment
def getTimeVec(lngrps, numWells, lgtdt, metDat, machi):
    
    timVecs = {}
    
    if machi == 'MRI':    
        if metDat[str(0)][lgtdt[str(0)][0]][6] == 1:
        
            for j in range(lngrps):
                
                timVecs[str(j)] = {}
                
                for i in range(numWells):
                    dds = lgtdt[str(j)][0+i:len(lgtdt[str(j)]):numWells]

                    fols = [int(dds[s]) for s in range(len(dds))]
                    tmptim = [metDat[str(j)][str(fols[k])][-1] for k in range(len(fols))]
                    
                    tm0 = [metDat[str(j)][str(lgtdt[str(j)][0])][-1] for k in range(len(fols))]
                    tm00 = [float(tm0[0].split(' ')[2].split(":")[i]) for i in range(3)]
                    tm000 = (tm00[0]*(60*60)+tm00[1]*60+tm00[2])
                    
                    timstp0 = [float(tmptim[0].split(' ')[2].split(":")[i]) for i in range(3)]
                    timstp0_cor = (timstp0[0]*(60*60)+timstp0[1]*60+timstp0[2])
                    
                    timstp = [[float(tmptim[k].split(' ')[2].split(":")[i]) for i in range(3)] for k in range(len(tmptim))]
                    timstp_cor = [(timstp[k][0]*(60*60)+timstp[k][1]*60+timstp[k][2]) for k in range(len(tmptim))]
                    
                    timVecs[str(j)][str(i)] = [timstp_cor[k]-tm000 for k in range(len(tmptim))]
                    
        else:
            for j in range(lngrps):
                timVecs[str(j)] = {}
                timVecs[str(j)]["0"] = list(np.arange(0,metDat["0"][lgtdt["0"][0]][4]*float(metDat["0"][lgtdt["0"][0]][6]), metDat["0"][lgtdt["0"][0]][4]))
    elif machi == 'Pulsar':
        timVecs['0'] = {}
        timVecs['0']['0'] = np.arange(0, len(lgtdt['0'])*5, 5)
            
    return(timVecs)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to ask the user the factor used in the exponential for appodisation
def getLB():

    LB = -1
    while (LB > 100 or LB < 0):
        LB = simpledialog.askstring("Wells", "Enter Appodisation Factor (0 to 100):")
        try:
            LB = float(LB)
        except:
            LB = -1
            
    return(LB)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Function to extract the raw spectroscopic data, rearrange it into an FID, do MRI phase corrections and generate a fourier transform
def ProcDatMult_Long(lgtdt, lngrps, metDat, numWells, foldPath, lineborad, stck = False,  xin = 76):
    datFut = {}
    datFutmag = {}
    datco3 = {}
    foldsDis = {}

    ppms = np.linspace(metDat["0"][lgtdt["0"][0]][3]-metDat["0"][lgtdt["0"][0]][2]/2, metDat["0"][lgtdt["0"][0]][3]+metDat["0"][lgtdt["0"][0]][2]/2, metDat["0"][lgtdt["0"][0]][5])

    for j in range(lngrps):
        datFut[str(j)] = {}
        datFutmag[str(j)] = {}
        datco3[str(j)] = {}
        foldsDis[str(j)] = {}
        
        
        foldNum = int(lgtdt[str(j)][0])
        foldLas = int(lgtdt[str(j)][-1])    
        if stck == False: # THis is if scans are in separated scans
            for i in range(numWells):
                fols = lgtdt[str(j)][0+i:len(lgtdt[str(j)]):numWells]
                foldsDis[str(j)][str(i)] = fols
                
                tmpdat = np.empty((len(fols), metDat[str(j)][str(foldNum)][5]*2))
                re = np.empty((len(fols), metDat[str(j)][str(foldNum)][5]))
                im = np.empty((len(fols), metDat[str(j)][str(foldNum)][5]))

                co = np.ndarray((len(fols), metDat[str(j)][str(foldNum)][5]), dtype=np.complex128)

                co3F = np.ndarray((len(fols), metDat[str(j)][str(foldNum)][5]), dtype=np.complex128)
                futF = np.ndarray((len(fols), metDat[str(j)][str(foldNum)][5]), dtype=np.complex128)
                futmagF = np.ndarray((len(fols), metDat[str(j)][str(foldNum)][5]), dtype=np.float32)

                insrmv = []
                for k in range(len(fols)):
                    try:
                        tmpdat[k,::] = np.fromfile(foldPath+'/'+fols[k]+'/pdata/1/fid_proc.64', np.float64)
                        
                        re[k,::] = tmpdat[k,0::2]
                        im[k,::] = tmpdat[k,1::2]

                        for m in range(0,len(re[k,::])):
                            co[k,m] = complex(re[k,m], im[k,m])
                            
                        co2 = np.reshape(co[k,::], (metDat[str(j)][str(fols[k])][6],metDat[str(j)][str(fols[k])][5]))
                        co3 = np.reshape(co[k,::], (metDat[str(j)][str(fols[k])][6],metDat[str(j)][str(fols[k])][5]))
                        
                        co2ec = np.concatenate([co2[0,xin::], co2[0,0:xin]])
                        co3ec = np.concatenate([co3[0,xin::], co3[0,0:xin]])
                        
                        
                        
                        lb = np.exp(-np.linspace(start=0, stop=lineborad, num=len(co2ec)))
                        co2ec[::] = co2ec[::]*lb
                        
                        co3F[k,::] = co2ec

                        RECO_rotate = [metDat[str(j)][str(fols[k])][0]]
                        dims = np.shape(co2ec)[0]
                        phase_matrix = np.ones(dims)
                        
                        f=np.array(range(dims))
                        phase_vector=np.exp(complex(0,1)*2*np.pi*RECO_rotate[0]*f)

                        fut = fft(co2ec*phase_vector) 
                        
                        futmag = [np.sqrt(np.real(fut[s])**2 + np.imag(fut[s])**2) for s in range(np.shape(co2[0])[0])]
                        
                        futF[k,::] = fut
                        futmagF[k,::] = futmag
                    except:
                        insrmv.append(k)

                for m in reversed(range(len(insrmv))):
                    del foldsDis[str(j)][str(i)][insrmv[m]]

                # print(fols[insrmv])
                datFut[str(j)][str(i)] = futF
                datFutmag[str(j)][str(i)] = futmagF
                datco3[str(j)][str(i)] = co3F
                
        else: # this is if scans are in the same scan (repetitions)
            file_content = np.fromfile(foldPath+'/'+str(foldNum)+'/pdata/1/fid_proc.64', np.float64)
            NR = metDat[str(j)][lgtdt[str(j)][0]][6]
            ACQ_Size = metDat[str(j)][lgtdt[str(j)][0]][5]
            
            fols = lgtdt[str(j)][0]
            foldsDis[str(j)][str(0)] = [fols]*NR
            
            re = file_content[0::2]
            im = file_content[1::2]
            co = np.ndarray(len(re), dtype=np.complex128)

            for i in range(0,len(re)):
                co[i] = complex(re[i], im[i])
                
            co2 = np.reshape(co, (NR,ACQ_Size))
            co3 = np.reshape(co, (NR,ACQ_Size))

            co2ec = np.zeros(np.shape(co2), dtype=np.complex128)

            for m in range(np.shape(co2)[0]):
                co2ec[m] = np.concatenate([co2[m,xin::], co2[m,0:xin]])

            lb = np.exp(-np.linspace(start=0, stop=lineborad, num=len(co2[0])))
            for l in range(0,NR):
                co2ec[l] = co2ec[l]*lb
                
            RECO_rotate = [metDat[str(j)][lgtdt[str(j)][0]][0]]
            dims = np.shape(co2ec)[1]
            phase_matrix = np.ones(dims)

            # for index in range(0,dims):
            f=np.array(range(dims))
            phase_vector=np.exp(complex(0,1)*2*np.pi*RECO_rotate[0]*f)

            # fut = [fft(co2ec[i,:]*phase_vector) for i in range(np.shape(co2ec)[0])]
            fut = [fft(co2ec[i,:]*phase_vector) for i in range(np.shape(co2ec)[0])]
            futmag = [np.sqrt(np.real(fut[i])**2 + np.imag(fut[i])**2) for i in range(np.shape(co2ec)[0])]
            
            datFut[str(j)][str(0)] = fut
            datFutmag[str(j)][str(0)] = futmag
            datco3[str(j)][str(0)] = co2ec
            
            
            
            
    return(ppms, datFutmag, datFut, datco3, foldsDis)
    

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to switch between the processing of MRI data and NMR benchtop data
def procDatGen(lgtdt, lngrps, metDat, numWells, foldPath, lineborad, machi, stck = False,  xin = 76):
    if machi == 'MRI':
        ppms, datFutmag, datFut, datco3, foldsDis = ProcDatMult_Long(lgtdt, lngrps, metDat, numWells, foldPath, lineborad, stck, xin)
    elif machi == 'Pulsar':
        ppms, datFutmag, datFut, datco3, foldsDis = ProcDatMult_Pulsar(lgtdt, metDat, foldPath, lineborad)
        
    return(ppms, datFutmag, datFut, datco3, foldsDis)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Function to extract the raw data from the Pulsar benchtop NMR and generate the furier transform
def ProcDatMult_Pulsar(lgtdt, metDat, foldPath, lineborad):
    datFut = {}
    datFutmag = {}
    datco3 = {}
    foldsDis = {}


    RecPnts, SF, O1, Filter = metDat['0']['raw1']
    freq = np.linspace(-Filter/2, Filter/2, RecPnts)
    timaqu = np.linspace(1, 1/(freq[1]-freq[0]), RecPnts)
    ppms_tmp = ((freq - O1) / SF)


    datFut[str(0)] = {}
    datFutmag[str(0)] = {}
    datco3[str(0)] = {}
    foldsDis[str(0)] = {}

    co = np.ndarray((len(metDat['0']), RecPnts), dtype=np.complex128)
    co3F = np.ndarray((len(metDat['0']), RecPnts), dtype=np.complex128)
    futF = np.ndarray((len(metDat['0']), RecPnts), dtype=np.complex128)
    futmagF = np.ndarray((len(metDat['0']), RecPnts), dtype=np.float32)

    mxmag = []

    for i in range(len(metDat['0'])):
        foldNum2 = lgtdt['0'][i]
        
        fls = os.listdir(foldPath+'/'+foldNum2+'/')
        fils = {}
        for j in range(len(fls)):
            if fls[j][-3::] == 'par':
                fils['pars'] = fls[j]
            elif fls[j][-3::] == 'fid':
                fils['dats'] = fls[j]
                
        rawtmpdat = np.fromfile(foldPath+'/'+foldNum2+'/'+fils['dats'], np.float32)

        dela = len(rawtmpdat)-RecPnts*2
        
        tmpdat = rawtmpdat[dela::]
        
        re = tmpdat[0::2]
        im = tmpdat[1::2]
        for m in range(0,len(re)):
            co[i,m] = complex(re[m], im[m])
        
        lb = np.exp(-np.linspace(start=0, stop=lineborad, num=len(re)))
        co3F[i,::] = co[i,::]*lb
        
        futF[i,::] = np.fft.fftshift(fft(co3F[i,::]))[::-1]
        
        futmagF[i,::] = np.sqrt(np.real(futF[i,::])**2 + np.imag(futF[i,::])**2)
        
        mxmag.append(np.max(futmagF[i,::]))
        
        
    # mxind = mxmag.index(max(mxmag))
    # mxpk = list(futmagF[mxind,::]).index(max(futmagF[mxind,::]))
    # crr = ppms_tmp[mxpk]-171
    ppms= ppms_tmp#-crr

    foldsDis['0']['0'] = lgtdt['0']
    datFut[str(0)][str(0)] = futF
    datFutmag[str(0)][str(0)] = futmagF
    datco3[str(0)][str(0)] = co3F
        
    
    return(ppms, datFutmag, datFut, datco3, foldsDis)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to change the ppm shift having as reference the most intense peak of all scans only for Pulsar benchtop NMR experiments
def shiftPPM(ppm_pre, datFutmag, machi = 'MRI'):
    if machi == 'Pulsar':
        LB = 2000
        while (LB > 1000 or LB < -1000):
            LB = simpledialog.askstring("PPM Shift", "At which ppm would you like the highest peak (overall) be?:")
            try:
                LB = float(LB)
            except:
                LB = -1
        mxmag = []
        for i in range(len(datFutmag['0']['0'])):
            mxmag.append(np.max(datFutmag['0']['0'][i]))
        
        mxind = mxmag.index(max(mxmag))
        mxpk = list(datFutmag['0']['0'][mxind]).index(max(datFutmag['0']['0'][mxind]))
        crr = ppm_pre[mxpk]-LB
        ppms= ppm_pre-crr        
                
    elif machi=='MRI':
        ppms = ppm_pre
        
        
    return(ppms)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------




def getMetaDat(foldPath, foldNum):
    rot = []
    with open(foldPath+'/'+foldNum+'/pdata/1/reco', 'r+') as f:
        for line in f:
            if 'RECO_rotate' in line.strip():
                # print(line.strip())
                # print(f.next())
                # print(next(f,''))
                rot = float(next(f,''))
                break

    RECO_ft_mode = []
    with open(foldPath+'/'+foldNum+'/pdata/1/reco', 'r+') as f:
        for line in f:
            if 'RECO_ft_mode' in line.strip():
                # print(line.strip())
                # print(f.next())
                # print(next(f,''))
                RECO_ft_mode = next(f,'')
                break     

    bw = []
    with open(foldPath+'/'+foldNum+'/method', 'r+') as f:
        for line in f:
            if 'PVM_SpecSW=' in line.strip():
                # print(line.strip())
                # print(f.next())
                # print(next(f,''))
                bw = float(next(f,''))
                break


    bwHz = []
    with open(foldPath+'/'+foldNum+'/method', 'r+') as f:
        for line in f:
            if 'PVM_SpecSWH=' in line.strip():
                # print(line.strip())
                # print(f.next())
                # print(next(f,''))
                bwHz = float(next(f,''))
                break


    bwc = []
    with open(foldPath+'/'+foldNum+'/method', 'r+') as f:
        for line in f:
            if 'PVM_FrqWorkPpm' in line.strip():
                # print(line.strip())
                # print(f.next())
                # print(np.array(next(f,'')).astype(np.float))
                bwc = next(f,'')
                break
    bwc = float(bwc[0:bwc.index(' ')])

    ACQ_repetition_time = []
    with open(foldPath+'/'+foldNum+'/acqp', 'r+') as f:
        for line in f:
            if '##$ACQ_repetition_time=' in line.strip():
                ACQ_repetition_time = float(next(f,''))/1000
                break

    ACQ_Size = []
    with open(foldPath+'/'+foldNum+'/method', 'r+') as f:
        for line in f:
            if '##$PVM_SpecMatrix=' in line.strip():
                ACQ_Size = int(next(f,''))
                break

    NR = []
    with open(foldPath+'/'+foldNum+'/acqp', 'r+') as f:
        for line in f:
            if '##$NR=' in line.strip():
                # NR = float(next(f,''))
                NR = int(line.strip()[6:])
                break

    ti = []
    with open(foldPath+'/'+str(foldNum)+'/acqp', 'r+') as f:
        for line in f:
            if '$$ Write Options' in line.strip():
                # print(line.strip())
                # print(f.next())
                # print(next(f,''))
                ti = next(f,'')
                break

    return(rot, RECO_ft_mode, bw, bwc, ACQ_repetition_time, ACQ_Size, NR, bwHz, ti)




# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Function to get the necessary parameters from the NMR acquisition done in the pulsar
def getMetaDat_Pulsar(foldPath, foldNum2):
    fls = os.listdir(foldPath+'/'+foldNum2+'/')
    fils = {}
    for i in range(len(fls)):
        if fls[i][-3::] == 'par':
            fils['pars'] = fls[i]
        elif fls[i][-3::] == 'fid':
            fils['dats'] = fls[i]
        
    RecPnts = []
    with open(foldPath+'/'+foldNum2+'/'+fils['pars'], 'r+') as f:
        for line in f:
            if 'ReceiverPoints' in line.strip():
                RecPnts = int(line.strip()[15::].split()[0])
                break
            
    SF = []
    with open(foldPath+'/'+foldNum2+'/'+fils['pars'], 'r+') as f:
        for line in f:
            if 'SF ' in line.strip():
                SF = float(line.strip()[3::].split()[0])
                break
            
    O1 = []
    with open(foldPath+'/'+foldNum2+'/'+fils['pars'], 'r+') as f:
        for line in f:
            if 'O1 ' in line.strip():
                O1 = float(line.strip()[3::].split()[0])
                break
            
    Filter = []
    with open(foldPath+'/'+foldNum2+'/'+fils['pars'], 'r+') as f:
        for line in f:
            if 'Filter ' in line.strip():
                Filter = int(line.strip()[6::].split()[0])
                break
    return(RecPnts, SF, O1, Filter)



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Functions for an interactive plot to visualise the magnitude data while eltting the user select the subexperiment, slab selected and scan and if they want to save the plots or not. 

def plotMagSp(ppms, datFutmag, foldPath, exp, wll, scn, sav=False, zoomlow = [], zoomlowY = 0):
    fig, ax = plt.subplots(figsize=(7, 4), dpi=170)
    
    yy = datFutmag[str(exp-1)][str(wll-1)][scn-1][abs(zoomlow[1]):abs(zoomlow[0])]
    
    ax.plot(ppms[abs(zoomlow[1]):abs(zoomlow[0])], yy, linewidth = 0.5, color='#80000aff')
    ax.invert_xaxis()
    ax.set_ylim(min(datFutmag[str(exp-1)][str(wll-1)][scn-1])-(max(datFutmag[str(exp-1)][str(wll-1)][scn-1])*0.1), max(datFutmag[str(exp-1)][str(wll-1)][scn-1])+(max(datFutmag[str(exp-1)][str(wll-1)][scn-1])*0.1))
    ax.set_ylabel("Intensity")
    ax.set_xlabel(r"$^{13}$C Chemical Shift (ppm)")
    ax.set_title("Magnitude Spectra Sub.Exp. "+str(exp)+", Slab "+str(wll)+", Scan "+str(scn))
    ax.spines[['right', 'top']].set_visible(False)
    
    # if zoomlowY <= 99:
    #     ax.set_ylim(min(yy)-(max(yy)*0.05), max(yy)+((max(yy)*0.05)-max(yy)*zoomlowY/100))
    # else:
    #     ax.set_ylim(min(yy)-(min(yy)*0.1), max(yy)+((max(yy)*0.05)-max(yy)*zoomlowY/150))
        
    mm = np.linspace(104, 150, 5100)
    mm1 = np.logspace(104, 104.999999, 5100)
    mm2 = 104 + ((mm1-min(mm1))*(104.999999 - 104))   /(max(mm1)-min(mm1))

    ts =(max(yy)+((max(yy)*0.05)-max(yy)*zoomlowY/100))*0.05
    if zoomlowY <= 99:
        ax.set_ylim(min(yy)-(ts), max(yy)+((max(yy)*0.05)-max(yy)*zoomlowY/100))
    elif zoomlowY >= 99 and zoomlowY <= 104:
        ax.set_ylim(min(yy)-(ts), max(yy)+((max(yy)*0.05)-max(yy)*zoomlowY/100))
    else:
        a = (mm >= zoomlowY)
        ind = np.where(a == True)[0][0]
        ts =(max(yy)+((max(yy)*0.05)-max(yy)*mm2[ind]/100))*0.05
        ax.set_ylim(min(yy)-(abs(ts)), max(yy)+((max(yy)*0.05)-max(yy)*mm2[ind]/100))
        
        
    
    if sav == True:
        plt.savefig(foldPath+'/Plots/MagnitudeSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.png')
        plt.savefig(foldPath+'/Plots/MagnitudeSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.svg')

    
def interMag(datFutmag, ppms, foldPath):    
    maxspc = max([len(datFutmag[str(m)][str(s)]) for m in range(len(datFutmag)) for s in range(len(datFutmag[str(m)]))])
    maxdta = max([len(datFutmag[str(m)][str(s)][0]) for m in range(len(datFutmag)) for s in range(len(datFutmag[str(m)]))])
    mm=interact(
        plotMagSp,
        exp = widgets.Dropdown(options=[list(range(len(datFutmag)))[g]+1 for g in range(len(datFutmag))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datFutmag["0"])))[g]+1 for g in range(len(datFutmag["0"]))], value = 1, description = "Slab:"),
        scn = widgets.Select(options=[list(range(maxspc))[g]+1 for g in range(maxspc)], value = 1, description = "Scan:"),
        ppms=fixed(ppms),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        zoomlow=widgets.IntRangeSlider(min=-maxdta,max=0,step=1,value=[-maxdta, 0], readout = False, description ='Zoom X', layout=Layout(width='1000px')),
        zoomlowY=widgets.IntSlider(min=0,max=150,step=0.01,value=0, readout = False, description ='Zoom Y', layout=Layout(width='1000px')),
        datFutmag=fixed(datFutmag),
        foldPath=fixed(foldPath))


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Functions for an interactive plot to visualise the furier transformed spectra (real or imaginary part) data while eltting the user select the subexperiment, slab selected and scan and if they want to save the plots or not. 

def plotFTSp(ppms, datFut, foldPath, exp, wll, scn, datyp = "real", ppmSelec = 0, sav=False, phsd=False, zoomlow = [], zoomlowY = 0):
    
    fig, ax = plt.subplots(figsize=(7, 4), dpi=170)
    if datyp == "imag":
        dtt = np.imag(datFut[str(exp-1)][str(wll-1)][scn-1][abs(zoomlow[1]):abs(zoomlow[0])])
    else:
        dtt = np.real(datFut[str(exp-1)][str(wll-1)][scn-1][abs(zoomlow[1]):abs(zoomlow[0])])

    ax.plot([ppmSelec,ppmSelec],[min(dtt), max(dtt)], color = "Gray", alpha = 0.5)
    ax.plot(ppms[abs(zoomlow[1]):abs(zoomlow[0])], dtt, linewidth = 0.5, color='#80000aff')
    
    ax.invert_xaxis()
    ax.set_ylim(min(dtt)-(max(dtt)*0.1), max(dtt)+(max(dtt)*0.1))
    ax.set_ylabel("Intensity")
    ax.set_xlabel(r"$^{13}$C Chemical Shift (ppm)")
    if datyp == "imag":
        if phsd == False:
            ax.set_title("Imaginary FFT Spectra Sub.Exp. "+str(exp)+", Slab "+str(wll)+", Scan "+str(scn)) 
        else:
            ax.set_title("Imaginary Phased FFT Spectra Sub.Exp. "+str(exp)+", Slab "+str(wll)+", Scan "+str(scn)) 
    else: 
        if phsd == False:
            ax.set_title("Real FFT Spectra Sub.Exp. "+str(exp)+", Slab "+str(wll)+", Scan "+str(scn)) 
        else:
            ax.set_title("Real Phased FFT Spectra Sub.Exp. "+str(exp)+", Slab "+str(wll)+", Scan "+str(scn))
            
    ax.spines[['right', 'top']].set_visible(False)
    
    
    
    mm = np.linspace(104, 150, 5100)
    mm1 = np.logspace(104, 104.999999, 5100)
    mm2 = 104 + ((mm1-min(mm1))*(104.999999 - 104))   /(max(mm1)-min(mm1))

    ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))*0.05
    if zoomlowY <= 99:
        ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    elif zoomlowY >= 99 and zoomlowY <= 104:
        ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    else:
        a = (mm >= zoomlowY)
        ind = np.where(a == True)[0][0]
        ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))*0.05
        ax.set_ylim(min(dtt)-(abs(ts)), max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))
        
    
    if sav == True:
        if datyp == "imag":
            if phsd == False:
                plt.savefig(foldPath+'/Plots/ImagFFTSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.png')
                plt.savefig(foldPath+'/Plots/ImagFFTSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.svg')
            else:
                plt.savefig(foldPath+'/Plots/ImagFFTSpecPhased_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.png')
                plt.savefig(foldPath+'/Plots/ImagFFTSpecPhased_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.svg')
        else:
            if phsd == False:
                plt.savefig(foldPath+'/Plots/RealFFTSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.png')
                plt.savefig(foldPath+'/Plots/RealFFTSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.svg')
            else:
                plt.savefig(foldPath+'/Plots/RealFFTSpecPhased_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.png')
                plt.savefig(foldPath+'/Plots/RealFFTSpecPhased_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.svg')
        


def interFTSpc(datFut, ppms, foldPath, datyp, phsd=False):    
    maxspc = max([len(datFut[str(m)][str(s)]) for m in range(len(datFut)) for s in range(len(datFut[str(m)]))])
    maxdta = max([len(datFut[str(m)][str(s)][0]) for m in range(len(datFut)) for s in range(len(datFut[str(m)]))])
    mm=interact(
        plotFTSp,
        exp = widgets.Dropdown(options=[list(range(len(datFut)))[g]+1 for g in range(len(datFut))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datFut["0"])))[g]+1 for g in range(len(datFut["0"]))], value = 1, description = "Slab:"),
        scn = widgets.Select(options=[list(range(maxspc))[g]+1 for g in range(maxspc)], value = 1, description = "Scan:"),
        ppms=fixed(ppms),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        datFut=fixed(datFut),
        datyp=fixed(datyp),
        phsd=fixed(phsd),
        ppmSelec = widgets.FloatSlider(min = min(ppms), max = max(ppms), value = max(ppms), layout=Layout(width='1400px'), step = 0.0001, description = "PPM:"),
        zoomlow=widgets.IntRangeSlider(min=-maxdta,max=0,step=1,value=[-maxdta, 0], readout = False, description ='Zoom X', layout=Layout(width='1000px')),
        zoomlowY=widgets.IntSlider(min=0,max=150,step=0.01,value=0, readout = False, description ='Zoom Y', layout=Layout(width='1000px')),
        foldPath=fixed(foldPath))


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Functions for an interactive plot to visualise the FID data while eltting the user select the subexperiment, slab selected and scan and if they want to save the plots or not. 

def plotFIDs(foldsDis, datco3, foldPath, exp, wll, scn, metDat, machi, sav=False):
    
    dtt = datco3[str(exp-1)][str(wll-1)][scn-1]
    
    if machi == 'MRI':
        ACQ_Time = []
        with open(foldPath+'/'+foldsDis[str(exp-1)][str(wll-1)][scn-1]+'/method', 'r+') as f:
            for line in f:
                if '##$PVM_SpecAcquisitionTime=' in line.strip():
                    ACQ_Time = float(line.strip()[27:])
                    break
        
        xax = np.linspace(0,ACQ_Time,len(dtt))
    elif machi == 'Pulsar':
        RecPnts, SF, O1, Filter = metDat['0']['raw1']
        freq = np.linspace(-Filter/2, Filter/2, RecPnts)
        xax = np.linspace(0, 1/(freq[1]-freq[0]), RecPnts)*1000

            
    # fig, (ax1, ax2) = plt.subplots(figsize=(7, 4), dpi=170)
    fig, (ax1, ax2) = plt.subplots(2,1, dpi=174, figsize=(7, 4))
    
    

    ax1.plot(xax, np.real(dtt), linewidth = 0.5, color='#80000aff', label='real')
    ax2.plot(xax, np.imag(dtt), linewidth = 0.5, color='#650080ff', label='imag')
    ax2.set_xlabel('time (ms)')
    ax2.set_ylabel('Intensity')
    ax1.set_ylabel('Intensity')
    ax1.legend()
    ax2.legend()
    ax1.set_title("FID") 
    ax1.spines[['right', 'top']].set_visible(False)
    ax2.spines[['right', 'top']].set_visible(False)
    # plt.show()
    
    if sav == True:
        plt.savefig(foldPath+'/Plots/FIDSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.png')
        plt.savefig(foldPath+'/Plots/FIDSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scan'+str(scn)+'.svg')

# plotFIDs(foldsDis, datco3, foldPath, exp, wll, scn, sav=False)

def interFIDs(datco3, foldsDis, foldPath, machi, metDat):    
    maxspc = max([len(datco3[str(m)][str(s)]) for m in range(len(datco3)) for s in range(len(datco3[str(m)]))])
    mm=interact(
        plotFIDs,
        exp = widgets.Dropdown(options=[list(range(len(datco3)))[g]+1 for g in range(len(datco3))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datco3["0"])))[g]+1 for g in range(len(datco3["0"]))], value = 1, description = "Slab:"),
        scn = widgets.Select(options=[list(range(maxspc))[g]+1 for g in range(maxspc)], value = 1, description = "Scan:"),
        foldsDis=fixed(foldsDis),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        datco3=fixed(datco3),
        machi=fixed(machi),
        metDat=fixed(metDat),
        foldPath=fixed(foldPath))
    
 
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
    
# Function for the selection of data to plot
def plotDatGen(ppms, foldsDis, datco3, datFutmag, datFut, foldPath, metDat=[], machi = [], datyp = "FID"):
    if datyp == "FID":
        interFIDs(datco3, foldsDis, foldPath, machi, metDat)
    elif datyp == "mag":
        interMag(datFutmag, ppms, foldPath)
    elif datyp == "real":
        interFTSpc(datFut, ppms, foldPath, datyp)
    elif datyp == "imag":
        interFTSpc(datFut, ppms, foldPath, datyp)
    elif datyp == "phsd":
        interFTSpc(datFut, ppms, foldPath, datyp, phsd=True)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# This are the necessary functions to interactively do phase correction for all sub experiments and slab selections. 

def phasecorr_Long(datFut, phcorr0, phcorr1, pivot, phasePars, foldPath, exp, wll, scn, zoomlowY, zoomlow = []):


        dataset = datFut[str(exp-1)][str(wll-1)][scn-1][abs(zoomlow[1]):abs(zoomlow[0])]
        
        fig, ax = plt.subplots(figsize=(10, 4), dpi=170)
        phaseddata = dataset * np.exp(
            1j * (phcorr0 + phcorr1 * (
                np.arange(-pivot, -pivot+dataset.size)/dataset.size)))

        ax.axvline(pivot, color='Gray', alpha=0.5, linewidth=2)
        ax.plot(np.real(phaseddata), color='#80000aff', linewidth = 0.5)

        # ax.set(ylim=(np.min(np.real(phaseddata))*2, np.max(np.real(phaseddata))*1.2))
        # ax.set_ylabel("Intensity")
        # ax.set_xlabel(r"$^{13}$C Chemical Shift (ppm)")
        ax.spines[['right', 'top']].set_visible(False)
        
        
        mm = np.linspace(104, 150, 5100)
        mm1 = np.logspace(104, 104.999999, 5100)
        mm2 = 104 + ((mm1-min(mm1))*(104.999999 - 104))   /(max(mm1)-min(mm1))
        dtt = np.real(phaseddata)
        
        ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))*0.05
        if zoomlowY <= 99:
            ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
        elif zoomlowY >= 99 and zoomlowY <= 104:
            ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
        else:
            a = (mm >= zoomlowY)
            ind = np.where(a == True)[0][0]
            ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))*0.05
            ax.set_ylim(min(dtt)-(abs(ts)), max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))
        
        
        # if zoomlowY <= 99:
        #     ax.set_ylim(min(np.real(phaseddata))-(max(np.real(phaseddata))*0.05), max(np.real(phaseddata))+((max(np.real(phaseddata))*0.05)-max(np.real(phaseddata))*zoomlowY/100))
        # else:
        #     ax.set_ylim(min(np.real(phaseddata))-(min(np.real(phaseddata))*0.1), max(np.real(phaseddata))+((max(np.real(phaseddata))*0.05)-max(np.real(phaseddata))*zoomlowY/150))
        
        
        plt.gca().invert_xaxis()

        plt.show()

        
        phasePars[str(exp-1)][str(wll-1)] = [phcorr0, phcorr1, pivot]
                
        # np.save(foldPath+'/PhasePars.npy', phasePars)

        with open(foldPath+'/Results/PhasePars.p', 'wb') as fp:
                pickle.dump(phasePars, fp, protocol=pickle.HIGHEST_PROTOCOL)

# phasecorr_Long(ppms, datFut, 1, 1, 1, foldPath, 1, 1, 1)

def interactPhase(foldPath, datFut):
        phagan = True
        if os.path.isfile(foldPath+'/Results/PhasePars.p'):
                root = tk.Tk()
                root.attributes('-topmost',True)
                root.iconify()
                phagan = tk.messagebox.askyesno("WARNING", "You already phased this data. If you proceed the phase values will be deleted and you will have to re-do it. Are you sure you want to proceed?")
                root.destroy()
                
        if phagan == True:
                phasePars={}
                for i in range(len(datFut)):
                        phasePars[str(i)]={}
                        for j in range(len(datFut[str(i)])):
                                phasePars[str(i)][str(j)] = [0,0,0]
                phcorr0a, phcorr1a, pivota = 0,0,0
        
                maxspc = max([len(datFut[str(m)][str(s)]) for m in range(len(datFut)) for s in range(len(datFut[str(m)]))])
                maxdta = max([len(datFut[str(m)][str(s)][0]) for m in range(len(datFut)) for s in range(len(datFut[str(m)]))])

                mm=interact(
                        phasecorr_Long,
                        datFut=fixed(datFut),
                        phasePars=fixed(phasePars),
                        exp = widgets.Dropdown(options=[list(range(len(datFut)))[g]+1 for g in range(len(datFut))], value = 1, description = "Sub.Exp.:"),
                        wll = widgets.Dropdown(options=[list(range(len(datFut["0"])))[g]+1 for g in range(len(datFut["0"]))], value = 1, description = "Slab:"),
                        scn = widgets.Select(options=[list(range(maxspc))[g]+1 for g in range(maxspc)], value = 1, description = "Scan:"),
                        phcorr0=widgets.FloatSlider(min=-np.pi*2, max=np.pi*2, step=0.0001, value=phcorr0a, layout=Layout(width='1000px')),
                        phcorr1=widgets.FloatSlider(min=-30*np.pi, max=30*np.pi, step=0.0001, value=phcorr1a, layout=Layout(width='1000px')),
                        pivot=widgets.IntSlider(min=0, max=maxdta, step=0.001, value=pivota, layout=Layout(width='1000px')),
                        zoomlowY=widgets.IntSlider(min=0,max=150,step=0.01,value=0, readout = False, description ='Zoom Y', layout=Layout(width='1000px')),
                        zoomlow=widgets.IntRangeSlider(min=-maxdta,max=0,step=1,value=[-maxdta, 0], readout = False, description ='Zoom', layout=Layout(width='1000px')),
                        foldPath=fixed(foldPath))





# This function applyes the phase correction (and baseline correction) to your dataset iteratively
def appPhsBsl(foldPath, datFut, ppms, rdbc = False): # Set last parameter to true to re-do the baseline correction
    with open(foldPath+'/Results/PhasePars.p', 'rb') as fp:
        phasePars = pickle.load(fp)
        
    if os.path.isfile(foldPath+'/Results/datBaselines.p') and rdbc == False:
        with open(foldPath+'/Results/datFutPhBSl.p', 'rb') as fp:
            datFutPhBSl = pickle.load(fp)
        with open(foldPath+'/Results/datFutPh.p', 'rb') as fp:
            datFutPh = pickle.load(fp)
        with open(foldPath+'/Results/datBaselines.p', 'rb') as fp:
            datBaselines = pickle.load(fp)
    else:    
        datFutPh = {}
        datFutPhBSl = {}
        datBaselines = {}
        for i in range(len(datFut)):
            datFutPh[str(i)] = {}
            datFutPhBSl[str(i)] = {}
            datBaselines[str(i)] = {}   
            for j in range(len(datFut[str(i)])):        
                datFutPh[str(i)][str(j)] = np.zeros(np.shape(datFut[str(i)][str(j)]), dtype=np.complex128)
                datFutPhBSl[str(i)][str(j)] = np.zeros(np.shape(datFut[str(i)][str(j)]))
                datBaselines[str(i)][str(j)] = np.zeros(np.shape(datFut[str(i)][str(j)]))
                phcorr0, phcorr1, pivot = phasePars[str(i)][str(j)]
                for m in range(len(datFut[str(i)][str(j)])):
                    phaseddata = datFut[str(i)][str(j)][m] * np.exp(
                        1j * (phcorr0 + phcorr1 * (
                            np.arange(-int(pivot), -int(pivot)+datFut[str(i)][str(j)][m].size)/datFut[str(i)][str(j)][m].size)))
                    
                    baseline_fitter = pybaselines.Baseline(x_data=ppms) # Baseline correction with pybaseline
                    stiff_baseline = baseline_fitter.arpls(phaseddata)[0]
                    datFutPhBSl[str(i)][str(j)][m] = np.real(phaseddata)-stiff_baseline#-baseline(np.real(phaseddata), deg = 70) # This includes baseline correction!
                    datFutPh[str(i)][str(j)][m] = phaseddata
                    datBaselines[str(i)][str(j)][m] = stiff_baseline
                    
        with open(foldPath+'/Results/datFutPhBSl.p', 'wb') as fp:
                    pickle.dump(datFutPhBSl, fp, protocol=pickle.HIGHEST_PROTOCOL)  
        with open(foldPath+'/Results/datFutPh.p', 'wb') as fp:
                    pickle.dump(datFutPh, fp, protocol=pickle.HIGHEST_PROTOCOL)
        with open(foldPath+'/Results/datBaselines.p', 'wb') as fp:
                    pickle.dump(datBaselines, fp, protocol=pickle.HIGHEST_PROTOCOL)
    
        
    return(phasePars, datFutPh, datFutPhBSl, datBaselines)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Functions to plot a stack of the spectra acquired in a 3D plot for better visualisation. Users can move the plot around, change colors and all the functions from the previous plots. 

def plot3DStack(ppms, datFutPhBSl, exp, wll, foldPath, scns=1, col = "#80000aff", el=25, az=290, grad = False, sav=[], zoomlow = [-1, 0], zoomlowY = 0):
    
    
    dtt = datFutPhBSl[str(exp-1)][str(wll-1)][0:scns]
    
    mxl = [max(dtt[i]) for i in range(len(dtt))]
    mxm = max([max(dtt[i]) for i in range(len(dtt))])
    mxi = mxl.index(mxm)

    
    fig = plt.figure(figsize=(10,15))
    ax = plt.axes(projection='3d')
    xline = ppms[abs(zoomlow[1]):abs(zoomlow[0])]
    
    viridis = mpl.colormaps["viridis"].resampled(len(dtt))
    for i in range(len(dtt)):
        dtt_tmp = dtt[len(dtt)-1-i]
        if grad == False:
            ax.plot3D(xline, np.ones(len(xline))+(len(dtt)-1-i), dtt_tmp[abs(zoomlow[1]):abs(zoomlow[0])], color=col)
        else:
            ax.plot3D(xline, np.ones(len(xline))+(len(dtt)-1-i), dtt_tmp[abs(zoomlow[1]):abs(zoomlow[0])], color=viridis.colors[len(dtt)-1-i])
    ax.view_init(el, az)
    ax.invert_xaxis()
    
    ax.set_xlabel(r"$^{13}$C Chemical Shift (ppm)", fontsize=15, rotation=-10)
    ax.set_ylabel("Scan", fontsize=15, rotation=60)
    ax.set_zlabel("Intensity", fontsize=15, rotation=-10)
    
    
    
    # if zoomlowY <= 99:
    #     ax.set_zlim(min(dtt[mxi])-(max(dtt[mxi])*0.05), max(dtt[mxi])+((max(dtt[mxi])*0.05)-max(dtt[mxi])*zoomlowY/100))
    # else:
    #     ax.set_zlim(min(dtt[mxi])-(min(dtt[mxi])*0.1), max(dtt[mxi])+((max(dtt[mxi])*0.05)-max(dtt[mxi])*zoomlowY/150))
        
        
    mm = np.linspace(104, 150, 5100)
    mm1 = np.logspace(104, 104.999999, 5100)
    mm2 = 104 + ((mm1-min(mm1))*(104.999999 - 104))   /(max(mm1)-min(mm1))

    ts =(max(dtt[mxi])+((max(dtt[mxi])*0.05)-max(dtt[mxi])*zoomlowY/100))*0.05
    if zoomlowY <= 99:
        ax.set_zlim(min(dtt[mxi])-(ts), max(dtt[mxi])+((max(dtt[mxi])*0.05)-max(dtt[mxi])*zoomlowY/100))
    elif zoomlowY >= 99 and zoomlowY <= 104:
        ax.set_zlim(min(dtt[mxi])-(ts), max(dtt[mxi])+((max(dtt[mxi])*0.05)-max(dtt[mxi])*zoomlowY/100))
    else:
        a = (mm >= zoomlowY)
        ind = np.where(a == True)[0][0]
        ts =(max(dtt[mxi])+((max(dtt[mxi])*0.05)-max(dtt[mxi])*mm2[ind]/100))*0.05
        ax.set_zlim(min(dtt[mxi])-(abs(ts)), max(dtt[mxi])+((max(dtt[mxi])*0.05)-max(dtt[mxi])*mm2[ind]/100))
        
    
    
    if sav == True:
        plt.savefig(foldPath+'/Plots/Stack3D_SubExp'+str(exp)+'_Slab'+str(wll)+'.png')
        plt.savefig(foldPath+'/Plots/Stack3D_SubExp'+str(exp)+'_Slab'+str(wll)+'.svg')
    
    
def interStack3D(foldPath, datFutPhBSl, ppms, clls = ["#80000aff"], grad=False):
    
    maxspc = max([len(datFutPhBSl[str(m)][str(s)]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    maxdta = max([len(datFutPhBSl[str(m)][str(s)][0]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    mm=interact(
        plot3DStack,
        col = widgets.Dropdown(options=clls, value = "#80000aff", description = "Colour:"),
        grad = fixed(grad),
        az = widgets.FloatSlider(min=0,max=360,step=0.01,value=290, description ='azim', layout=Layout(width='1000px')),
        el = widgets.FloatSlider(min=0,max=100,step=0.01,value=25, description ='elev', layout=Layout(width='1000px')),
        exp = widgets.Dropdown(options=[list(range(len(datFutPhBSl)))[g]+1 for g in range(len(datFutPhBSl))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datFutPhBSl["0"])))[g]+1 for g in range(len(datFutPhBSl["0"]))], value = 1, description = "Slab:"),
        scns = widgets.IntSlider(min=1,max=maxspc,step=1,value=maxspc, description ='Scans:', layout=Layout(width='1000px')),
        ppms=fixed(ppms),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        zoomlow=widgets.IntRangeSlider(min=-maxdta,max=0,step=1,value=[-maxdta, 0], readout = False, description ='Zoom X', layout=Layout(width='1000px')),
        zoomlowY=widgets.IntSlider(min=0,max=150,step=0.01,value=0, readout = False, description ='Zoom Y', layout=Layout(width='1000px')),
        datFutPhBSl=fixed(datFutPhBSl),
        foldPath=fixed(foldPath)) 



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Functions to plot a stack of the spectra acquired in a 2D plot (MNova Style) for better visualisation. Users can move the plot around, change colors and all the functions from the previous plots. 


def plotStack2D(foldPath, datFutPhBSl, ppms, exp=1, wll=1, scns=1, col = "#80000aff", multip = 1, zoomlowY=0, sav=False, zoomlow = [-1, 0]):
    
    dtt = datFutPhBSl[str(exp-1)][str(wll-1)][0:scns]
    nrmfc_up = max([max(dtt[i]) for i in range(len(dtt))])
    ss = [max(dtt[i]) for i in range(len(dtt))].index(nrmfc_up)
    nrmfc_dw = min(dtt[ss])

    dtt_nrm = [(dtt[i]-nrmfc_dw)/(nrmfc_up-nrmfc_dw) for i in range(len(dtt))]
    
    tmprs = (100-zoomlowY)/100
    for i in range(len(dtt_nrm)):
        for j in range(len(dtt_nrm[i])):
            if dtt_nrm[i][j] > tmprs:
                dtt_nrm[i][j] = np.NaN
            
                
    fig, ax = plt.subplots(figsize=(10, 4), dpi=170)
    for i in range(len(dtt)):
        ax.plot(ppms[abs(zoomlow[1]):abs(zoomlow[0])], dtt_nrm[i][abs(zoomlow[1]):abs(zoomlow[0])]*multip+i+1, color=col, linewidth = 0.5)

    ax.set_ylabel("Scan")
    ax.set_xlabel(r"$^{13}$C Chemical Shift (ppm)")
    ax.spines[['right', 'top']].set_visible(False)
    
    plt.gca().invert_xaxis()
    
    
    
    if sav == True:
        plt.savefig(foldPath+'/Plots/Stack2D_SubExp'+str(exp)+'_Slab'+str(wll)+'.png')
        plt.savefig(foldPath+'/Plots/Stack2D_SubExp'+str(exp)+'_Slab'+str(wll)+'.svg')

    plt.show()
    
# plotStack2D(foldPath, datFutPhBSl, ppms, exp=1, wll=1, scns=3, zoomlow=[-1024,0], multip=1)    


def interStack2D(foldPath, datFutPhBSl, ppms, clls = ["#80000aff"]):
    
    maxspc = max([len(datFutPhBSl[str(m)][str(s)]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    maxdta = max([len(datFutPhBSl[str(m)][str(s)][0]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    mm=interact(
        plotStack2D,
        col = widgets.Dropdown(options=clls, value = "#80000aff", description = "Colour:"),
        exp = widgets.Dropdown(options=[list(range(len(datFutPhBSl)))[g]+1 for g in range(len(datFutPhBSl))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datFutPhBSl["0"])))[g]+1 for g in range(len(datFutPhBSl["0"]))], value = 1, description = "Slab:"),
        scns = widgets.IntSlider(min=1,max=maxspc,step=1,value=maxspc, description ='Scans:', layout=Layout(width='1000px')),
        ppms=fixed(ppms),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        zoomlow=widgets.IntRangeSlider(min=-maxdta,max=0,step=1,value=[-maxdta, 0], readout = False, description ='Zoom X', layout=Layout(width='1000px')),
        multip=widgets.FloatSlider(min=1,max=200,step=0.01,value=0, readout = False, description ='Zoom Y', layout=Layout(width='1000px')),
        zoomlowY=widgets.FloatSlider(min=0,max=99.99,step=0.0001,value=0, readout = False, description ='Zoom Y2', layout=Layout(width='1000px')),
        datFutPhBSl=fixed(datFutPhBSl),
        foldPath=fixed(foldPath)) 


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to ask the user how many Spectra wants for the sum up. It is the same experiment, so for reproducibility it should be the same number for all
def getSumUpNum(maxSpex):

    SUS = -1
    while (SUS > maxSpex or SUS < 1):
        SUS = simpledialog.askstring("Sum Up", "Enter Number of Spectra you want for the Sum Up (maximum is "+str(maxSpex)+"):")
        try:
            SUS = int(SUS)
        except:
            SUS = -1
            
    return(SUS)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Function to generate the sum up spectra, asking the user how many spectra to take for it. 
def genSumUpDat(foldPath, datFutPhBSl):
    
    maxSpcnm = max([len(datFutPhBSl[str(i)][str(j)]) for i in range(len(datFutPhBSl)) for j in range(len(datFutPhBSl[str(i)]))])
    SUS = getSumUpNum(maxSpcnm)

    datFutPhBSl_Sum = {}
    for i in range(len(datFutPhBSl)):
        datFutPhBSl_Sum[str(i)] = {}
        for j in range(len(datFutPhBSl[str(i)])):
            datFutPhBSl_Sum[str(i)][str(j)] = sum(datFutPhBSl[str(i)][str(j)][0:SUS])

    with open(foldPath+'/Results/datFutPhBSl_Sum.p', 'wb') as fp:
            pickle.dump(datFutPhBSl_Sum, fp, protocol=pickle.HIGHEST_PROTOCOL)  
            
    return(datFutPhBSl_Sum, SUS)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Functions to visualise the sum up spectra in an interactive plot
def plotSumUpSp(foldPath, ppms, datFutPhBSl_Sum, exp, wll, sav=False, zoomlow = [], zoomlowY = 0, col = "#80000aff"):
    
    fig, ax = plt.subplots(figsize=(7, 4), dpi=170)

    dtt = datFutPhBSl_Sum[str(exp-1)][str(wll-1)][abs(zoomlow[1]):abs(zoomlow[0])]

    ax.plot(ppms[abs(zoomlow[1]):abs(zoomlow[0])], dtt, linewidth = 0.5, color=col)
    
    ax.invert_xaxis()
    ax.set_ylim(min(dtt)-(max(dtt)*0.1), max(dtt)+(max(dtt)*0.1))
    ax.set_ylabel("Intensity")
    ax.set_xlabel(r"$^{13}$C Chemical Shift (ppm)")

    ax.set_title("Sum Up Sub.Exp. "+str(exp)+", Slab "+str(wll)) 
    
    ax.spines[['right', 'top']].set_visible(False)
    
    
    
    # if zoomlowY <= 99:
    #     ax.set_ylim(min(dtt)-(max(dtt)*0.05), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    # else:
    #     ax.set_ylim(min(dtt)-(min(dtt)*0.1), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/150))
        
    mm = np.linspace(104, 150, 5100)
    mm1 = np.logspace(104, 104.999999, 5100)
    mm2 = 104 + ((mm1-min(mm1))*(104.999999 - 104))   /(max(mm1)-min(mm1))

    ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))*0.05
    if zoomlowY <= 99:
        ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    elif zoomlowY >= 99 and zoomlowY <= 104:
        ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    else:
        a = (mm >= zoomlowY)
        ind = np.where(a == True)[0][0]
        ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))*0.05
        ax.set_ylim(min(dtt)-(abs(ts)), max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))
        
    
    if sav == True:
        plt.savefig(foldPath+'/Plots/SumUpSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'.png')
        plt.savefig(foldPath+'/Plots/SumUpSpec_SubExp'+str(exp)+'_Slab'+str(wll)+'.svg')
            
            

def interSumUpSp(foldPath, datFutPhBSl_Sum, ppms, clls = ["#db6030ff"]):
    
    maxspc = max([len(datFutPhBSl_Sum[str(m)][str(s)]) for m in range(len(datFutPhBSl_Sum)) for s in range(len(datFutPhBSl_Sum[str(m)]))])
    maxdta = max([len(datFutPhBSl_Sum[str(m)][str(s)]) for m in range(len(datFutPhBSl_Sum)) for s in range(len(datFutPhBSl_Sum[str(m)]))])
    mm=interact(
        plotSumUpSp,
        col = widgets.Dropdown(options=clls, value = "#db6030ff", description = "Colour:"),
        exp = widgets.Dropdown(options=[list(range(len(datFutPhBSl_Sum)))[g]+1 for g in range(len(datFutPhBSl_Sum))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datFutPhBSl_Sum["0"])))[g]+1 for g in range(len(datFutPhBSl_Sum["0"]))], value = 1, description = "Slab:"),
        ppms=fixed(ppms),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        zoomlow=widgets.IntRangeSlider(min=-maxdta,max=0,step=1,value=[-maxdta, 0], readout = False, description ='Zoom X', layout=Layout(width='1000px')),
        zoomlowY=widgets.IntSlider(min=0,max=150,step=0.01,value=0, readout = False, description ='Zoom Y', layout=Layout(width='1000px')),
        datFutPhBSl_Sum=fixed(datFutPhBSl_Sum),
        foldPath=fixed(foldPath)) 
    



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Functions to display spectra to estimate the width of any desire peak at any point
def plotBWCheck(ppms, datFutPhBSl, datFutPhBSl_Sum, metDat, exp, wll, scn, BW, machi = [],  zoomlow = [-1,0], zoomlowY = 0):
    
    fig, ax = plt.subplots(figsize=(7, 4), dpi=170)

    if scn == "Sum Up":
        dtt = datFutPhBSl_Sum[str(exp-1)][str(wll-1)][abs(zoomlow[1]):abs(zoomlow[0])]
    else:
        dtt = datFutPhBSl[str(exp-1)][str(wll-1)][scn-1][abs(zoomlow[1]):abs(zoomlow[0])]

    ax.plot(ppms[abs(zoomlow[1]):abs(zoomlow[0])], dtt, linewidth = 0.5, color="#80000aff")
    
    ax.invert_xaxis()
    ax.set_ylim(min(dtt)-(max(dtt)*0.1), max(dtt)+(max(dtt)*0.1))
    ax.set_ylabel("Intensity")
    ax.set_xlabel(r"$^{13}$C Chemical Shift (ppm)")

    drs = list(metDat[str(exp-1)].keys())
    ppmBW = str(abs(BW[0]) - abs(BW[1]))
    if machi == 'MRI':
        hzBW = str(round((metDat[str(exp-1)][str(drs[wll-1])][7]/(max(ppms)-min(ppms)))*(abs(BW[0]) - abs(BW[1])), 4))
    elif machi == 'Pulsar':
        RecPnts, SF, O1, Filter = metDat['0']['raw1']        
        freq1 = (BW[0]*SF)+O1
        freq2 = (BW[1]*SF)+O1
        
        hzBW = str(freq2-freq1)
     

    ax.set_title("Bandwidth: "+hzBW+" Hz, "+ppmBW+" ppm") 
    
    ax.spines[['right', 'top']].set_visible(False)
    
    
    
    # if zoomlowY <= 99:
    #     ax.set_ylim(min(dtt)-(max(dtt)*0.05), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    # else:
    #     ax.set_ylim(min(dtt)-(min(dtt)*0.1), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/150))
        
    mm = np.linspace(104, 150, 5100)
    mm1 = np.logspace(104, 104.999999, 5100)
    mm2 = 104 + ((mm1-min(mm1))*(104.999999 - 104))   /(max(mm1)-min(mm1))

    ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))*0.05
    if zoomlowY <= 99:
        ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    elif zoomlowY >= 99 and zoomlowY <= 104:
        ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    else:
        a = (mm >= zoomlowY)
        ind = np.where(a == True)[0][0]
        ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))*0.05
        ax.set_ylim(min(dtt)-(abs(ts)), max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))
    
    
    ax.plot([abs(BW[1]), abs(BW[1])],[0, np.max(dtt)], color = 'black', alpha = 0.4)
    ax.plot([abs(BW[0]), abs(BW[0])],[0, np.max(dtt)], color = 'black', alpha = 0.4)
    ax.fill_between(ppms, 0, np.max(dtt), where= (ppms > abs(BW[1])) & (ppms < abs(BW[0])),
                    facecolor='red', alpha=0.2)
    
    
    


def getFWHM(foldPath, datFutPhBSl, datFutPhBSl_Sum, ppms, metDat, machi = []):
    
    maxspc = max([len(datFutPhBSl[str(m)][str(s)]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    maxdta = max([len(datFutPhBSl[str(m)][str(s)][0]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    lstsh = ["Sum Up"] + [list(range(maxspc))[g]+1 for g in range(maxspc)]
        
    mm=interact(
        plotBWCheck,
        exp = widgets.Dropdown(options=[list(range(len(datFutPhBSl_Sum)))[g]+1 for g in range(len(datFutPhBSl_Sum))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datFutPhBSl_Sum["0"])))[g]+1 for g in range(len(datFutPhBSl_Sum["0"]))], value = 1, description = "Slab:"),
        scn =  widgets.Select(options=lstsh, value = "Sum Up", description = "Scan:"),
        ppms=fixed(ppms),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        zoomlow=widgets.IntRangeSlider(min=-maxdta,max=0,step=1,value=[-maxdta, 0], readout = False, description ='Zoom X', layout=Layout(width='1000px')),
        zoomlowY=widgets.IntSlider(min=0,max=150,step=0.01,value=0, readout = False, description ='Zoom Y', layout=Layout(width='1000px')),
        datFutPhBSl_Sum=fixed(datFutPhBSl_Sum),
        datFutPhBSl = fixed(datFutPhBSl),
        metDat = fixed(metDat),
        BW = widgets.FloatRangeSlider(min=-np.max(ppms),max=-np.min(ppms),step=0.000001,value=[-round(ppms[round(len(ppms)/2)]), -round(ppms[round(len(ppms)/2)])+1], readout = False, layout=Layout(width='1400px')),
        machi = fixed(machi),
        foldPath=fixed(foldPath)) 



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to extract important info about the processing to remind the user the next time it tries to process it again. It also asks the user if they want to add any additional comment for the future. Careful since every time you run this function, the previous file will be overwritten. 
def getDatFuture(lngrps, numWells, foldsDis, lineborad, SUS, foldPath):
    lines = []
    lines.append('The following directory contains '+str(lngrps)+' subexperiments. \n')
    lines.append('Each sub experiment contains '+str(numWells)+' sequential Slabs. \n')
    ss = [len(foldsDis[str(i)]["0"]) for i in range(len(foldsDis))]
    for j in range(lngrps):
        lines.append('The subexperiment '+str(j+1)+' has '+str(ss[j])+' scans per slab. \n')
        
    lines.append('In the last processing, you selected an apodisation factor of '+str(lineborad)+'. \n')
    lines.append('In the last processing, you selected '+str(SUS)+' scans for the sum up final spectra. \n')
    lines.append('\n \n')
    lines.append('User Comments: \n')
    
    Coms = simpledialog.askstring("Comments", "Add any comment or important information you want to be reminded of next time you load this directory:")
    lines.append(Coms)
    
    with open(foldPath+'/Results/ExpComments.txt', 'w') as f:
        for line in lines:
            f.write(line)
        f.close()



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to extract and print the info from getDatFuture. 
def printInfo(foldPath):
    if os.path.isfile(foldPath+'/Results/ExpComments.txt'):
        infl = []
        f = open(foldPath+'/Results/ExpComments.txt', "r")
        for x in f:
            infl.append(x)
        f.close()
            
        print("EXPERIMENT PROCESSED INFORMATION:")
        print("---------------------------------")
        print("\n")
        for i in infl:
            print(i)
    else:
        print("First Time You Process This Directory")


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Functions to display spectra to estimate the width of any desire peak at any point
def plotSNRCheck(ppms, datFutPhBSl, datFutPhBSl_Sum, exp, wll, scn, Noise=[0,1], Peak=[0,1], zoomlow = [-1,0], zoomlowY = 0):
    
    fig, ax = plt.subplots(figsize=(7, 4), dpi=170)

    if scn == "Sum Up":
        dtt = datFutPhBSl_Sum[str(exp-1)][str(wll-1)][abs(zoomlow[1]):abs(zoomlow[0])]
    else:
        dtt = datFutPhBSl[str(exp-1)][str(wll-1)][scn-1][abs(zoomlow[1]):abs(zoomlow[0])]

    ax.plot(ppms[abs(zoomlow[1]):abs(zoomlow[0])], dtt, linewidth = 0.5, color="#80000aff")
    
    ax.invert_xaxis()
    ax.set_ylim(min(dtt)-(max(dtt)*0.1), max(dtt)+(max(dtt)*0.1))
    ax.set_ylabel("Intensity")
    ax.set_xlabel(r"$^{13}$C Chemical Shift (ppm)")

    
    ax.spines[['right', 'top']].set_visible(False)
    
    
    
    # if zoomlowY <= 99:
    #     ax.set_ylim(min(dtt)-(max(dtt)*0.05), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    # else:
    #     ax.set_ylim(min(dtt)-(min(dtt)*0.1), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/150))
        
    mm = np.linspace(104, 150, 5100)
    mm1 = np.logspace(104, 104.999999, 5100)
    mm2 = 104 + ((mm1-min(mm1))*(104.999999 - 104))   /(max(mm1)-min(mm1))

    ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))*0.05
    if zoomlowY <= 99:
        ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    elif zoomlowY >= 99 and zoomlowY <= 104:
        ax.set_ylim(min(dtt)-(ts), max(dtt)+((max(dtt)*0.05)-max(dtt)*zoomlowY/100))
    else:
        a = (mm >= zoomlowY)
        ind = np.where(a == True)[0][0]
        ts =(max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))*0.05
        ax.set_ylim(min(dtt)-(abs(ts)), max(dtt)+((max(dtt)*0.05)-max(dtt)*mm2[ind]/100))
    
    
    ax.plot([abs(Noise[1]), abs(Noise[1])],[0, np.max(dtt)], color = 'black', alpha = 0.4)
    ax.plot([abs(Noise[0]), abs(Noise[0])],[0, np.max(dtt)], color = 'black', alpha = 0.4)

    ax.plot([abs(Peak[1]), abs(Peak[1])],[0, np.max(dtt)], color = 'black', alpha = 0.4)
    ax.plot([abs(Peak[0]), abs(Peak[0])],[0, np.max(dtt)], color = 'black', alpha = 0.4)
    
    ax.fill_between(ppms, 0, np.max(dtt), where= (ppms > abs(Noise[1])) & (ppms < abs(Noise[0])),
                    facecolor='red', alpha=0.2)
    ax.fill_between(ppms, 0, np.max(dtt), where= (ppms > abs(Peak[1])) & (ppms < abs(Peak[0])),
                    facecolor='green', alpha=0.2)
    
    
    pH = np.where(ppms == ppms[ppms > abs(Peak[0])][0])[0][0]
    pL = np.where(ppms == ppms[ppms > abs(Peak[1])][0])[0][0]

    nH = np.where(ppms == ppms[ppms > abs(Noise[0])][0])[0][0]
    nL = np.where(ppms == ppms[ppms > abs(Noise[1])][0])[0][0]
    SNR = np.max(dtt[pL:pH])/np.std(dtt[nL:nH])
    
    ax.set_title("SNR: "+str(round(SNR, 4)))    


def getSNR(foldPath, datFutPhBSl, datFutPhBSl_Sum, ppms):
    
    maxspc = max([len(datFutPhBSl[str(m)][str(s)]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    maxdta = max([len(datFutPhBSl[str(m)][str(s)][0]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    lstsh = ["Sum Up"] + [list(range(maxspc))[g]+1 for g in range(maxspc)]
        
    mm=interact(
        plotSNRCheck,
        exp = widgets.Dropdown(options=[list(range(len(datFutPhBSl_Sum)))[g]+1 for g in range(len(datFutPhBSl_Sum))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datFutPhBSl_Sum["0"])))[g]+1 for g in range(len(datFutPhBSl_Sum["0"]))], value = 1, description = "Slab:"),
        scn =  widgets.Select(options=lstsh, value = "Sum Up", description = "Scan:"),
        ppms=fixed(ppms),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        zoomlow=widgets.IntRangeSlider(min=-maxdta,max=0,step=1,value=[-maxdta, 0], readout = False, description ='Zoom X', layout=Layout(width='1000px')),
        zoomlowY=widgets.IntSlider(min=0,max=150,step=0.01,value=0, readout = False, description ='Zoom Y', layout=Layout(width='1000px')),
        datFutPhBSl_Sum=fixed(datFutPhBSl_Sum),
        datFutPhBSl = fixed(datFutPhBSl),
        Noise = widgets.FloatRangeSlider(min=-np.max(ppms),max=-np.min(ppms),step=0.01,value=[-np.max(ppms)+5, -np.max(ppms)+6], readout = False, layout=Layout(width='1400px')),
        Peak = widgets.FloatRangeSlider(min=-np.max(ppms),max=-np.min(ppms),step=0.01,value=[-round(ppms[round(len(ppms)/2)]), -round(ppms[round(len(ppms)/2)])+1], readout = False, layout=Layout(width='1400px')),
        foldPath=fixed(foldPath)) 


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# THis functions are used to select the integration regions (upt to 7 for now) for all experiments and subexperiments

def plotIntegrat(foldPath, ppms, datFutPhBSl, intPars, exp, wll, scns, multip = 1, zoomlowY=0 , zoomlow = [-1,0], int1=[0,0], int2=[0,0], int3=[0,0], int4=[0,0], int5=[0,0], int6=[0,0], int7=[0,0]):
    
    intAll = [int1, int2, int3, int4, int5, int6, int7]
    
    fig, ax = plt.subplots(figsize=(7, 4), dpi=170)

    dtt = datFutPhBSl[str(exp-1)][str(wll-1)][0:scns]
    nrmfc_up = max([max(dtt[i]) for i in range(len(dtt))])
    ss = [max(dtt[i]) for i in range(len(dtt))].index(nrmfc_up)
    nrmfc_dw = min(dtt[ss])

    dtt_nrm = [(dtt[i]-nrmfc_dw)/(nrmfc_up-nrmfc_dw) for i in range(len(dtt))]
    
    tmprs = (100-zoomlowY)/100
    for i in range(len(dtt_nrm)):
        for j in range(len(dtt_nrm[i])):
            if dtt_nrm[i][j] > tmprs:
                dtt_nrm[i][j] = np.NaN

    # ax.plot(ppms[abs(zoomlow[1]):abs(zoomlow[0])], dtt, linewidth = 0.5, color="#80000aff")
    
    mxvls = []
    for i in range(len(dtt)):
        ax.plot(ppms[abs(zoomlow[1]):abs(zoomlow[0])], dtt_nrm[i][abs(zoomlow[1]):abs(zoomlow[0])]*multip+i+1, color="#80000aff", linewidth = 0.5)
        mxvls.append(max(dtt_nrm[i][abs(zoomlow[1]):abs(zoomlow[0])]*multip+i+1))
    
    ax.invert_xaxis()
    ax.set_ylabel("Scan")
    ax.set_xlabel(r"$^{13}$C Chemical Shift (ppm)")

    
    ax.spines[['right', 'top']].set_visible(False)
        
    # cols = ['#F94144','#F3722C','#F8961E','#F9C74F','#90BE6D','#43AA8B','#577590']
    cols = ['#bc0a00ff','#bc9c00ff','#84bc00ff','#00adbcff','#2700bcff','#bc00b3ff','#000000ff']
    
    for i in range(len(intAll)):
        if intAll[i] != [0,0]:
            ax.plot([abs(intAll[i][1]), abs(intAll[i][1])],[1, max(mxvls)], color = 'black', alpha = 0.4)
            ax.plot([abs(intAll[i][0]), abs(intAll[i][0])],[1, max(mxvls)], color = 'black', alpha = 0.4)
            ax.fill_between(ppms, 1, max(mxvls), where= (ppms > abs(intAll[i][1])) & (ppms < abs(intAll[i][0])),
                            facecolor=cols[i], alpha=0.2)

            intPars[str(exp-1)][str(wll-1)][i] = intAll[i]
            
    with open(foldPath+'/Results/IntegrPars.p', 'wb') as fp:
        pickle.dump(intPars, fp, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    
    
    
def interactIntegrals(foldPath, datFutPhBSl, lngrps, numWells, ppms):
    phagan = True
    if os.path.isfile(foldPath+'/Results/IntegrPars.p'):
        root = tk.Tk()
        root.attributes('-topmost',True)
        root.iconify()
        phagan = tk.messagebox.askyesno("WARNING", "You already selected the integration regions for this data. If you proceed the these values will be deleted and you will have to re-do it. Are you sure you want to proceed?")
        root.destroy()
            
    if phagan == True:
        intPars = {}
        for i in range(lngrps):
            intPars[str(i)] = {}
            for j in range(numWells):
                intPars[str(i)][str(j)] = np.zeros((7,2))
    
        maxspc = max([len(datFutPhBSl[str(m)][str(s)]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
        maxdta = max([len(datFutPhBSl[str(m)][str(s)][0]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
        mm=interact(
            plotIntegrat,
            exp = widgets.Dropdown(options=[list(range(len(datFutPhBSl)))[g]+1 for g in range(len(datFutPhBSl))], value = 1, description = "Sub.Exp.:"),
            wll = widgets.Dropdown(options=[list(range(len(datFutPhBSl["0"])))[g]+1 for g in range(len(datFutPhBSl["0"]))], value = 1, description = "Slab:"),
            scns = widgets.IntSlider(min=1,max=maxspc,step=1,value=maxspc, description ='Scans:', layout=Layout(width='1000px')),
            ppms=fixed(ppms),
            sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
            zoomlow=widgets.IntRangeSlider(min=-maxdta,max=0,step=1,value=[-maxdta, 0], readout = False, description ='Zoom X', layout=Layout(width='1000px')),
            multip=widgets.FloatSlider(min=1,max=100,step=0.01,value=0, readout = False, description ='Zoom Y', layout=Layout(width='1000px')),
            zoomlowY=widgets.FloatSlider(min=0,max=99.99,step=0.0001,value=0, readout = False, description ='Zoom Y2', layout=Layout(width='1000px')),
            datFutPhBSl=fixed(datFutPhBSl),
            intPars = fixed(intPars),
            int1 = widgets.FloatRangeSlider(min=-max(ppms),max=-min(ppms),step=0.001,value=[0, 0], readout = False, layout=Layout(width='1000px'), description ='Integral 1:'),
            int2 = widgets.FloatRangeSlider(min=-max(ppms),max=-min(ppms),step=0.001,value=[0, 0], readout = False, layout=Layout(width='1000px'), description ='Integral 2:'),
            int3 = widgets.FloatRangeSlider(min=-max(ppms),max=-min(ppms),step=0.001,value=[0, 0], readout = False, layout=Layout(width='1000px'), description ='Integral 3:'),
            int4 = widgets.FloatRangeSlider(min=-max(ppms),max=-min(ppms),step=0.001,value=[0, 0], readout = False, layout=Layout(width='1000px'), description ='Integral 4:'),
            int5 = widgets.FloatRangeSlider(min=-max(ppms),max=-min(ppms),step=0.001,value=[0, 0], readout = False, layout=Layout(width='1000px'), description ='Integral 5:'),
            int6 = widgets.FloatRangeSlider(min=-max(ppms),max=-min(ppms),step=0.001,value=[0, 0], readout = False, layout=Layout(width='1000px'), description ='Integral 6:'),
            int7 = widgets.FloatRangeSlider(min=-max(ppms),max=-min(ppms),step=0.001,value=[0, 0], readout = False, layout=Layout(width='1000px'), description ='Integral 7:'),
            foldPath=fixed(foldPath)) 
    

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to extract the integration regions assigned before and compute the AUC for each. 

def getIntegrals(foldPath, datFutPhBSl, ppms):
    
    with open(foldPath+'/Results/IntegrPars.p', 'rb') as fp:
        intPars = pickle.load(fp)

    integVals = {}    
    integPPMS = {}
    integInte = {}
    for i in range(len(intPars)): # Number of Subexperiments
        integVals[str(i)] = {}
        integPPMS[str(i)] = {}
        integInte[str(i)] = {}
        for j in range(len(intPars[str(i)])): # Number of Slabs
            integVals[str(i)][str(j)] = np.zeros((7, len(datFutPhBSl["0"]["0"])))
            integPPMS[str(i)][str(j)] = {}
            integInte[str(i)][str(j)] = {}
            for k in range(len(intPars[str(i)][str(j)])): # Number of integration regions (7 for now)
                if intPars[str(i)][str(j)][k][0] != intPars[str(i)][str(j)][k][1]:
                    ppmSe = ppms[(ppms > abs(intPars[str(i)][str(j)][k][1])) & (ppms < abs(intPars[str(i)][str(j)][k][0]))]
                    integPPMS[str(i)][str(j)][str(k)] = ppmSe
                    integInte[str(i)][str(j)][str(k)] = {}
                    for m in range(len(datFutPhBSl[str(i)][str(j)])): # Number of scans
                        intSe = datFutPhBSl[str(i)][str(j)][m][(ppms > abs(intPars[str(i)][str(j)][k][1])) & (ppms < abs(intPars[str(i)][str(j)][k][0]))]
                        integInte[str(i)][str(j)][str(k)][str(m)] = intSe
                        integVals[str(i)][str(j)][k,m] = auc(ppmSe, intSe)
            
            pd.DataFrame(np.transpose(integVals[str(i)][str(j)]), columns = ["Integral_1","Integral_2","Integral_3","Integral_4","Integral_5","Integral_6","Integral_7"]).to_csv(foldPath+'/Results/Integrals_SubExp'+str(i+1)+'_Slab'+str(j+1)+'.csv')             
                        
                        
    with open(foldPath+'/Results/integVals.p', 'wb') as fp:
            pickle.dump(integVals, fp, protocol=pickle.HIGHEST_PROTOCOL)
    with open(foldPath+'/Results/integPPMS.p', 'wb') as fp:
            pickle.dump(integPPMS, fp, protocol=pickle.HIGHEST_PROTOCOL)
    with open(foldPath+'/Results/integInte.p', 'wb') as fp:
            pickle.dump(integInte, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    return(integVals, integPPMS, integInte, intPars)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to interactively display the integration resutls for all subexperiments

def plotIntegrat2(foldPath, timVecs, integVals, exp, wll, scns, sec=True, nmsIn = ['','','','','','',''], multip = [1,1,1,1,1,1,1], sav=False):
    
    cols = ['#bc0a00ff','#bc9c00ff','#84bc00ff','#00adbcff','#2700bcff','#bc00b3ff','#000000ff']
    
    lgn = nmsIn.split('\n')
    
    if type(multip) == str:
        mpf = [float(multip.split('\n')[s]) for s in range(len(multip.split('\n')))]
    else:
        mpf = multip
        
    
    fig, ax = plt.subplots(figsize=(7, 4), dpi=170)

    ax.set_ylabel("AUC Intensity (A.U.)")
    if sec == True:
        ax.set_xlabel("time (s)")
        xx = 1
    else:
        ax.set_xlabel("time (min)")
        xx = 60
        
    for i in range(len(integVals[str(exp-1)][str(wll-1)])):
        if sum(integVals[str(exp-1)][str(wll-1)][i]) != 0:
            xd = [timVecs[str(exp-1)][str(wll-1)][p]/xx for p in range(len(timVecs[str(exp-1)][str(wll-1)]))][0:scns]
            yd = integVals[str(exp-1)][str(wll-1)][i][0:scns]*mpf[i]
            if mpf[i] ==1:
                if len(yd) == 1:
                    ax.scatter(xd, yd, label = lgn[i], color = cols[i], linewidth=4)
                else:
                    ax.plot(xd, yd, label = lgn[i], color = cols[i], linewidth=4)
            else:
                if len(yd) == 1:
                    ax.scatter(xd, yd, label = lgn[i]+"*"+str(mpf[i]), color = cols[i], linewidth=4)
                else:
                    ax.plot(xd, yd, label = lgn[i]+"*"+str(mpf[i]), color = cols[i], linewidth=4)
    
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend()
    
    if sav == True:
        plt.savefig(foldPath+'/Plots/Integrals_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scans1-'+str(scns)+'.png')
        plt.savefig(foldPath+'/Plots/Integrals_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scans1-'+str(scns)+'.svg')
        


def interIntegrals(foldPath, datco3, datFutPhBSl, integVals, timVecs):    
    maxspc = max([len(datFutPhBSl[str(m)][str(s)]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    mm=interact(
        plotIntegrat2,
        exp = widgets.Dropdown(options=[list(range(len(datco3)))[g]+1 for g in range(len(datco3))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datco3["0"])))[g]+1 for g in range(len(datco3["0"]))], value = 1, description = "Slab:"),
        scns = widgets.IntSlider(min=1,max=maxspc,step=1,value=maxspc, description ='Scans:', layout=Layout(width='1000px')),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        integVals=fixed(integVals),
        timVecs = fixed(timVecs),
        nmsIn = widgets.Textarea(value='Pyr\n\n\n\n\n\n', placeholder='Type something', description='Legend:', disabled=False),
        multip = widgets.Textarea(value='1\n1\n1\n1\n1\n1\n1', placeholder='Type something', description='Multip:', disabled=False),
        sec = widgets.Checkbox(value = True, description = "Seconds/Minutes"),
        foldPath=fixed(foldPath))
    
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Functions to interactively visualise AUC of data and ratios
def plotIntegratRatios(foldPath, timVecs, integVals, exp, wll, SUS, nmsIn = ['','','','','','',''], SumUp = True, rati1 = 0, rati2 = 0, sav=False):
    
    rati = [rati1, rati2]
    
    cols = ['#bc0a00ff','#bc9c00ff','#84bc00ff','#00adbcff','#2700bcff','#bc00b3ff','#000000ff']
    
    lgn = nmsIn.split('\n')
    
    fig, ax = plt.subplots(figsize=(7, 4), dpi=170)

    ax.set_ylabel("AUC Intensity (A.U.)")
    ax.set_xlabel("Metabolite")
   
        
    tmp = 0    
    if rati == [0,0]:
        
        for i in range(len(integVals[str(exp-1)][str(wll-1)])):
            if sum(integVals[str(exp-1)][str(wll-1)][i]) != 0:
                if SumUp == True:
                    yd = [sum(integVals[str(exp-1)][str(wll-1)][i][0:SUS])]
                    xd = [1]
                    ax.set_title("Sum Up Sub.Exp. "+str(exp)+", Slab "+str(wll))
                else:
                    ax.set_title("AUC of selected first "+str(SUS)+" scans Sub.Exp. "+str(exp)+", Slab "+str(wll))
                    if SUS > 1:
                        yd_pre = integVals[str(exp-1)][str(wll-1)][i][0:SUS]
                        xd_pre = [timVecs[str(exp-1)][str(wll-1)][p] for p in range(len(timVecs[str(exp-1)][str(wll-1)]))][0:SUS]
                        yd = auc(xd_pre, yd_pre)
                    else:
                        yd = integVals[str(exp-1)][str(wll-1)][i][0:SUS]


                ax.bar([i], yd, label = lgn[i], color = cols[i], linewidth=4, alpha = 0.2)
                ax.bar([i], yd, linewidth=4, alpha = 1, edgecolor = cols[i], color = 'None')
                
                print(yd)
                
                tmp += 1
        ax.set_xticks(range(tmp), range(1, tmp+1))
        ax.legend()
                
    else: 
        try: 
            if SumUp == True:
                rtss = np.zeros(len(integVals[str(exp-1)]))
                for s in range(len(rtss)):
                    yd = sum(integVals[str(exp-1)][str(s)][rati[0]-1][0:SUS]) / sum(integVals[str(exp-1)][str(s)][rati[1]-1][0:SUS])
                    xd = [1]
                    rtss[s] = yd
                ax.set_title("Sum Up Sub.Exp. "+str(exp))
                ax.set_xlabel("Slab")
            else:
                rtss = np.zeros(len(integVals[str(exp-1)]))
                ax.set_title("AUC Ratio of selected first "+str(SUS)+" scans Sub.Exp. "+str(exp))
                
                for s in range(len(rtss)):
                    yd_pre1 = integVals[str(exp-1)][str(s)][rati[0]-1][0:SUS]
                    yd_pre2 = integVals[str(exp-1)][str(s)][rati[1]-1][0:SUS]
                    xd_pre = [timVecs[str(exp-1)][str(s)][p] for p in range(len(timVecs[str(exp-1)][str(wll-1)]))][0:SUS]
                    yd = auc(xd_pre, yd_pre1)/auc(xd_pre, yd_pre2)
                    rtss[s] = yd
                
            print(rtss)
                
            ax.bar(range(len(rtss)), rtss, label = lgn[1], color = cols[0:len(rtss)], linewidth=4, alpha = 0.2)
            ax.bar(range(len(rtss)), rtss, linewidth=4, alpha = 1, edgecolor = cols[0:len(rtss)], color = 'None')
            
            ax.set_xticks(range(len(rtss)), range(1, len(rtss)+1))
            
            ax.set_ylabel("AUC "+lgn[rati[0]-1]+" / AUC "+lgn[rati[1]-1])
            # ax.legend()
            
        except:
            print("index selected for metabolite is empty, you did not define an integration region! :(")        
            
    
    ax.spines[['right', 'top']].set_visible(False)
    
    if sav == True:
        if rati == [0,0]:
            plt.savefig(foldPath+'/Plots/IntegralRatios_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scans1-'+str(SUS)+'.png')
            plt.savefig(foldPath+'/Plots/IntegralRatios_SubExp'+str(exp)+'_Slab'+str(wll)+'_Scans1-'+str(SUS)+'.svg')
        else:
            plt.savefig(foldPath+'/Plots/IntegralRatios_SubExp'+str(exp)+'_Scans1-'+str(SUS)+'.png')
            plt.savefig(foldPath+'/Plots/IntegralRatios_SubExp'+str(exp)+'_Scans1-'+str(SUS)+'.svg')
        


def interIntegralsRatios(foldPath, datco3, datFutPhBSl, integVals, timVecs, SUS):    
    maxspc = max([len(datFutPhBSl[str(m)][str(s)]) for m in range(len(datFutPhBSl)) for s in range(len(datFutPhBSl[str(m)]))])
    mm=interact(
        plotIntegratRatios,
        exp = widgets.Dropdown(options=[list(range(len(datco3)))[g]+1 for g in range(len(datco3))], value = 1, description = "Sub.Exp.:"),
        wll = widgets.Dropdown(options=[list(range(len(datco3["0"])))[g]+1 for g in range(len(datco3["0"]))], value = 1, description = "Slab:"),
        SUS = fixed(SUS),
        sav = widgets.Checkbox(value = False, description = "Save Plot ON"),
        integVals=fixed(integVals),
        timVecs = fixed(timVecs),
        SumUp = widgets.Checkbox(value = False, description = "Use Sum Up: "),
        rati1 = widgets.Dropdown(options=list(range(0,7)), value = 0, description = "Top Ratio: "),
        rati2 = widgets.Dropdown(options=list(range(0,7)), value = 0, description = "Bottom Ratio: "),
        nmsIn = widgets.Textarea(value='Pyr\n\n\n\n\n\n', placeholder='Type something', description='Legend:', disabled=False),
        foldPath=fixed(foldPath))

    
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------






















































