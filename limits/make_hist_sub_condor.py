import os

categories=[
    "mumu_OS", 
    "mumu_SS",
    "ee_OS",
    "ee_SS",
    "mue_OS",
    "mue_SS",
    "emu_OS",
    "emu_SS",
    ]

ncats = len(categories)

regions = ["A", "B", "C", "D"]
nregions = len(regions)

procsbkg = []
procsHNL = []
procsData = []
years = ["2016","2017","2018"]#["2016","2018"]#, "2017", "2018"]

import yaml
with open("../config/samples.yml") as samples_file:
    samples_dict = yaml.load(samples_file, Loader=yaml.FullLoader)
'''
samples_dict = [
    'HNL_dirac_ntau_ctau1p0e-01_massHNL3p0_Vall1p388e-01',
    'HNL_dirac_ntau_ctau1p0e-01_massHNL4p5_Vall4p549e-02',
    'HNL_dirac_ntau_ctau1p0e-01_massHNL6p0_Vall2p054e-02',
    'HNL_dirac_ntau_ctau1p0e-02_massHNL4p5_Vall1p438e-01',
    'HNL_dirac_ntau_ctau1p0e-02_massHNL6p0_Vall6p496e-02',
    'HNL_dirac_ntau_ctau1p0e-02_massHNL8p0_Vall2p996e-02',
    'HNL_dirac_ntau_ctau1p0e-03_massHNL10p0_Vall5p262e-02',
    'HNL_dirac_ntau_ctau1p0e-03_massHNL6p0_Vall2p054e-01',
    'HNL_dirac_ntau_ctau1p0e-03_massHNL8p0_Vall9p475e-02',
    'HNL_dirac_ntau_ctau1p0e-04_massHNL10p0_Vall1p664e-01',
    'HNL_dirac_ntau_ctau1p0e-04_massHNL12p0_Vall1p035e-01',
    'HNL_dirac_ntau_ctau1p0e-04_massHNL8p0_Vall2p996e-01',
    'HNL_dirac_ntau_ctau1p0e-05_massHNL12p0_Vall3p272e-01',
    'HNL_dirac_ntau_ctau1p0e-05_massHNL14p0_Vall2p193e-01',
    'HNL_dirac_ntau_ctau1p0e-05_massHNL16p0_Vall1p551e-01',
    'HNL_dirac_ntau_ctau1p0e-05_massHNL18p0_Vall1p144e-01',
    'HNL_dirac_ntau_ctau1p0e-05_massHNL20p0_Vall8p709e-02',
    'HNL_dirac_ntau_ctau1p0e00_massHNL2p0_Vall1p286e-01',
    'HNL_dirac_ntau_ctau1p0e00_massHNL3p0_Vall4p390e-02',

    'HNL_majorana_ntau_ctau1p0e-01_massHNL3p0_Vall9p825e-02',
    'HNL_majorana_ntau_ctau1p0e-01_massHNL4p5_Vall3p213e-02',
    'HNL_majorana_ntau_ctau1p0e-01_massHNL6p0_Vall1p454e-02',
    'HNL_majorana_ntau_ctau1p0e-02_massHNL4p5_Vall1p016e-01',
    'HNL_majorana_ntau_ctau1p0e-02_massHNL6p0_Vall4p597e-02',
    'HNL_majorana_ntau_ctau1p0e-02_massHNL8p0_Vall2p119e-02',
    'HNL_majorana_ntau_ctau1p0e-03_massHNL10p0_Vall3p721e-02',
    'HNL_majorana_ntau_ctau1p0e-03_massHNL6p0_Vall1p454e-01',
    'HNL_majorana_ntau_ctau1p0e-03_massHNL8p0_Vall6p702e-02',
    'HNL_majorana_ntau_ctau1p0e-04_massHNL10p0_Vall1p177e-01',
    'HNL_majorana_ntau_ctau1p0e-04_massHNL12p0_Vall7p319e-02',
    'HNL_majorana_ntau_ctau1p0e-04_massHNL8p0_Vall2p119e-01',
    'HNL_majorana_ntau_ctau1p0e-05_massHNL12p0_Vall2p314e-01',
    'HNL_majorana_ntau_ctau1p0e-05_massHNL14p0_Vall1p551e-01',
    'HNL_majorana_ntau_ctau1p0e-05_massHNL16p0_Vall1p097e-01',
    'HNL_majorana_ntau_ctau1p0e-05_massHNL18p0_Vall8p092e-02',
    'HNL_majorana_ntau_ctau1p0e-05_massHNL20p0_Vall6p165e-02',
    'HNL_majorana_ntau_ctau1p0e00_massHNL2p0_Vall9p078e-02',
    'HNL_majorana_ntau_ctau1p0e00_massHNL3p0_Vall3p107e-02'
]
'''



for l in samples_dict:
    
    
    if ("HNL" in l):
    
        if l.find("_all_")>=0:
            if (l.replace("_all_","_pt20_") in samples_dict.keys()) or (l.replace("_all_","_ntau_") in samples_dict.keys()):
                #print ("skip",l)
                continue
            else:
                print ("keep",l)
        '''
        if l.find("_pt20_")>=0:
            if (l.replace("_pt20_","_ntau_") in samples_dict.keys()):
                #print ("skip",l)
                continue
            else:
                print ("keep",l)
        '''
        '''
        if l.find("HNL_dirac_all_ctau1p0e00_massHNL10p0")>=0:
            procsHNL.append(l)
        if l.find("HNL_dirac_pt20_ctau1p0e00_massHNL10p0")>=0:
            procsHNL.append(l)
        if l.find("HNL_dirac_ntau_ctau1p0e00_massHNL10p0")>=0:
            procsHNL.append(l)
        
        if l.find("HNL_majorana_ntau_ctau1p0e-03_massHNL14p0")>=0:
            procsHNL.append(l)
        if l.find("HNL_dirac_ntau_ctau1p0e-03_massHNL14p0")>=0:
            procsHNL.append(l)
        '''
        procsHNL.append(l)
            
    elif l == "muon" or l == "electron": 
        procsData.append(l)
    else:
        continue
        #procsbkg.append(l)

nprocsbkg = len(procsbkg)
nprocsHNL = len(procsHNL)
nprocsData = len(procsData)

header = """Executable   = exec_histo.sh
Output       = %s/log.$(Process).out
Error        = %s/log.$(Process).err
Log          = /dev/null

Request_CPUs   = 1
Request_Memory = 1.95GB
+RequestRuntime = 10440
Should_Transfer_Files = NO

Requirements = ( OpSysAndVer == "CentOS7" && Machine != "batch1055.desy.de" && Machine != "batch1053.desy.de" )

queue arguments from (
"""
dagJobs = []
for year in years:
    jobs = []

    '''
    for proc in procsbkg:
        for region in regions:
            for category in categories:
                jobs.append(f'--proc {proc} --category {category} --region {region} --year {year}')
    '''
    
    for couplingRange in [range(1,14),range(14,27),range(27,41),range(41,54),range(54,68)]:
        optCoupling = ""
        for coupling in couplingRange:
            optCoupling += "--coupling %i "%coupling
        for proc in procsHNL:
            for region in regions : 
                for category in categories:
                    jobs.append(f'--proc {proc} --category {category} --region {region} --year {year} {optCoupling} --suffix _coupl'+str(couplingRange[0]))
    
    for proc in procsData:
        for region in regions:
            for category in categories:
                if proc == "muon" and category.startswith("mu"):
                    jobs.append(f'--data --proc {proc} --category {category} --region {region} --year {year}')
                elif proc == "electron" and category.startswith("e"):
                    jobs.append(f'--data --proc {proc} --category {category} --region {region} --year {year}')
    
    njobs = len(jobs)

    print(f"The number of jobs is {njobs}.")

    chunk_size = 4500
    for i,ichunk in enumerate(range(0, len(jobs), chunk_size)): 
        logFolder = f"log/{year}/chunk{i}"
        try:
            os.makedirs(logFolder)
        except:
            pass
    
        jobsPerChunk = jobs[ichunk:ichunk+chunk_size]

        dagJobs.append(f"subHists_{year}_{i}.condor")
        with open(f"subHists_{year}_{i}.condor", "w") as f:
            print("   - make job with "+str(len(jobsPerChunk)))
            f.write(header%(logFolder,logFolder))
            for job in jobsPerChunk:
                f.write('"'+job+'"\n')
            
            f.write(")\n")
            
with open("dag.condor","w") as f:
    for i,dagJob in enumerate(dagJobs):
        f.write(f"Job J{i+1} {dagJob}\n")
    for i in range(len(dagJobs)-1):
        f.write(f"PARENT J{i+1} CHILD J{i+2}\n")
    for i in range(len(dagJobs)):
        f.write(f"RETRY J{i+1} 2\n")    
            
