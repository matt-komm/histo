from hepdata_lib import *
from hepdata_lib.c_file_reader import CFileReader
import numpy as np
import scipy.interpolate
import tarfile

submission = Submission()
#TODO
submission.read_abstract("abstract.txt")
#TODO
submission.add_link("Webpage with all figures and tables", "https://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-21-013/")
#TODO
submission.add_link("arXiv", "http://arxiv.org/abs/arXiv:XXXXXXX")
#TODO
submission.add_record_id(0, "inspire")

itable = 1

def findCoupling2AtMass(mass,massCoupling2):
    
    crossings = []
    eps=1e-4
    for i in range(massCoupling2.shape[0]-1):
        if (massCoupling2[i,0]<(mass+eps) and massCoupling2[i+1,0]>(mass-eps)) or (massCoupling2[i,0]>(mass-eps) and massCoupling2[i+1,0]<(mass+eps)):
            coupl = np.interp(mass, massCoupling2[i:i+2,0], massCoupling2[i:i+2,1])
            #ctauAtCoupling2 = 10**np.poly1d(np.polyfit(np.log10(massCoupling2Ctau[selectMass][:,1]),np.log10(massCoupling2Ctau[selectMass][:,2]),1))(np.log10(coupling2AtMass))
            if len(crossings)==0:
                crossings.append(coupl)
            elif abs(crossings[-1]-coupl)/coupl>0.2:
                crossings.append(coupl)
    return crossings
    
def getCoupling2Values(arrMass,massCoupling2):
    
    arrCoupling2Crossings = -1*np.ones((len(arrMass),3),dtype=np.float32)
    for i in range(len(arrMass)):
        crossings = findCoupling2AtMass(arrMass[i],massCoupling2)
        if len(crossings)==0:
            pass
        elif len(crossings)==1:
            arrCoupling2Crossings[i,0]=crossings[0]    
        elif len(crossings)==2:
            arrCoupling2Crossings[i,0]=crossings[0]
            arrCoupling2Crossings[i,2]=crossings[1]
        elif len(crossings)==3:
            arrCoupling2Crossings[i,0]=crossings[0]
            arrCoupling2Crossings[i,1]=crossings[1]
            arrCoupling2Crossings[i,2]=crossings[2]
        elif len(crossings)==4 or len(crossings)==5:
            for j in range(3):
                idx = np.argmin(np.abs(arrCoupling2Crossings[i-1,j]-crossings))
                if abs(arrCoupling2Crossings[i-1,j]-crossings[idx])/arrCoupling2Crossings[i-1,j]<4.:
                    arrCoupling2Crossings[i,j] = crossings[idx]
                else:
                    print ("skip: ",arrMass[i],arrCoupling2Crossings[i-1,j],crossings)
            
        else:
            raise Exception("Unknown number of crossings: "+str(len(crossings)),str(crossings))
    return arrCoupling2Crossings
    
def buildMassCouplingTable(arrMass,dataName,descr="",loc="",coupling=""):
    global itable
    
    data = np.load(dataName)
    
    massCoupling2ToCtau = scipy.interpolate.interp2d(
        np.log(data['massCoupling2Ctau'][:,0]), 
        np.log(data['massCoupling2Ctau'][:,1]), 
        np.log(data['massCoupling2Ctau'][:,2]),
        kind='linear',
        bounds_error=True
    )
    def getCtau(mass,coupling2):
        selNeg = coupling2>0.
        return np.exp(massCoupling2ToCtau(
            np.log(mass),
            np.log(coupling2) 
        ))

    obsCoupling2 = getCoupling2Values(arrMass,data['obsMassCoupling2'])
    expCoupling2 = getCoupling2Values(arrMass,data['expMassCoupling2'])
    exp68UpCoupling2 = getCoupling2Values(arrMass,data['exp68UpMassCoupling2'])
    exp68DownCoupling2 = getCoupling2Values(arrMass,data['exp68DownMassCoupling2'])
    exp95UpCoupling2 = getCoupling2Values(arrMass,data['exp95UpMassCoupling2'])
    exp95DownCoupling2 = getCoupling2Values(arrMass,data['exp95DownMassCoupling2'])
   
   
    
    
    def addVar(name,label,arr,dep=False,units=""):
        var = Variable(name, is_independent=dep, is_binned=False, units=units)
        var.values = arr
        if not dep:
            var.add_qualifier("SQRT(S)", 13, "TeV")
            var.add_qualifier("LUMINOSITY", 138, "fb$^{-1}$")
            var.add_qualifier("QUANTITY", coupling)
            var.add_qualifier("Limit", label)
        return var
        
    def addPart(table,coupling2,name):
        ipart = 1
        if (coupling2[:,1]<0.).all() and (coupling2[:,2]<0.).all():
            table.add_variable(addVar(name,name,coupling2[:,0]))
        else:
            for i in range(0,3):
                if (coupling2[:,i]>0.).any():
                    table.add_variable(addVar(name+" (part "+str(ipart)+")",name,coupling2[:,i]))
                    ipart+=1


    table = Table("Table "+str(itable))
    itable+=1
    table.description = descr
    table.location = loc
    table.add_variable(addVar("HNL mass","HNL mass",arrMass,dep=True,units="GeV"))
    '''
    print (dataName)
    for i in range(len(arrMass)):
        print (arrMass[i],obsCoupling2[i,0],getCtau(arrMass[i],obsCoupling2[i,0]))
    '''
    
    addPart(table,obsCoupling2,"Observed")
    addPart(table,expCoupling2,"Expected")
    addPart(table,exp68UpCoupling2,"Expected 68% Up")
    addPart(table,exp68DownCoupling2,"Expected 68% Down")
    addPart(table,exp95UpCoupling2,"Expected 95% Up")
    addPart(table,exp95DownCoupling2,"Expected 95% Down")
        

    table.add_image(dataName.replace("_limit_new.npz","_xsec.pdf"))
        
    return table

    
def buildTriangleCouplingTable(dataName,descr="",loc="",quantity=""):
    global itable

    data = np.load(dataName)
    obsTrig = data['trig']
    
    
    def addVar(name,label,arr,dep=False,units=""):
        var = Variable(name, is_independent=dep, is_binned=False, units=units)
        var.values = arr
        if not dep:
            var.add_qualifier("SQRT(S)", 13, "TeV")
            var.add_qualifier("LUMINOSITY", 138, "fb$^{-1}$")
            var.add_qualifier("QUANTITY", quantity)
            var.add_qualifier("UNIT", units)
            var.add_qualifier("Limit", label)
        elif units!="":
            var.add_qualifier("UNIT", units)
        return var
        


    table = Table("Table "+str(itable))
    itable+=1
    table.description = descr
    table.location = loc
    
    #print (obsTrig)
    
    
    table.add_variable(addVar("Relative electron coupling","$|V_{\\text{e}\\text{N}}|^{2}/\\sum_{\\ell}|V_{\\ell\\text{N}}|^{2}$",obsTrig[:,0],dep=True,units=""))
    table.add_variable(addVar("Relative muon coupling","$|V_{\\mu\\text{N}}|^{2}/\\sum_{\\ell}|V_{\\ell\\text{N}}|^{2}$",obsTrig[:,1],dep=True,units=""))
    table.add_variable(addVar("Relative tau lepton coupling","$|V_{\\tau\\text{N}}|^{2}/\\sum_{\\ell}|V_{\\ell\\text{N}}|^{2}$",obsTrig[:,2],dep=True,units=""))
    
    if quantity=="mass":
        table.add_variable(addVar("Observed","Observed",obsTrig[:,3],units="GeV"))
    else:
        table.add_variable(addVar("Observed","Observed",obsTrig[:,3],units="mm"))
    
    table.add_image(dataName.replace(".npz",".pdf"))
        
    return table
    
    



def buildHistTable(filePath,descr="",loc=""):
    global itable
    
    reader = RootFileReader(filePath)
    bkg = reader.read_hist_1d("bkg")
    bkgErr = reader.read_hist_1d("err")
    data = reader.read_hist_1d("data")
    signal4p5 = reader.read_hist_1d("signal4p5")
    signal10p0 = reader.read_hist_1d("signal10p0")


    def addVar(name,arr,unc=None,uncName="",dep=False):
        var = Variable(name, is_independent=dep, is_binned=False)
        var.values = arr
        if unc is not None:
            uncVar = Uncertainty(uncName, is_symmetric=True)
            uncVar.values = unc
            var.add_uncertainty(uncVar)
        if not dep:
            var.add_qualifier("SQRT(S)", 13, "TeV")
            var.add_qualifier("LUMINOSITY", 138, "fb$^{-1}$")
        return var
        
    table = Table("Table "+str(itable))
    itable+=1
    table.description = descr
    table.location = loc

    signCatVar = Variable("Lepton charge combination", is_independent=True, is_binned=False, units="")
    signCatVar.values = ["OS"]*12+["SS"]*12
    table.add_variable(signCatVar)

    displCatVar = Variable("Second lepton displacement", is_independent=True, is_binned=False, units="")
    displCatVar.values = (["$d_{xy}^\\text{sig}<3$"]*4+["$3<d_{xy}^\\text{sig}<10$"]*4+["$d_{xy}^\\text{sig}>10$"]*4)*2
    table.add_variable(displCatVar)

    flavCatVar = Variable("Lepton flavour combination", is_independent=True, is_binned=False, units="")
    flavCatVar.values = ["$\\mu\\mu$","$\\text{ee}$","$\\mu\\text{e}$","$\\text{e}\\mu$"]*6
    table.add_variable(flavCatVar)


    table.add_variable(addVar(
        "Number of background events", bkg["y"],
        unc=bkgErr["dy"], uncName="Uncertainty"
    ))
    
    table.add_variable(addVar(
        "Data", data["y"],
        unc=data["dy"], uncName="Poisson errors"
    ))
    
    table.add_variable(addVar(
        "Signal scenario 1 ($m=4.5~\\text{GeV}$, $c\\tau=100~\\text{mm}$)",signal4p5["y"],
        unc=signal4p5["dy"],  uncName="Uncertainty"
    ))
    table.add_variable(addVar(
        "Signal scenario 2 ($m=10~\\text{GeV}$, $c\\tau=1~\\text{mm}$)",signal10p0["y"],
        unc=signal10p0["dy"],  uncName="Uncertainty"
    ))
    
    return table

tableResolved = buildHistTable(
    "resolved.root",
    descr="Observed number of events and predicted number of background events per category for resolved categories",
    loc="Figure 5, left"
)
    
tableResolved.add_image("resolved.pdf")
submission.add_table(tableResolved)

tableBoosted = buildHistTable(
    "boosted.root",
    descr="Observed number of events and predicted number of background events per category for boosted categories",
    loc="Figure 5, right"
)
tableBoosted.add_image("boosted.pdf")
submission.add_table(tableBoosted)


arrMass = np.linspace(2.,20.,18*4+1)
arrMassTau = np.linspace(3.,20.,17*4+1)

### pure ee ###
eeMajorana = buildMassCouplingTable(
    arrMass,
    "ee_majorana_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Majorana HNL pure electron coupling scenario",
    loc="Figure 6, upper left",
    coupling="$|V_{\\text{e}\\text{N}}|^{2}$",
)

submission.add_table(eeMajorana)

eeDirac = buildMassCouplingTable(
    arrMass,
    "ee_dirac_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Dirac HNL pure electron coupling scenario",
    loc="Figure 6, upper right",
    coupling="$|V_{\\text{e}\\text{N}}|^{2}$"
)
submission.add_table(eeDirac)

### pure mumu ###
mumuMajorana = buildMassCouplingTable(
    arrMass,
    "mumu_majorana_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Majorana HNL pure muon coupling scenario",
    loc="Figure 6, middle left",
    coupling="$|V_{\\mu\\text{N}}|^{2}$"
)
submission.add_table(mumuMajorana)

mumuDirac = buildMassCouplingTable(
    arrMass,
    "mumu_dirac_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Dirac HNL pure muon coupling scenario",
    loc="Figure 6, middle right",
    coupling="$|V_{\\mu\\text{N}}|^{2}$"
)
submission.add_table(mumuDirac)


### pure tautau ###
tautauMajorana = buildMassCouplingTable(
    arrMassTau,
    "tautau_majorana_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Majorana HNL pure tau lepton coupling scenario",
    loc="Figure 6, lower left",
    coupling="$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(tautauMajorana)

tautauDirac = buildMassCouplingTable(
    arrMassTau,
    "tautau_dirac_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Dirac HNL pure tau lepton coupling scenario",
    loc="Figure 6, lower right",
    coupling="$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(tautauDirac)




### emu ###
emuMajorana = buildMassCouplingTable(
    arrMass,
    "emu_majorana_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Majorana HNL mixed electron-muon coupling scenario",
    loc="Figure 7, upper left",
    coupling="$|V_{\\text{e}\\text{N}}|^{2}=|V_{\\mu\\text{N}}|^{2}$"
)
submission.add_table(emuMajorana)

emuDirac = buildMassCouplingTable(
    arrMass,
    "emu_dirac_combined_limit_new.npz",
    descr="DTwo-dimensional exclusion limits for Dirac HNL mixed electron-muon coupling scenario",
    loc="Figure 7, upper right",
    coupling="$|V_{\\text{e}\\text{N}}|^{2}=|V_{\\mu\\text{N}}|^{2}$"
)
submission.add_table(emuDirac)


### etau ###
etauMajorana = buildMassCouplingTable(
    arrMass,
    "etau_majorana_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Majorana HNL mixed electron-tau lepton coupling scenario",
    loc="Figure 7, middle left",
    coupling="$|V_{\\text{e}\\text{N}}|^{2}$=$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(etauMajorana)

etauDirac = buildMassCouplingTable(
    arrMass,
    "etau_dirac_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Dirac HNL mixed electron-tau lepton coupling scenario",
    loc="Figure 7, middle right",
    coupling="$|V_{\\text{e}\\text{N}}|^{2}$=$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(etauDirac)


### mutau ###
mutauMajorana = buildMassCouplingTable(
    arrMass,
    "mutau_majorana_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Majorana HNL mixed muon-tau lepton coupling scenario",
    loc="Figure 7, lower left",
    coupling="$|V_{\\mu\\text{N}}|^{2}$=$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(mutauMajorana)

mutauDirac = buildMassCouplingTable(
    arrMass,
    "mutau_dirac_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Dirac HNL mixed muon-tau lepton coupling scenario",
    loc="Figure 7, lower right",
    coupling="$|V_{\\mu\\text{N}}|^{2}$=$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(mutauDirac)


### emutau ###
emutauMajorana = buildMassCouplingTable(
    arrMass,
    "emutau_majorana_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Majorana HNL democratic lepton coupling scenario",
    loc="Figure 8, left",
    coupling="$|V_{\\text{e}\\text{N}}|^{2}$=$|V_{\\mu\\text{N}}|^{2}$=$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(emutauMajorana)

emutauDirac = buildMassCouplingTable(
    arrMass,
    "emutau_dirac_combined_limit_new.npz",
    descr="Two-dimensional exclusion limits for Dirac HNL democratic lepton coupling scenario",
    loc="Figure 8, right",
    coupling="$|V_{\\text{e}\\text{N}}|^{2}$=$|V_{\\mu\\text{N}}|^{2}$=$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(emutauDirac)







### ctau=0.1 ###
ctau01Majorana = buildTriangleCouplingTable(
    "limit_Majorana_ctau0p1.npz",
    descr="Lower limits on Majorana HNL mass for fixed proper decay length of $c\\tau=0.1\\,\\text{mm}$",
    loc="Figure 9, upper left",
    quantity="mass"
)
submission.add_table(ctau01Majorana)

ctau01Dirac = buildTriangleCouplingTable(
    "limit_Dirac_ctau0p1.npz",
    descr="Lower limits on Dirac HNL mass for fixed proper decay length of $c\\tau=0.1\\,\\text{mm}$",
    loc="Figure 9, upper right",
    quantity="mass"
)
submission.add_table(ctau01Dirac)


### ctau=1.0 ###
ctau10Majorana = buildTriangleCouplingTable(
    "limit_Majorana_ctau1p0.npz",
    descr="Lower limits on Majorana HNL mass for fixed proper decay length of $c\\tau=1\\,\\text{mm}$",
    loc="Figure 9, lower left",
    quantity="mass"
)
submission.add_table(ctau10Majorana)

ctau10Dirac = buildTriangleCouplingTable(
    "limit_Dirac_ctau1p0.npz",
    descr="Lower limits on Dirac HNL mass for fixed proper decay length of $c\\tau=1\\,\\text{mm}$",
    loc="Figure 9, lower right",
    quantity="mass"
)
submission.add_table(ctau10Dirac)




### mass=4.5 ###
mass45Majorana = buildTriangleCouplingTable(
    "limit_Majorana_mass4p5.npz",
    descr="Lower limits on Majorana HNL proper decay length for fixed mass of 4.5 GeV",
    loc="Figure 10, upper left",
    quantity="$c\\tau$"
)
submission.add_table(mass45Majorana)

mass45Dirac = buildTriangleCouplingTable(
    "limit_Dirac_mass4p5.npz",
    descr="Lower limits on Dirac HNL proper decay length for fixed mass of 4.5 GeV",
    loc="Figure 10, upper right",
    quantity="$c\\tau$"
)
submission.add_table(mass45Dirac)


### mass=8.0 ###
mass8Majorana = buildTriangleCouplingTable(
    "limit_Majorana_mass8p0.npz",
    descr="Lower limits on Majorana HNL proper decay length for fixed mass of 8 GeV",
    loc="Figure 10, lower left",
    quantity="$c\\tau$"
)
submission.add_table(mass8Majorana)

mass8Dirac = buildTriangleCouplingTable(
    "limit_Dirac_mass8p0.npz",
    descr="Lower limits on Dirac HNL proper decay length for fixed mass of 8 GeV",
    loc="Figure 10, lower right",
    quantity="$c\\tau$"
)
submission.add_table(mass8Dirac)


submission.create_files("output",remove_old=True)






