from hepdata_lib import *
import numpy as np
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


def findCoupling2AtMass(mass,massCoupling2):
    crossings = []
    eps=1e-4
    for i in range(massCoupling2.shape[0]-1):
        if (massCoupling2[i,0]<(mass+eps) and massCoupling2[i+1,0]>(mass-eps)) or (massCoupling2[i,0]>(mass-eps) and massCoupling2[i+1,0]<(mass+eps)):
            coupl = np.interp(mass, massCoupling2[i:i+2,0], massCoupling2[i:i+2,1])
            if len(crossings)==0:
                crossings.append(coupl)
            elif abs(crossings[-1]-coupl)/coupl>0.2:
                crossings.append(coupl)
    return crossings
    
def getCoupling2Values(arrMass,massCoupling2):
    arrCrossings = -1*np.ones((len(arrMass),3),dtype=np.float32)
    for i in range(len(arrMass)):
        crossings = findCoupling2AtMass(arrMass[i],massCoupling2)
        if len(crossings)==0:
            pass
        elif len(crossings)==1:
            arrCrossings[i,0]=crossings[0]    
        elif len(crossings)==2:
            arrCrossings[i,0]=crossings[0]
            arrCrossings[i,2]=crossings[1]
        elif len(crossings)==3:
            arrCrossings[i,0]=crossings[0]
            arrCrossings[i,1]=crossings[1]
            arrCrossings[i,2]=crossings[2]
        elif len(crossings)==4 or len(crossings)==5:
            for j in range(3):
                idx = np.argmin(np.abs(arrCrossings[i-1,j]-crossings))
                if abs(arrCrossings[i-1,j]-crossings[idx])/arrCrossings[i-1,j]<4.:
                    arrCrossings[i,j] = crossings[idx]
                else:
                    print ("skip: ",arrMass[i],arrCrossings[i-1,j],crossings)
            
        else:
            raise Exception("Unknown number of crossings: "+str(len(crossings)),str(crossings))
    return arrCrossings
    
def buildMassCouplingTable(name,arrMass,data,descr="",coupling=""):
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
            var.add_qualifier("UNIT", coupling)
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


    table = Table(name)
    table.description = descr
    table.add_variable(addVar("HNL mass","HNL mass",arrMass,dep=True,units="GeV"))
    
    addPart(table,obsCoupling2,"Observed")
    addPart(table,expCoupling2,"Expected")
    addPart(table,exp68UpCoupling2,"Expected 68% Up")
    addPart(table,exp68DownCoupling2,"Expected 68% Down")
    addPart(table,exp95UpCoupling2,"Expected 95% Up")
    addPart(table,exp95DownCoupling2,"Expected 95% Down")
        
    return table

    

arrMass = np.linspace(2.,20.,18*4+1)
arrMassTau = np.linspace(3.,20.,17*4+1)

### pure ee ###
eeMajorana = buildMassCouplingTable(
    "Table 1",arrMass,
    np.load("ee_majorana_combined_limit_new.npz"),
    descr="Majorana HNL pure electron coupling scenario",
    coupling="$|V_{e\\text{N}}|^{2}$"
)
submission.add_table(eeMajorana)

eeDirac = buildMassCouplingTable(
    "Table 2",arrMass,
    np.load("ee_dirac_combined_limit_new.npz"),
    descr="Dirac HNL pure electron coupling scenario",
    coupling="$|V_{e\\text{N}}|^{2}$"
)
submission.add_table(eeDirac)

### pure mumu ###
mumuMajorana = buildMassCouplingTable(
    "Table 3",arrMass,
    np.load("mumu_majorana_combined_limit_new.npz"),
    descr="Majorana HNL pure muon coupling scenario",
    coupling="$|V_{\\mu\\text{N}}|^{2}$"
)
submission.add_table(mumuMajorana)

mumuDirac = buildMassCouplingTable(
    "Table 4",arrMass,
    np.load("mumu_dirac_combined_limit_new.npz"),
    descr="Dirac HNL pure muon coupling scenario",
    coupling="$|V_{\\mu\\text{N}}|^{2}$"
)
submission.add_table(mumuDirac)


### pure tautau ###
tautauMajorana = buildMassCouplingTable(
    "Table 5",arrMassTau,
    np.load("tautau_majorana_combined_limit_new.npz"),
    descr="Majorana HNL pure tau lepton coupling scenario",
    coupling="$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(tautauMajorana)

tautauDirac = buildMassCouplingTable(
    "Table 6",arrMassTau,
    np.load("tautau_dirac_combined_limit_new.npz"),
    descr="Dirac HNL pure tau lepton coupling scenario",
    coupling="$|V_{\\tau\\text{N}}|^{2}$"
)
submission.add_table(tautauDirac)


submission.create_files("output",remove_old=True)






