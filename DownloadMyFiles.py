
import DownloadObservations as D
from astroquery.alma import Alma
from matplotlib import pyplot as plt



path = 'C:\\Users\\vmart\\Documents\\KandidatArbete\\Data\\'
projectCodes, data = D.initdatatable()
Alma.login("vmart")
D.downloadFiles(projectCodes,path , True, 'cont')
fitsfiles = D.getAllFits(path)
print(fitsfiles)

for f in fitsfiles:
   if 'pbcor' and 'cont' in f and not 'mask' in f:
        D.astro_plot(f,zoom = -1,sigma=[3],color='purple',save=False)
        plt.show()
