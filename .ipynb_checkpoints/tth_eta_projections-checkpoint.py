from scipy import stats
import numpy as np
def tth_eta_projections(list_of_images, cen_pos, tth_nbins, eta_nbins):
    def xy2rphi(x, y, x0, y0):
        r = np.sqrt((x-x0)**2+(y-y0)**2)
        phi = np.arctan((y-y0)/(x-x0))*180/np.pi
        return r, phi

    def rphi2xy(r, phi):
        x=r*np.cos(phi)
        y=r*np.sin(phi)
        return x, y
    
    img = list_of_images[0]
    y = np.linspace(0, np.shape(img)[0]-1, np.shape(img)[0])
    x = np.linspace(0, np.shape(img)[1]-1, np.shape(img)[1])
    yy,xx=np.meshgrid(y,x, sparse=False, indexing='ij')
    rr, phiphi= xy2rphi(xx,yy,cen_pos[1],cen_pos[0])

    rr_1d=np.reshape(rr,np.shape(rr)[0]*np.shape(rr)[1])
    phi_1d=np.reshape(phiphi,np.shape(rr)[0]*np.shape(rr)[1])
    tth_omega = []
    eta_omega = []
    for img in list_of_images:
        img_1d=np.reshape(img,np.shape(rr)[0]*np.shape(rr)[1])
        r_binned, r_edges, _ = stats.binned_statistic(rr_1d,img_1d, statistic='mean', bins=tth_nbins)
        phi_binned, phi_edges, _ = stats.binned_statistic(phi_1d,img_1d, statistic='mean', bins=eta_nbins)
        tth_omega.append(r_binned)
        eta_omega.append(phi_binned)
        
    r_plot = r_edges[:-1]+np.diff(r_edges)/2
    phi_plot = phi_edges[:-1]+np.diff(phi_edges)/2
    
    return tth_omega, eta_omega