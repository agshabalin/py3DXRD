import numpy as np
from skimage.transform import warp_polar
from skimage.util import img_as_float


def calc_projections(imgs, mask, cenpos):
    maxrad = max( list(np.asarray(imgs[0].shape)-np.asarray(cenpos)) + cenpos )
    if type(mask) !=np.ndarray:
        mask = np.ones(imgs[0].shape, dtype=np.float32)
    imgs_sum = mask*imgs[0]
    img = warp_polar(img_as_float(mask*imgs[0]), center=cenpos, radius=maxrad)
    eta_tth = img
    omg_eta = np.zeros([len(imgs), img.shape[0]], dtype=np.float32)
    omg_tth = np.zeros([len(imgs), img.shape[1]], dtype=np.float32)
    omg_eta[0,:] = eta_tth.sum(1)
    omg_tth[0,:] = eta_tth.sum(0)

    for i in range(1,len(imgs)):
        imgs_sum += mask*imgs[i]
        img = warp_polar(img_as_float(mask*imgs[i]), center=cenpos, radius=maxrad)
        eta_tth += img
        omg_eta[i,:] = img.sum(1)
        omg_tth[i,:] = img.sum(0)
    
    del cenpos, maxrad, img

    return imgs_sum, eta_tth, omg_eta, omg_tth