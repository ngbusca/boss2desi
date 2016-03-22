import scipy as sp
from astropy.io import fits
import glob
import sys
from scipy import interpolate
import os

lmin=3600.
dl=1.
nbins=10000-3600
lmax=lmin+nbins*dl
lam=sp.linspace(lmin,lmax,num=nbins)

fi=glob.glob(sys.argv[1]+"/spec*.fits")
da=sp.zeros([len(fi),nbins])
iv=sp.zeros([len(fi),nbins])
ndiag=11
re=sp.zeros([len(fi),ndiag,nbins])
count=0
objtype=[]
targcat=[]
targetid=[]
mag=[]
fil=[]
specid=[]
positioner=[]
fiber=[]
lambdaref=[]
ra=[]
dec=[]
ra_obs=[]
dec_obs=[]
xt=[]
yt=[]
xfv=[]
yfv=[]
dxfv=[]
dyfv=[]
night=[]
expid=[]
index=[]
for i in fi:
	print i
	h=fits.open(i)
	objtype.append(h[2].data.CLASS[0])
	targcat.append(h[2].data.CLASS[0])
	targetid.append(h[2].data.THING_ID[0])
	mag.append(h[2].data.FIBERMAG[0])
	fil.append("UGRIZ")
	specid.append(h[2].data.THING_ID[0])
	positioner.append(0)
	fiber.append(0)
	lambdaref.append(1215.67)
	ra.append(h[2].data.RA[0])
	dec.append(h[2].data.DEC[0])
	ra_obs.append(h[2].data.RA[0])
	dec_obs.append(h[2].data.DEC[0])
	xt.append(0)
	yt.append(0)
	xfv.append(0)
	yfv.append(0)
	dxfv.append(0)
	dyfv.append(0)
	night.append(0)
	expid.append(0)
	index.append(0)


	flux=interpolate.interp1d(10**h[1].data.loglam,h[1].data.flux)
	ivar=interpolate.interp1d(10**h[1].data.loglam,h[1].data.ivar)
	da[count,:]=flux(lam)
	iv[count,:]=ivar(lam)
	res=sp.exp(-(sp.arange(ndiag)-ndiag/2)**2)[:,None]
	res/=sp.sum(res)
	re[count,:,:]=res

cols=[]
cols.append(fits.Column(name="OBJTYPE",format='20A',array=objtype))
cols.append(fits.Column(name="TARGETCAT",format='20A',array=targcat))
cols.append(fits.Column(name="TARGETID",format='int64',array=targetid))
cols.append(fits.Column(name="TARGET_MASK0",format='int64',array=targetid))
cols.append(fits.Column(name="MAG",format='5E',array=mag))
cols.append(fits.Column(name="FILTER",format='50A',array=fil))
cols.append(fits.Column(name="SPECID",format='int64',array=specid))
cols.append(fits.Column(name="POSITIONER",format='int64',array=positioner))
cols.append(fits.Column(name="FIBER",format='int32',array=fiber))
cols.append(fits.Column(name="LAMBDAREF",format='E',array=lambdaref))
cols.append(fits.Column(name="RA_TARGET",format='D',array=ra))
cols.append(fits.Column(name="DEC_TARGET",format='D',array=dec))
cols.append(fits.Column(name="RA_OBS",format='D',array=ra_obs))
cols.append(fits.Column(name="DEC_OBS",format='D',array=dec_obs))
cols.append(fits.Column(name="X_TARGET",format='D',array=xt))
cols.append(fits.Column(name="Y_TARGET",format='D',array=yt))
cols.append(fits.Column(name="X_FVCOBS",format='D',array=xfv))
cols.append(fits.Column(name="Y_FVCOBS",format='D',array=yfv))
cols.append(fits.Column(name="X_FVCERR",format='E',array=dxfv))
cols.append(fits.Column(name="Y_FVCERR",format='E',array=dyfv))
cols.append(fits.Column(name="NIGHT",format='int32',array=night))
cols.append(fits.Column(name="EXPID",format='int32',array=expid))
cols.append(fits.Column(name="INDEX",format='int32',array=index))


tbhdu=fits.BinTableHDU.from_columns(cols)
wb=lam<(lmin+lmax)/2

hdu0=fits.PrimaryHDU(da[:,wb])
hdu1=fits.ImageHDU(iv[:,wb])
hdu2=fits.ImageHDU(lam[wb])
hdu3=fits.ImageHDU(re[:,:,wb])
hdu0.update_ext_name("FLUX")
hdu1.update_ext_name("IVAR")
hdu2.update_ext_name("WAVELENGTH")
hdu3.update_ext_name("RESOLUTION")

hdulist=fits.HDUList([hdu0,hdu1,hdu2,hdu3,tbhdu])
hdulist.writeto("brick-b-"+sys.argv[2]+".fits",clobber=True)

wr=lam>=(lmin+lmax)/2

hdu0=fits.PrimaryHDU(da[:,wr])
hdu1=fits.ImageHDU(iv[:,wr])
hdu2=fits.ImageHDU(lam[wr])
hdu3=fits.ImageHDU(re[:,:,wr])
hdu0.update_ext_name("FLUX")
hdu1.update_ext_name("IVAR")
hdu2.update_ext_name("WAVELENGTH")
hdu3.update_ext_name("RESOLUTION")

hdulist=fits.HDUList([hdu0,hdu1,hdu2,hdu3,tbhdu])
hdulist.writeto("brick-r-"+sys.argv[2]+".fits",clobber=True)
