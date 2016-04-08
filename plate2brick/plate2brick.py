import scipy as sp
from scipy import random
from astropy.io import fits
import glob
import sys
from scipy import interpolate
import os

desi_spec=dict()
desi_spec["b"]=(3600,5930)
desi_spec["r"]=(5660,7720)
desi_spec["z"]=(7470,9800)


ndiag=11

col_type=dict()
col_type["OBJTYPE"]  	= "10A" 	##	char[10]	 	 
col_type["TARGETCAT"]	= "20A"		##	char[20]	 	 
col_type["TARGETID"] 	= "K"		##	int64	 	 
col_type["TARGET_MASK0"]= "K" 		##	int64
col_type["MAG"]		= "5E"		##	float32[5]	 	 
col_type["FILTER"]	= "50A"		##	char[50]	 	 
col_type["SPECTROID"]	= "K"		##	int64	 	 
col_type["POSITIONER"]	= "K"		##	int64	 	 
col_type["FIBER"]	= "J"		##	int32	 	 
col_type["LAMBDAREF"]	= "E"		##	float32	 	 
col_type["RA_TARGET"]	= "D"		##	float64	 	 
col_type["DEC_TARGET"]	= "D"		##	float64	 	 
col_type["RA_OBS"]	= "D"		##	float64	 	 
col_type["DEC_OBS"]	= "D"		##	float64	 	 
col_type["X_TARGET"]	= "D"		##	float64	 	 
col_type["Y_TARGET"]	= "D"		##	float64	 	 
col_type["X_FVCOBS"]	= "D"		##	float64	 	 
col_type["Y_FVCOBS"]	= "D"		##	float64	 	 
col_type["Y_FVCERR"]	= "E"		##	float32	 	 
col_type["X_FVCERR"]	= "E"		##	float32	 	 
col_type["NIGHT"]	= "J"		##	int32	 	 
col_type["EXPID"]	= "J"		##	int32	 	 
col_type["INDEX"]	= "J"		##	int32
col_type["EBOSS_CLASS"]	= "50A"		##	int32
col_type["EBOSS_Z"]	= "E"		##	int32

def edge(l,l0,l1,typ,sigma=4):
	eff=l*0+1
	if typ=="blue":
		w=l<l1
		eff[w]-=sp.exp(-(l[w]-l0)**2/2*sigma**2/(l1-l0)**2)
	elif typ=="red":
		w=l>l0
		eff[w]=sp.exp(-(l[w]-l0)**2/2*sigma**2/(l1-l0)**2)

	return eff

class DESIDatum:
	lmin=3599
	lmax=10000
	l=sp.linspace(lmin,lmax,num=lmax-lmin)
	lam=dict()
	for band in desi_spec:
		w=(l>desi_spec[band][0]) & (l<desi_spec[band][1])
		lam[band]=l[w]

	eff=dict()
	eff["b"]=edge(lam["b"],desi_spec["r"][0],desi_spec["b"][1],"red")
	eff["r"]=edge(lam["r"],desi_spec["r"][0],desi_spec["b"][1],"blue")*edge(lam["r"],desi_spec["z"][0],desi_spec["r"][1],"red")
	eff["z"]=edge(lam["z"],desi_spec["z"][0],desi_spec["r"][1],"blue")


	def __init__(self,row,loglam,flux,ivar,wdisp):
		flux=interpolate.interp1d(10**loglam,flux)
		ivar=interpolate.interp1d(10**loglam,ivar)
		wdisp=interpolate.interp1d(10**loglam,wdisp)

		self.header=dict()
		self.header["OBJTYPE"] = row["CLASS"]
		self.header["TARGETCAT"] = row["CLASS"]
		self.header["TARGETID"] = row["THING_ID"]
		self.header["TARGET_MASK0"]=0
		self.header["MAG"] = row["FIBERMAG"]
		self.header["FILTER"] = "UGRIZ"
		self.header["SPECTROID"] = row["THING_ID"]
		self.header["POSITIONER"] = 0
		self.header["FIBER"] = 0
		self.header["LAMBDAREF"] = 1215.67
		self.header["RA_TARGET"] = row["RA"]
		self.header["DEC_TARGET"] = row["DEC"]
		self.header["RA_OBS"] = row["RA"]
		self.header["DEC_OBS"] = row["DEC"]
		self.header["X_TARGET"] = 0
		self.header["Y_TARGET"] = 0
		self.header["X_FVCOBS"] = 0
		self.header["Y_FVCOBS"] = 0
		self.header["Y_FVCERR"] = 0
		self.header["X_FVCERR"] = 0
		self.header["NIGHT"] = 0
		self.header["EXPID"] = 0
		self.header["INDEX"] = 0
		self.header["EBOSS_CLASS"]=row["CLASS"]
		self.header["EBOSS_Z"]=row["Z"]

		self.fl=dict()
		self.iv=dict()
		self.re=dict()

		for band in desi_spec:
			noise=sp.random.normal(size=len(self.lam[band]))
			self.iv[band]=ivar(self.lam[band])*self.eff[band]
			self.fl[band]=flux(self.lam[band])+noise*sp.sqrt((1-self.eff[band])/self.iv[band])
			self.re[band]=sp.exp(-(sp.arange(ndiag)-ndiag/2)[:,None]**2/2./wdisp(self.lam[band])**2)
			self.re[band]/=sp.sum(self.re[band],axis=0)

class DESIData:

	def __init__(self,spall,plate_dir,plate):

		spa = fits.open(spall)
		fi=glob.glob(plate_dir+"/spPlate-"+str(plate)+"*.fits")
		spPlate=fits.open(fi[0])
		
		w_0 = (spa[1].data.PLATE==plate) & (spa[1].data.CLASS=='GALAXY') & (spa[1].data.ZWARNING_NOQSO==0)  
		bt1=40+sp.arange(3)

		w_elg = sp.zeros(len(spa[1].data.PLATE),dtype=bool)
		
		for b in bt1:
			w_elg = w_elg | (w_0 & ((spa[1].data.EBOSS_TARGET1 & b) > 0))

		bt2=40+sp.arange(8)


		for b in bt2:
			w_elg = w_elg | (w_0 & ((spa[1].data.EBOSS_TARGET2 & b) > 0))
			
		print "found: ",len(spa[1].data.FIBERID[w_elg])," elgs in plate ",plate

		self.data=[]

		for row in spa[1].data[w_elg]:
			fid=row["FIBERID"]
			flux=spPlate[0].data[fid-1,:]
			ivar=spPlate[1].data[fid-1,:]
			amask=spPlate[2].data[fid-1,:]
			#w=amask!=0
			#ivar[w]=0
			wdisp=spPlate[4].data[fid-1,:]
			w=wdisp==0
			wdisp[w]=100.
			c0=spPlate[0].header["COEFF0"]
			c1=spPlate[0].header["COEFF1"]
			loglam=c0+c1*sp.arange(len(flux))
			dat=DESIDatum(row,loglam,flux,ivar,wdisp)
			self.data.append(dat)

		spa.close()
		spPlate.close()

	@staticmethod
	def export(data,sufix=None):
		for band in desi_spec:
			fl=sp.zeros([len(data),desi_spec[band][1]-desi_spec[band][0]])
			iv=sp.zeros([len(data),desi_spec[band][1]-desi_spec[band][0]])
			re=sp.zeros([len(data),ndiag,desi_spec[band][1]-desi_spec[band][0]])

			for (i,d) in enumerate(data):
				fl[i,:]=d.fl[band]
				iv[i,:]=d.iv[band]
				re[i,:]=d.re[band]

			hdu0=fits.PrimaryHDU(fl)
			hdu1=fits.ImageHDU(iv)
			hdu2=fits.ImageHDU(DESIDatum.lam[band])
			hdu3=fits.ImageHDU(re)

			hdu0.update_ext_name("FLUX")
			hdu1.update_ext_name("IVAR")
			hdu2.update_ext_name("WAVELENGTH")
			hdu3.update_ext_name("RESOLUTION")

			cols=[]
			for key in col_type:
				a=[]
				for (i,d) in enumerate(data):
					if key=="INDEX":d.header[key]=i
					if key=="FILTER":d.header[key]=band.upper()
					a.append(d.header[key])
				cols.append(fits.Column(name=key,format=col_type[key],array=a))

			hdu4=fits.BinTableHDU.from_columns(cols)
			hdu4.update_ext_name("FIBERMAP")
			hdulist=fits.HDUList([hdu0,hdu1,hdu2,hdu3,hdu4])
			hdulist.writeto("brick-"+band+"-"+sufix+".fits",clobber=True)

