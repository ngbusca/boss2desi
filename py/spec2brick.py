import scipy as sp
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
col_type["TARGETID"] 	= "int64"	##	int64	 	 
col_type["TARGET_MASK0"]= "int64" 	##	int64
col_type["MAG"]		= "5E"		##	float32[5]	 	 
col_type["FILTER"]	= "50A"		##	char[50]	 	 
col_type["SPECTROID"]	= "int64"	##	int64	 	 
col_type["POSITIONER"]	= "int64"	##	int64	 	 
col_type["FIBER"]	= "int32"	##	int32	 	 
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
col_type["NIGHT"]	= "int32"	##	int32	 	 
col_type["EXPID"]	= "int32"	##	int32	 	 
col_type["INDEX"]	= "int32"	##	int32

fi=glob.glob(sys.argv[1]+"/spec*.fits")

class DESIData:
	lam=dict()
	for band in desi_spec:
		lam[band]=sp.linspace(desi_spec[band][0],desi_spec[band][1],num=desi_spec[band][1]-desi_spec[band][0])

	def __init__(self,spec_file):
		h=fits.open(spec_file)

		flux=interpolate.interp1d(10**h[1].data.loglam,h[1].data.flux)
		ivar=interpolate.interp1d(10**h[1].data.loglam,h[1].data.ivar)

		self.header=dict()
		self.header["OBJTYPE"] = h[2].data.CLASS[0]
		self.header["TARGETCAT"] = h[2].data.CLASS[0]
		self.header["TARGETID"] = h[2].data.THING_ID[0]
		self.header["TARGET_MASK0"]=0
		self.header["MAG"] = h[2].data.FIBERMAG[0]
		self.header["FILTER"] = "UGRIZ"
		self.header["SPECTROID"] = h[2].data.THING_ID[0]
		self.header["POSITIONER"] = 0
		self.header["FIBER"] = 0
		self.header["LAMBDAREF"] = 1215.67
		self.header["RA_TARGET"] = h[2].data.RA[0]
		self.header["DEC_TARGET"] = h[2].data.DEC[0]
		self.header["RA_OBS"] = h[2].data.RA[0]
		self.header["DEC_OBS"] = h[2].data.DEC[0]
		self.header["X_TARGET"] = 0
		self.header["Y_TARGET"] = 0
		self.header["X_FVCOBS"] = 0
		self.header["Y_FVCOBS"] = 0
		self.header["Y_FVCERR"] = 0
		self.header["X_FVCERR"] = 0
		self.header["NIGHT"] = 0
		self.header["EXPID"] = 0
		self.header["INDEX"] = 0

		self.fl=dict()
		self.iv=dict()
		self.re=dict()

		for band in desi_spec:
			self.fl[band]=flux(self.lam[band])
			self.iv[band]=ivar(self.lam[band])
			self.re[band]=sp.exp(-(sp.arange(ndiag)-ndiag/2)**2)[:,None]
			self.re[band]/=sp.sum(self.re[band])

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
			hdu2=fits.ImageHDU(DESIData.lam[band])
			hdu3=fits.ImageHDU(re)
			hdu0.update_ext_name("FLUX")
			hdu1.update_ext_name("IVAR")
			hdu2.update_ext_name("WAVELENGTH")
			hdu3.update_ext_name("RESOLUTION")

			cols=[]
			for key in col_type:
				a=[]
				for d in data:
					a.append(d.header[key])
				cols.append(fits.Column(name=key,format=col_type[key],array=a))

			hdu4=fits.BinTableHDU.from_columns(cols)
			hdu4.update_ext_name("FIBERMAP")
			hdulist=fits.HDUList([hdu0,hdu1,hdu2,hdu3,hdu4])
			hdulist.writeto("brick-"+band+"-"+sufix+".fits",clobber=True)


data=[]
for (i,f) in enumerate(fi):
	print f,i,len(fi)
	data.append(DESIData(f))

DESIData.export(data,"test")

