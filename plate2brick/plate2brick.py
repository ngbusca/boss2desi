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

class DESIDatum:
	lam=dict()
	for band in desi_spec:
		lam[band]=sp.linspace(desi_spec[band][0],desi_spec[band][1],num=desi_spec[band][1]-desi_spec[band][0])

	def __init__(self,row,loglam,flux,ivar):
		flux=interpolate.interp1d(loglam,flux)
		ivar=interpolate.interp1d(10**h[1].data.loglam,h[1].data.ivar)

		self.header=dict()
		self.header["OBJTYPE"] = row.CLASS[0]
		self.header["TARGETCAT"] = row.CLASS[0]
		self.header["TARGETID"] = row.THING_ID[0]
		self.header["TARGET_MASK0"]=0
		self.header["MAG"] = row.FIBERMAG[0]
		self.header["FILTER"] = "UGRIZ"
		self.header["SPECTROID"] = row.THING_ID[0]
		self.header["POSITIONER"] = 0
		self.header["FIBER"] = 0
		self.header["LAMBDAREF"] = 1215.67
		self.header["RA_TARGET"] = row.RA[0]
		self.header["DEC_TARGET"] = row.DEC[0]
		self.header["RA_OBS"] = row.RA[0]
		self.header["DEC_OBS"] = row.DEC[0]
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

class DESIData:

	def __init__(self,spall,plate_dir,plate):

		sp = fits.open(spall)
		fi=glob.glob(plate_dir+"/spPlate-"+str(plate)+"*.fits")
		spPlate=fits.open(fi[0])
		
		w_elg = (sp[1].data.PLATE==plate) & (sp[1].data.CLASS=='GALAXY') & (sp[1].data.ZWARNING_NOQSO==0)
		bt1=40+sp.arange(3)
		for b in bt1:
			w_elg = w_elg | ((sp[1].data.EBOSS_TARGET1 & b) == 1)

		bt2=40+sp.arange(8)
		for b in bt2:
			w_elg = w_elg | ((sp[1].data.EBOSS_TARGET2 & b) == 1)

			
		print "found: ",len(sp[1].data.FIBERID[w_elg])," elgs in plate ",plate

		data=[]

		for row in sp[1].data[w_elg]:
			fid=row["FIBERID"]
			flux=spPlate[0].data[fid-1,:]
			ivar=spPlate[1].data[fid-1,:]
			c0=spPlate[0].header["COEFF0"]
			c1=spPlate[0].header["COEFF1"]
			loglam=c0+c1*sp.arange(len(flux))
			data.append(DESIDatum(row,loglam,flux,ivar))


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

