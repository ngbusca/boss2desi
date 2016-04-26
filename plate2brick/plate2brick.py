import scipy as sp
from scipy import random
from astropy.io import fits
import fitsio
import glob
import sys
from scipy import interpolate
import os

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

class DESIDatum:
	def __init__(self,row):
		self.header=dict()
		self.header["OBJTYPE"] = row["CLASS"]
		self.header["TARGETCAT"] = row["OBJTYPE"]
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
		self.header["BOSS_CLASS"]=row["CLASS"]
		self.header["BOSS_SUBCLASS"]=row["SUBCLASS"]
		self.header["BOSS_Z"]=row["Z"]

		self.fl=dict()
		self.iv=dict()
		self.wd=dict()

	def add_band(self,band,exp,lam,ivar,flux,wdisp):
		band_exp=(band,exp)

		self.iv[band_exp]=interpolate.interp1d(lam,ivar)
		self.fl[band_exp]=interpolate.interp1d(lam,flux)
		self.wd[band_exp]=interpolate.interp1d(lam,wdisp)

	def get_resolution(self,band,exp,lam):
		band_exp = (band,exp)
		re=sp.exp(-(sp.arange(ndiag)-ndiag/2)[:,None]**2/2./self.wd[band_exp](lam)**2)
		re/=sp.sum(re,axis=0)
		return re

class DESIData:

	def __init__(self,spall,plate_dir,plate=None):

		spa = fits.open(spall)
		w_0 = (spa[1].data.CLASS=='GALAXY') & (spa[1].data.ZWARNING_NOQSO==0) & (spa[1].data.THING_ID>0)
		if plate is not None:
			w_0 = w_0 & (spa[1].data.PLATE==plate)
		bt1=[61]

		w_elg = sp.zeros(len(spa[1].data.PLATE),dtype=bool)
		
		for b in bt1:
			w_elg = w_elg | (w_0 & ((spa[1].data.ANCILLARY_TARGET1 & 2**b) > 0))

		bt2=[18,34,39]

		for b in bt2:
			w_elg = w_elg | (w_0 & ((spa[1].data.ANCILLARY_TARGET2 & 2**b) > 0))
			
		print "found: ",len(spa[1].data.FIBERID[w_elg])," elgs "

		plates=spa[1].data.PLATE[w_elg]
		mjds=spa[1].data.MJD[w_elg]
		rows=spa[1].data[w_elg]

		name_row = dict()
		
		count=0
		for (plate,mjd) in zip(plates,mjds):
			key = str(plate)+"-"+str(mjd)
			if not name_row.has_key(key):
				name_row[key] = []

			name_row[key].append(rows[count])
			count+=1

		names = [str(plate)+"-"+str(mjd) for (plate,mjd) in zip(plates,mjds)]
		names = sp.unique(names)

		self.data=[]

		for name in names:
			spPlate = fitsio.FITS(plate_dir+"/spPlate-"+name+".fits")
			print "reading plate "+name," ELGs: ",len(name_row[name])

			for row in name_row[name]:
				fid = row["FIBERID"]
				flux=spPlate[0][fid-1,:].flatten()
				ivar=spPlate[1][fid-1,:].flatten()
				amask=spPlate[2][fid-1,:].flatten()
				wdisp=spPlate[4][fid-1,:].flatten()
				w=(amask!=0) | (wdisp==0)
				ivar[w]=1e-10
				wdisp[w]=100.
				head=spPlate[0].read_header()
				c0=head["COEFF0"]
				c1=head["COEFF1"]
				loglam=c0+c1*sp.arange(len(flux))
				dat=DESIDatum(row)
				dat.add_band("coadd",0,10**loglam,ivar,flux,wdisp)
				self.data.append(dat)

			spPlate.close()

		spa.close()

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

