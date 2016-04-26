import scipy as sp
from scipy import random
from astropy.io import fits
import fitsio
import glob
import sys
from scipy import interpolate
import os

ndiag=11

boss_bands=["b","r"]

FID_LREF=250
LMIN_BOSS=3600.
LMAX_BOSS=10000.

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
col_type["BOSS_CLASS"]	= "50A"		##	int32
col_type["BOSS_SUBCLASS"]= "50A"	##	int32
col_type["BOSS_Z"]	= "E"		##	float32
col_type["BOSS_PLATE"]	= "J"		##	int32
col_type["BOSS_MJD"]	= "J"		##	int32
col_type["BOSS_ZWARNOQSO"]= "J"		##	int32

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
		self.header["FIBER"] = row["FIBERID"]
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
		self.header["BOSS_PLATE"]=row["PLATE"]
		self.header["BOSS_MJD"]=row["MJD"]
		self.header["BOSS_ZWARNOQSO"]=row["ZWARNING_NOQSO"]

		self.fl=dict()
		self.iv=dict()
		self.wd=dict()
		self.ll=dict()
		self.nexp=dict()
		for band in boss_bands:
			self.nexp[band]=0

	def __call__(self,exp,band,lam):
		ll = self.ll[band][exp]
		fl = self.fl[band][exp]
		iv = self.iv[band][exp]
		wd = self.wd[band][exp]

		la=10**ll
		i = sp.searchsorted(la,lam)
		w=i>=len(ll)
		i[w]=len(ll)-1
		w=i==0
		i[w]=1

		flux = (la[i]-lam)*fl[i-1]*iv[i-1] + (lam-la[i-1])*fl[i]*iv[i]
		ivar = (la[i]-lam)*iv[i-1] + (lam-la[i-1])*iv[i]

		w=(iv[i-1]==0) | (iv[i]==0)
		ivar[w]=0
		flux[w]=0

		w=ivar>0
		flux[w]/=ivar[w]
		ivar/=la[i]-la[i-1]


		wdisp=(la[i]-lam)*wd[i-1] + (lam-la[i-1])*wd[i]
		wdisp/=la[i]-la[i-1]
		re = sp.exp(-(sp.arange(ndiag)-ndiag/2)[:,None]**2/2./wdisp**2)
		re/=sp.sum(re,axis=0)
		return flux,ivar,re

	def add_exp(self,band,ll,fl,iv,wd):
		self.nexp[band]+=1
		if not band in self.fl:
			self.fl[band]=[]
			self.iv[band]=[]
			self.wd[band]=[]
			self.ll[band]=[]

		self.iv[band].append(iv)
		self.fl[band].append(fl)
		self.wd[band].append(wd)
		self.ll[band].append(ll)

	def get_resolution(self,band,exp,lam):
		re=sp.exp(-(sp.arange(ndiag)-ndiag/2)[:,None]**2/2./self.wd[band][exp](lam)**2)
		re/=sp.sum(re,axis=0)
		return re

class DESIData:

	def __init__(self,spall,plate_dir,plates=None,ssample=None):

		spa = fits.open(spall)
#		w_0 = (spa[1].data.CLASS=='GALAXY') & (spa[1].data.ZWARNING_NOQSO==0) & (spa[1].data.THING_ID > 0)
		w_0 = spa[1].data.THING_ID > 0
		if not ssample is None:
			r=sp.random.uniform(size=len(spa[1].data.THING_ID))
			w_0 = w_0 & (r<ssample)
			print "retained: ",len(spa[1].data.THING_ID[w_0])," of ",len(spa[1].data.THING_ID)
		if not plates is None:
			w_plates = sp.zeros(len(spa[1].data.PLATE),dtype=bool)
			for plate in plates:
				w_plates = w_plates | (spa[1].data.PLATE==plate)
			w_0 = w_0 & w_plates

		bt1=[61]

		w_elg = sp.zeros(len(spa[1].data.PLATE),dtype=bool)
		
		for b in bt1:
			w_elg = w_elg | (w_0 & ((spa[1].data.ANCILLARY_TARGET1 & 2**b) > 0))

		bt2=[18,34,39]

		for b in bt2:
			w_elg = w_elg | (w_0 & ((spa[1].data.ANCILLARY_TARGET2 & 2**b) > 0))
			
		print "found: ",len(spa[1].data.FIBERID[w_elg])," elgs "

		rows=spa[1].data[w_elg]
		if plates is None:
			plates = spa[1].data.PLATE[w_elg]


		self.plates = plates

		self.data = dict()
		self.lref = dict()

		for row in rows:
			dat=DESIDatum(row)
			key = str(row["PLATE"])+"-"+str(row["MJD"])+"-"+str(row["FIBERID"])
			self.data[key]=dat


		uplates = sp.unique(plates)

		for plate in uplates:
			fids = [row["FIBERID"] for row in rows if row["PLATE"]==plate]
			mjds = [row["MJD"] for row in rows if row["PLATE"]==plate]
			for band in boss_bands:
				cframes = glob.glob(plate_dir+"/"+str(plate)+"/spCFrame-"+band+"?-*.fits.gz")
				for cframe in cframes:
					print "reading "+cframe
					h=fitsio.FITS(cframe)
					fl=h[0][:,:]
					iv=h[1][:,:]
					ma=h[2][:,:]
					la=h[3][:,:]
					wd=h[4][:,:]
					w=wd==0
					wd[w]=100.
					fibers = [fib for fib in h[5]["FIBERID"][:]]
					if not plate in self.lref:
						self.lref[plate]=dict()

					if not band in self.lref[plate]:
						lam = 10**la[FID_LREF-1]
						w=(LMIN_BOSS<lam) & (lam<LMAX_BOSS)
						self.lref[plate][band]=lam[w]

					print "read"
					for (mjd,fid) in zip(mjds,fids):
						pmf=str(plate)+"-"+str(mjd)+"-"+str(fid)
						if fid in fibers:
							print "adding exposure: ",cframe,band,pmf
							i = fibers.index(fid)
							self.data[pmf].add_exp(band,la[i,:],fl[i,:],iv[i,:]*(ma[i,:]==0),wd[i,:])

					h.close()

		spa.close()

	@staticmethod
	def export(data,dirout="./"):

		for p in sp.unique(data.plates):
			objs = [d for d in data.data.values() if d.header["BOSS_PLATE"]==p]
			print "exporting ",len(objs),"objects in plate ",p
			for band in boss_bands:
				nspec = 0
				for d in objs:
					nspec+=d.nexp[band]

				nbins = len(data.lref[p][band])
				fl=sp.zeros([nspec,nbins])
				iv=sp.zeros([nspec,nbins])
				re=sp.zeros([nspec,ndiag,nbins])

				count = 0
				for d in objs:
					for exp in range(d.nexp[band]):
						fl[count,:],iv[count,:],re[count,:]=d(exp,band,data.lref[p][band])
						count+=1

				hdu0=fits.PrimaryHDU(fl)
				hdu1=fits.ImageHDU(iv)
				hdu2=fits.ImageHDU(data.lref[p][band])
				hdu3=fits.ImageHDU(re)

				hdu0.update_ext_name("FLUX")
				hdu1.update_ext_name("IVAR")
				hdu2.update_ext_name("WAVELENGTH")
				hdu3.update_ext_name("RESOLUTION")

				cols=[]
				index=0
				for key in col_type:
					a=[]
					for d in objs:
						for exp in range(d.nexp[band]):
							if key=="INDEX":
								d.header[key]=index
								index+=1
							if key=="FILTER":d.header[key]=band.upper()
							a.append(d.header[key])
					cols.append(fits.Column(name=key,format=col_type[key],array=a))

				hdu4=fits.BinTableHDU.from_columns(cols)
				hdu4.update_ext_name("FIBERMAP")
				hdulist=fits.HDUList([hdu0,hdu1,hdu2,hdu3,hdu4])
				hdulist.writeto(dirout+"/brick-"+band+"-"+str(p)+"p000.fits",clobber=True)
