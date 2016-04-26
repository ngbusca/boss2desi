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

class brickobj:
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
		norm = (la[i]-lam)*iv[i-1] + (lam-la[i-1])*iv[i]
		ivar = norm**2/(iv[i-1]*(la[i]-lam)**2 + iv[i]*(lam-la[i-1])**2)

		w=(iv[i-1]==0) | (iv[i]==0)
		norm[w]=0
		flux[w]=0
		ivar[w]=0

		w=norm>0
		flux[w]/=norm[w]

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
