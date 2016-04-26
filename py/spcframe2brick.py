import scipy as sp
from scipy import random
from astropy.io import fits
import fitsio
import glob
import sys
from scipy import interpolate
import os
from brickobj import *

class brick:

	def __init__(self,spall,plate_dir,plates=None,ssample=None):

		spa = fits.open(spall)
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
			dat=brickobj(row)
			key = str(row["PLATE"])+"-"+str(row["MJD"])+"-"+str(row["FIBERID"])
			self.data[key]=dat


		uplates = sp.unique(plates)

		for plate in uplates:
			fids = [row["FIBERID"] for row in rows if row["PLATE"]==plate]
			mjds = [row["MJD"] for row in rows if row["PLATE"]==plate]
			for band in boss_bands:
				cframes = []
				for di in plate_dir:
					spcfs=glob.glob(di+"/"+str(plate)+"/spCFrame-"+band+"?-*.fits.gz")
					cframes = sp.concatenate([cframes,glob.glob(di+"/"+str(plate)+"/spCFrame-"+band+"?-*.fits.gz")])
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
