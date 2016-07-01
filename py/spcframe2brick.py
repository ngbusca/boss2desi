import scipy as sp
from scipy import random
from astropy.io import fits
import fitsio
import glob
import sys
from scipy import interpolate
import os
from brickobj import *

delta_mjd = 2 ## difference of mjd values in a single plugin

class brick:

	def __init__(self,rows,plate_dir,mask=False):

		plates = [row["PLATE"] for row in rows]
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
					cframes = sp.concatenate([cframes,glob.glob(di+"/"+str(plate)+"/spCFrame-"+band+"?-*.fits*")])
				for cframe in cframes:
					print "reading "+cframe
					h=fitsio.FITS(cframe)
					head = h[0].read_header()
					spcframe_mjd = head["MJD"]
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
						if abs(mjd-spcframe_mjd)>delta_mjd:continue
						if fid in fibers:
							print "adding exposure: ",cframe,band,pmf
							i = fibers.index(fid)
							self.data[pmf].add_exp(band,la[i,:],fl[i,:],iv[i,:]*(ma[i,:]==0),wd[i,:])

					h.close()


	@staticmethod
	def export(data,dirout="./",mask=False,pext="p000"):

		inv_plate_map = None
		if mask:
			thid_map,plate_map = brick.mask_data(data)
			f=open(dirout+"/mask_map.dat","a")
			for th in thid_map:
				f.write(str(th)+" "+str(thid_map[th])+"\n")

			f.close()
			inv_plate_map = dict()
			for p in plate_map:
				inv_plate_map[plate_map[p]]=p




		for p in sp.unique(data.plates):
			objs = [d for d in data.data.values() if d.header["BOSS_PLATE"]==p]
			print "exporting ",len(objs),"objects in plate ",p
			plam = p
			if inv_plate_map !=None:
				plam = inv_plate_map[p]
			for band in boss_bands:
				nspec = 0
				for d in objs:
					nspec+=d.nexp[band]
				nbins = len(data.lref[plam][band])
				fl=sp.zeros([nspec,nbins])
				iv=sp.zeros([nspec,nbins])
				re=sp.zeros([nspec,ndiag,nbins])

				count = 0
				for d in objs:
					for exp in range(d.nexp[band]):
						fl[count,:],iv[count,:],re[count,:]=d(exp,band,data.lref[plam][band])
						count+=1

				hdu0=fits.PrimaryHDU(fl)
				hdu1=fits.ImageHDU(iv)
				hdu2=fits.ImageHDU(data.lref[plam][band])
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
				hdulist.writeto(dirout+"/brick-"+band+"-"+str(p)+pext+".fits",clobber=True)

	@staticmethod
	def mask_data(data):
		thids = [d.header["TARGETID"] for d in data.data.values()]
		thids = sp.unique(thids)
		thid_map = dict()

		for th in thids:
			old_th = str(th)
			new_th = ''
			for (i,oth) in enumerate(old_th):
				aux = random.randint(9)
				if aux==0 and i==0:aux=1
				new_th = new_th+str(aux)
			thid_map[th]=int(new_th)
		
		plates = [d.header["BOSS_PLATE"] for d in data.data.values()]
		plates = sp.unique(plates)
		print plates
		plate_map = dict()
		for (i,p) in enumerate(plates):
			plate_map[p]=i
		
		data.plates = plate_map.values()
		for d in data.data.values():
			d.header["OBJTYPE"] = "HIDDEN"
			d.header["TARGETCAT"]="HIDDEN"
			d.header["TARGETID"] = thid_map[d.header["TARGETID"]]
			d.header["MAG"]=23
			d.header["SPECTROID"] = d.header["TARGETID"]
			d.header["FIBERID"]=0
			d.header["RA_TARGET"]=3.14
			d.header["DEC_TARGET"]=3.14
			d.header["RA_OBS"]=3.14
			d.header["DEC_OBS"]=3.14
			d.header["BOSS_CLASS"]="HIDDEN"
			d.header["BOSS_SUBCLASS"]="HIDDEN"
			d.header["BOSS_Z"]=42
			d.header["BOSS_PLATE"]=plate_map[d.header["BOSS_PLATE"]]
			d.header["BOSS_MJD"] = 43689
			d.header["BOSS_ZWARNINGNOQSO"]=4

		return thid_map,plate_map
