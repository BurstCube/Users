import numpy as np
import matplotlib.pylab as plot
from astropy.io import ascii,fits
from scipy import interpolate
import grb_catalogs
from BurstCube.LocSim.Detector import *
from BurstCube.LocSim.Spacecraft import *
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.optimize import curve_fit
from astropy.table import Table
import healpy as hp
from gammaray_proposal_tools import *

### run code

def run(dir='/Users/jracusin/BurstCube/gitrep/Users/jracusin/',nsims=10000,minflux=0.5):

	burstcube, BCpointings, aeff_bc = setup_BC(dir=dir)
	fermi, GBMpointings, aeff_gbm=setup_GBM(dir=dir)

	## Aeff at 100 keV
	# bcaeff=loginterpol(aeff_bc['keV'],aeff_bc['aeff'],150.)
	# gbmaeff=loginterpol(aeff_gbm['energy'],aeff_gbm['aeff'],150.)
	# print(bcaeff,gbmaeff)

	#Aeff on same energy points
	eng=np.logspace(np.log10(50),np.log10(300),100)
	bcaeff=loginterpol(aeff_bc['keV'],aeff_bc['aeff'],eng)
	gbmaeff=loginterpol(aeff_gbm['energy'],aeff_gbm['aeff'],eng)

#	print(bcaeff/gbmaeff)

	trig,gbm=load_GBM_catalogs(dir=dir)
	s=np.where(gbm['T90']<=2.0)[0]
	sgbm=gbm[s]
	print(len(sgbm))
	# realgbmflux=sgbm['FLUX_BATSE_1024']
	# wreal=np.where(realgbmflux>0)[0]

	interval=1.0 #s
	bgrate=300. #cts/s in 50-300 keV
	gbmexposures, bcexposures, secondhighestgbm, secondhighestbc, randgbmexposures, randbcexposures=throw_grbs(fermi,burstcube,nsims=nsims)
#	simgbmcr,simbccr,simgbmpfsample,simbcpfsample,realpf,pinterval=grb_spectra(sgbm,gbmaeff,bcaeff,eng,nsims,interval=interval)
	gbmflux2counts,bcflux2counts,realpf=grb_spectra(sgbm,gbmaeff,bcaeff,eng,nsims,interval=interval)
	pf=logNlogS(bcaeff,gbmaeff,minflux=minflux,nsims=nsims,interval=interval)
	r=np.array(np.round(np.random.rand(nsims)*(len(realpf)-1)).astype('int'))


	simgbmcr=pf*gbmflux2counts[r]
	simbccr=pf*bcflux2counts[r]
	simgbmpfsample=pf
	simbcpfsample=pf
	pinterval=1.

#	simgbmcr,simbccr,simgbmpfsample,simbcpfsample=logNlogS(bcaeff,gbmaeff,minflux=minflux,nsims=nsims,interval=interval)

	realgbmflux=realpf
	wreal=np.where(realgbmflux>0)[0]

	pf=simgbmpfsample

	#Solve for the number of detected counts which will equal our source photons
	sourcegbm = simgbmcr*secondhighestgbm*pinterval
	sourcebc = simbccr*secondhighestbc*pinterval
	#randomize background rate around typical background of 300 cts/s (50-300 keV, GBM)
	bckgrd=np.random.poisson(bgrate,nsims)
	scaledgbmbckgrd = bckgrd*pinterval
	scaledbcbckgrd = bckgrd*np.median(bcaeff/gbmaeff)*pinterval
	#creating an array of zeros that I can manipulate to create an array of detected GRBs
	detectgbm = np.zeros(len(sourcegbm))
	detectbc = np.zeros(len(sourcebc))

	#calculate the significance of the second highest exposure detector. If the significance is greater than 4.5 sigma than the burst is detectable.
	for u in range(len(sourcegbm)):
		if sourcegbm[u]>0:
			sig = sourcegbm[u] / (np.sqrt(sourcegbm[u] + scaledgbmbckgrd[u]))
			if sig > 4.5:
			    detectgbm[u] = 1.0
			else:
			    detectgbm[u] = 0.0

	for j in range(len(sourcebc)):
		if sourcebc[j]>0:
			sig = sourcebc[j] / (np.sqrt(sourcebc[j] + scaledbcbckgrd[j]))
			if sig > 4.5:
			    detectbc[j] = 1.0
			else:
			    detectbc[j] = 0.0
		else: sig=0

	#Creating plot of peak flux versus counts for real and simulated GBM
	w=np.where(pf>0)[0]
	wg = np.where(simgbmcr*detectgbm>0.)[0]
	wbc = np.where(simbccr*detectbc>0.)[0]

	fig=plot.figure(figsize=(10,8))
	plot.subplot(2,2,1)
#	plot.hist(gbmcr[w],label='real GBM',bins=np.logspace(1,6,40),color='orange')
	plot.hist(simgbmcr[wg],label='GBM',bins=np.logspace(1,6,40),alpha=0.7,color='blue')
	plot.hist(simbccr[wbc],label='BurstCube',bins=np.logspace(1,6,40),alpha=0.7,color='green')
	plot.xlabel('Count Rate (50-300 keV; cts/s)')
	plot.xscale('log')
	plot.yscale('log')
	plot.xlim([10,5e4])
	plot.ylabel('N Simulated  sGRBs')

	plot.legend()
	plot.subplot(2,2,2)
	plot.hist(simgbmpfsample,label='Simulated total',bins=np.logspace(-1,4,40),alpha=1.0,color='C3')
	plot.hist(realgbmflux[wreal],label='real GBM',bins=np.logspace(-1,4,40),color='orange', alpha=0.7)
	# this is the simulated GBM
	plot.hist(simgbmpfsample[wg],label='GBM',bins=np.logspace(-1,4,40),alpha=0.5,color='blue')
	plot.hist(simbcpfsample[wbc],label='BC',bins=np.logspace(-1,4,40),alpha=0.5,color='green')
	plot.xlabel('Peak Flux (50-300 keV; ph/cm2/s)')
	#plot.hist(flux[w],label='BC',bins=np.logspace(-1,2,40),alpha=0.7,color='red')
	plot.xscale('log')
	plot.yscale('log')
	plot.xlim([.1,300])
	plot.legend()
	plot.ylabel('N Simulated  sGRBs')

#	plot.show()

	#solve for the detection fraction of BurstCube and Simulated GBM
	detgbm = np.where(detectgbm == 1)[0]
	ratiogbm = float(len(detgbm)) / float(len(detectgbm))
	print(ratiogbm)

	detbc = np.where(detectbc == 1)[0]
	ratiobc = float(len(detbc)) / float(len(detectbc))
	print(ratiobc)

	print('fraction of GBM sGRBs BC will detect = %0.2f'%(ratiobc/ratiogbm))
	#number of bursts BurstCube will see a year
	bcbursts = ratiobc/ratiogbm *40.
	print('bc rate = %.2f'%bcbursts+' sGRBs/yr')

	### Duty Cycle to detect 20 sGRBs/yr
	gbmduty=0.85
	duty=20./(bcbursts/gbmduty)
	print("duty cycle to detect 20 sGRBs/yr = %.2f" %duty)
	duty=10./(bcbursts/gbmduty)
	print("duty cycle to detect 10 sGRBs/yr = %.2f" %duty)

	### Min sensitivity to detect 10 per year
	nbursts10=bcbursts-10.
	nbursts20=bcbursts-20.
	so=np.argsort(simbcpfsample[wbc])
	gso=np.argsort(simgbmpfsample[wg])
	c=np.cumsum(np.ones(len(wbc)))/len(wbc)*bcbursts
	plot.subplot(2,2,3)
	plot.plot(simbcpfsample[wbc[so]],c)
	plot.xlabel(r'BurstCube 50-300 keV Peak Flux (ph cm$^{-2}$ s$^{-1}$)')
	plot.ylabel('Cumulative Number')
	plot.xscale('log')
	fluxlim10=loginterpol(c,simbcpfsample[wbc[so]],nbursts10)
	fluxlim20=loginterpol(c,simbcpfsample[wbc[so]],nbursts20)
	plot.plot([fluxlim10,fluxlim10],[nbursts10,nbursts10],marker='*',label='Limit for 10 sGRBs')
	plot.plot([fluxlim20,fluxlim20],[nbursts20,nbursts20],marker='*',label='Limit for 20 sGRBs')
	plot.xlim([1,100])

	print("flux limit to detect 10 sGRBs/yr = %.2f"%fluxlim10+' ph/cm2/s')
	print("flux limit to detect 20 sGRBs/yr = %.2f"%fluxlim20+' ph/cm2/s')
	print('expected minimum flux = '+"%.2f"%min(simbcpfsample[wbc[so]])+' ph/cm2/s')
	print('expected maximum flux = '+"%.2f"%max(simbcpfsample[wbc[so]])+' ph/cm2/s')
	print('expected 5% maximum flux = '+"%.2f"%simbcpfsample[wbc[so[int(0.05*len(so))]]]+' ph/cm2/s')
	print('expected 10% maximum flux = '+"%.2f"%simbcpfsample[wbc[so[int(0.1*len(so))]]]+' ph/cm2/s')
	print('expected 90% maximum flux = '+"%.2f"%simbcpfsample[wbc[so[int(0.9*len(so))]]]+' ph/cm2/s')
	print('expected 95% maximum flux = '+"%.2f"%simbcpfsample[wbc[so[int(0.95*len(so))]]]+' ph/cm2/s')

	# print('GBM')
	# print('expected minimum flux = '+"%.2f"%min(simgbmpfsample[wg[gso]])+' ph/cm2/s')
	# print('expected maximum flux = '+"%.2f"%max(simgbmpfsample[wg[gso]])+' ph/cm2/s')
	# print('expected 5% maximum flux = '+"%.2f"%simgbmpfsample[wg[gso[int(0.05*len(gso))]]]+' ph/cm2/s')
	# print('expected 10% maximum flux = '+"%.2f"%simgbmpfsample[wg[gso[int(0.1*len(gso))]]]+' ph/cm2/s')
	# print('expected 90% maximum flux = '+"%.2f"%simgbmpfsample[wg[gso[int(0.9*len(gso))]]]+' ph/cm2/s')
	# print('expected 95% maximum flux = '+"%.2f"%simgbmpfsample[wg[gso[int(0.95*len(gso))]]]+' ph/cm2/s')

	## FoV - adjusted exposure alt until total reached 20
	BCFoVrad = 90-0. # deg radius
	BCFoV=(1-np.cos(np.radians(BCFoVrad)))/2.*4.*np.pi
#	print("FoV for "+"%.1f" % BCFoV+' ster')

	## max distance of GW170817
	mpc2cm=3.086e24
	fgw=3.7 # ph/cm2/s
	fmax=min(simgbmpfsample[wg])
	dgw=42.9*mpc2cm
	dmax=np.sqrt(fgw*dgw**2/fmax)
	f=80.*mpc2cm/dmax
	print("%.2f" % (dmax/mpc2cm*f)+' Mpc - distance GBM for GW170817')

	fmax=min(simbcpfsample[wbc])
	dmax=np.sqrt(fgw*dgw**2/fmax)
	print("%.2f" % (dmax/mpc2cm*f)+' Mpc - distance BC for GW170817')

	### mission lifetime to detect 10 sGRBs
	print("Mission Duration to detect 10 sGRBs = " + "%.1f" % (10./bcbursts*12.)+' months')
	plot.legend()
	plot.show()

#	return realgbmflux,simgbmpfsample

def setup_BC(dir=''):

	burstcube, BCpointings, Aeff, index=load_mission('BurstCube')

	## read in BurstCube Aeff for various BC configurations
	bcaeffs=ascii.read(dir+'BC_eff_area_curves.ecsv',format='ecsv')
	w=np.where((bcaeffs['diameter']==90) & (bcaeffs['height']==19) )
	aeff_bc=bcaeffs[w]
#	eng_bc=bcaeffs['keV'][w]

	return burstcube, BCpointings, aeff_bc#, eng_bc

def setup_GBM(dir=''):

	fermi, GBMpointings, Aeff, index=load_mission('GBM')
	## read in the GBM Aeff
	aeff_gbm = np.genfromtxt(dir+'gbm_effective_area.dat',skip_header=2,names=('energy', 'aeff'))

	return fermi, GBMpointings, aeff_gbm

def load_GBM_catalogs(dir=''):
	#read in GBM Trigger Catalog
	trigfit=fits.open(dir+'gbmtrigcat.fits')
	trig=trigfit[1].data

	#read in GBM Burst Catalog
	gbmfit=fits.open(dir+'gbmgrbcat.fits')
	gbm=gbmfit[1].data

	return trig,gbm

# now that GBM and BurstCube's pointings are set up we will throw GRBs at it and determine the exposure for each GRB. 
#generate GRBs and throw them at GBM

def throw_grbs(fermi,burstcube,nsims=10000):
    
    ra,dec=random_sky(nsims)
    ra=np.array(ra)-180
    dec=np.array(dec)
   
    #GBM and BurstCube exposures for each random GRB.
    randgbmexposures = np.array([[detector.exposure(ra[i],dec[i], alt=-23.,index=0.78) for i in range(nsims)] for detector in fermi.detectors])
    randbcexposures = np.array([[detector.exposure(ra[i],dec[i], alt=-23.,index=0.6) for i in range(nsims)] for detector in burstcube.detectors])
    
    #Order randgbmexposures into descending order
    for column in randgbmexposures.T:
        newrandgbm = -np.sort(-randgbmexposures.T) 
    gbmexposures = np.transpose(newrandgbm)
    
    for col in randbcexposures.T:
        newrandbc = -np.sort(-randbcexposures.T) 
    bcexposures = np.transpose(newrandbc)

    #Select the second highest exposure value. 
    #We will use this to ensure the second highest exposure detector has a sig >4.5
    secondhighestgbm = gbmexposures[1,:]
    secondhighestbc = bcexposures[1,:]
        
    return gbmexposures, bcexposures, secondhighestgbm, secondhighestbc, randgbmexposures, randbcexposures

def logNlogS(aeff_bc,aeff_gbm,minflux=0.5,nsims=10000,interval=1.0):

    #1 sec 50-300 keV peak flux ph/cm2/s
	time = interval#1.0#0.064 # s
	f=np.logspace(np.log10(minflux),2.2,50)
	p=f**-0.9#1.5 # comes from fitting GBM sGRB logN-log peak flux
	pnorm=p/np.sum(p)
	r=np.random.choice(f,p=pnorm,size=nsims)

	# bg_gbm=bgrate*time
	# bg_bc=bgrate*np.max(aeff_bc)/np.max(aeff_gbm)*time  # scaling from GBM average background rate
	src_bc=r*np.max(aeff_bc)*time
	src_gbm=r*np.max(aeff_gbm)*time

	simgbmpfsample = np.array(r)
	simgbmcr = np.array(src_gbm/time)
	simbcpfsample = np.array(r)
	simbccr = np.array(src_bc/time)

	return r#simgbmcr,simbccr,simgbmpfsample,simbcpfsample

def grb_spectra(gbmbursts,gbmaeff,bcaeff,eng,nsims,interval=1.0):

	#Integrating the best fit spectrum for each GRB in the energy range of 50-300 keV to get max. observed photon flux. 
	#Doing the same but also folding in the effective area in order to get count rate.
	#This will give us the photon flux in units of ph/cm^2/s. 
	mo=gbmbursts['PFLX_BEST_FITTING_MODEL']
	pf=np.zeros(len(mo))
	gbmcr=np.zeros(len(mo))
	bccr=np.zeros(len(mo))
	pflux_interval=np.zeros(len(mo))
	realpf=np.zeros(len(mo))

	for i in range(len(mo)):
#	    for j in range(len(gbmbursts)):
#	        Aratio=(aeff_bc/aeff_gbm)
	        # this should give us an array of the maximum observed photon flux for GBM
		if mo[i]=='PFLX_PLAW':
			gbmcr[i]=np.trapz(gbmbursts['PFLX_PLAW_AMPL'][i]*grb_catalogs.pl(eng,gbmbursts['PFLX_PLAW_INDEX'][i])*gbmaeff,eng)
			pf[i]=np.trapz(gbmbursts['PFLX_PLAW_AMPL'][i]*grb_catalogs.pl(eng,gbmbursts['PFLX_PLAW_INDEX'][i]),eng)
			bccr[i]=np.trapz(gbmbursts['PFLX_PLAW_AMPL'][i]*grb_catalogs.pl(eng,gbmbursts['PFLX_PLAW_INDEX'][i])*bcaeff,eng)
			realpf[i]=gbmbursts['PFLX_PLAW_PHTFLUXB'][i]

		if mo[i]=='PFLX_COMP':
			gbmcr[i]=np.trapz(gbmbursts['PFLX_COMP_AMPL'][i]*grb_catalogs.comp(eng,gbmbursts['PFLX_COMP_INDEX'][i],gbmbursts['PFLX_COMP_EPEAK'][i])*gbmaeff,eng)
			pf[i]=np.trapz(gbmbursts['PFLX_COMP_AMPL'][i]*grb_catalogs.comp(eng,gbmbursts['PFLX_COMP_INDEX'][i],gbmbursts['PFLX_COMP_EPEAK'][i]),eng)
			bccr[i]=np.trapz(gbmbursts['PFLX_COMP_AMPL'][i]*grb_catalogs.comp(eng,gbmbursts['PFLX_COMP_INDEX'][i],gbmbursts['PFLX_COMP_EPEAK'][i])*bcaeff,eng)
			realpf[i]=gbmbursts['PFLX_COMP_PHTFLUXB'][i]

		if mo[i]=='PFLX_BAND':
			gbmcr[i]=np.trapz(gbmbursts['PFLX_BAND_AMPL'][i]*grb_catalogs.band(eng,gbmbursts['PFLX_BAND_ALPHA'][i],gbmbursts['PFLX_BAND_EPEAK'][i],gbmbursts['PFLX_BAND_BETA'][i])*gbmaeff,eng)
			pf[i]=np.trapz(gbmbursts['PFLX_BAND_AMPL'][i]*grb_catalogs.band(eng,gbmbursts['PFLX_BAND_ALPHA'][i],gbmbursts['PFLX_BAND_EPEAK'][i],gbmbursts['PFLX_BAND_BETA'][i]),eng)
			bccr[i]=np.trapz(gbmbursts['PFLX_BAND_AMPL'][i]*grb_catalogs.band(eng,gbmbursts['PFLX_BAND_ALPHA'][i],gbmbursts['PFLX_BAND_EPEAK'][i],gbmbursts['PFLX_BAND_BETA'][i])*bcaeff,eng)
			realpf[i]=gbmbursts['PFLX_BAND_PHTFLUXB'][i]

		if mo[i]=='PFLX_SBPL':
			gbmcr[i]=np.trapz(gbmbursts['PFLX_SBPL_AMPL'][i]*grb_catalogs.sbpl(eng,gbmbursts['PFLX_SBPL_INDX1'][i],gbmbursts['PFLX_SBPL_BRKEN'][i],gbmbursts['PFLX_SBPL_INDX2'][i])*gbmaeff,eng)
			pf[i]=np.trapz(gbmbursts['PFLX_SBPL_AMPL'][i]*grb_catalogs.sbpl(eng,gbmbursts['PFLX_SBPL_INDX1'][i],gbmbursts['PFLX_SBPL_BRKEN'][i],gbmbursts['PFLX_SBPL_INDX2'][i]),eng)
			bccr[i]=np.trapz(gbmbursts['PFLX_SBPL_AMPL'][i]*grb_catalogs.sbpl(eng,gbmbursts['PFLX_SBPL_INDX1'][i],gbmbursts['PFLX_SBPL_BRKEN'][i],gbmbursts['PFLX_SBPL_INDX2'][i])*bcaeff,eng)
			realpf[i]=gbmbursts['PFLX_SBPL_PHTFLUXB'][i]


		pflux_interval[i]=gbmbursts['PFLX_SPECTRUM_STOP'][i]-gbmbursts['PFLX_SPECTRUM_START'][i]

	flux=gbmbursts['FLUX_BATSE_1024']
	gbmflux2counts=gbmcr/pf
	bcflux2counts=bccr/pf
#	fluxwrong=flux/pf#*pflux_interval
	r=np.array(np.round(np.random.rand(nsims)*(len(mo)-1)).astype('int'))
	simgbmcr=gbmcr[r]*interval#*fluxwrong[r]#*pflux_interval[r]
	simbccr=bccr[r]*interval#*fluxwrong[r]#*pflux_interval[r]
	simpf=pf[r]*interval#*fluxwrong[r]#*pflux_interval[r]
	pinterval=pflux_interval[r]
	realflux=flux[r]

	return gbmflux2counts,bcflux2counts,realpf#simgbmcr,simbccr,simpf,simpf,realpf,pinterval