import numpy as nm
import matplotlib.pyplot as pl
import pandas as pd

# load in the spreadsheet
df_inst_inputs = pd.read_excel('instrument_definition_inputs.xlsx')

# always skip first data point, because f=0 causes problem with CMB integral

# load loading data
tmp = nm.loadtxt('am_runs/JJA_50.dat')
f_GHz  = tmp[1:,0]
atm_tx = tmp[1:,2]
T_atm  = tmp[1:,3]

# create a function to calculate the loading
def calculate_net_and_nep(f_GHz,system_efficiency,T_det,atmosphere_transmission, dish_diameter, P0_override=False):
    # constants [SI]
    h = 6.626068e-34
    k = 1.3806503e-23
    c = 299792458.0

    # the way this will work is a little different
    # to handle the integral properly, assume that each frequency in the f_GHz array is its own detector
    # with bandwidth df
    # calculate the net and nep and everything ALL the way through to the end, for each of these "detectors"
    # then at the very end, take the noise weighted coadded average of all of them
    # that should yield the same answer for a tophant band of ones and zeros
    # and should yield the correct answer for a general bandbass function

    # convert from GHz
    f = f_GHz*1.0e9

    df = nm.mean(nm.diff(f))

    # calculate power vs frequency for RJ source
    P0 = k*T_det*df

    shot = 2*k*T_det*h*f*df
    wave = 2*(P0)**2.0 / df
    print('Shot noise = '+str(nm.sqrt(nm.sum(shot))*1.0e18)+' aWrtHz')
    print('Wave noise = '+str(nm.sqrt(nm.sum(wave))*1.0e18)+' aWrtHz')
    # this nep will be "low" because bad optical efficiency has reduced the as-seen loading
    nep = nm.sqrt(shot+wave)# watts root second

    # do the NET_CMB integral (Lueker's Thesis, 2.23)
    # in this scaling, bad optical efficiency makes the NET high again
    # since NEP_shot drops as sqrt(OE), but this equation embiggens as OE
    # that means NET does scale up as sqrt(OE) as expected
    net_integrand = ((h*f)/(2.725))**2.0 * (1.0/k) * ( nm.exp( (h*f)/(k*2.725) ) / ((nm.exp( (h*f)/(k*2.725) ) - 1)**2.0) )
    net = (nep / (nm.sqrt(2.0) * system_efficiency * net_integrand * df) )*1.0e6 # uK root second 

    # nerj integrand, Lueker 2.21 without an extra 2
    nerj_integrand = nm.sqrt(2)*k*system_efficiency*df
    nerj = (nep / nerj_integrand)*1.0e6 # uK root second

    # scale the NEP 
    nep = nep*1.0e18 #attowatts root hz

    ## multiply both by the fact that detector noise = 0.3 * blip
    #nep = nep*nm.sqrt(1.3)
    #net = net*nm.sqrt(1.3)

    # calculate neft using the dish size
    nefd = (nep*1e-18*2.0e26*1e3) / ((nm.pi*(dish_diameter/2.0)**2)*df) # mJrtHz
    # and system efficieicny since nep is the at-detector number
    nefd = nefd/system_efficiency
    # and convert from rtHz to rts
    nefd = nefd/nm.sqrt(2.0)

    # scale NET and NEFD up by atmospheric transmission
    net = net/atmosphere_transmission
    nefd = nefd/atmosphere_transmission

    # net, nefd, and nep are all arrays now
    # https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
    # sigma^2 = 1 / sum(sigma^(-2))
    net = nm.sqrt(1.0 / nm.sum(net**(-2.0)))
    nerj = nm.sqrt(1.0 / nm.sum(nerj**(-2.0)))
    nefd = nm.sqrt(1.0 / nm.sum(nefd**(-2.0)))
    # nep is sum of squares
    nep = nm.sqrt(nm.sum(nep**2.0))
    # power just adds
    P0 = nm.sum(P0)

    return P0,net,nefd,nep,nerj

# make numpy arrays for the outputs
Ndets_array = nm.zeros(len(df_inst_inputs))
beams = nm.zeros(len(df_inst_inputs))
loading_at_det = nm.zeros(len(df_inst_inputs)) 
nep_at_det = nm.zeros(len(df_inst_inputs)) 
nefd_array = nm.zeros(len(df_inst_inputs))
nerj_array = nm.zeros(len(df_inst_inputs))

# go through and run the model for each row (except the last row that just adds up the number of detector wafers)
for i in range(len(df_inst_inputs)):
    fmin = df_inst_inputs['fmin'][i]
    fmax = df_inst_inputs['fmax'][i]
    fcent_nom = (fmin+fmax)/2.0

    passband = nm.zeros_like(f_GHz) + 1e-10 # literally zero will make noise infinite, which may mess up the sum later?
    passband[(f_GHz>fmin) & (f_GHz<fmax)] = 1.0

    pl.ion()
    pl.figure(100)
    pl.plot(f_GHz,passband)
    pl.plot(f_GHz,atm_tx,'k')
    pl.xlabel('[GHz]')
    pl.ylabel('Passband')

    # detector efficiency
    detector_optical_efficiency = 0.95*0.5 # filter bank

    # aperture efficiency for 1flambda horns
    aperture_efficiency = 0.35

    primary_mirror_optical_efficiency = 0.95

    # add excess instrument loading
    T_excess = 10.0 # K

    # scale this down by detector optical efficiency to make the temperature at the detector
    T_det = (T_atm+T_excess)*detector_optical_efficiency*passband*aperture_efficiency

    telescope_diameter = 50.0

    fwhm_arcsec = (50.0/telescope_diameter)*(280.0/fcent_nom)*5.0
    print('Nominal beam assumed to be '+str(fwhm_arcsec)+' arcseconds')
    beams[i] = fwhm_arcsec

    fnum = 2 # f number at focal plane array

    Ndet = df_inst_inputs['Nchan'][i]*df_inst_inputs['Nhorns'][i]
    Ndets_array[i] = Ndet


    # calculate system efficiency
    system_efficiency = detector_optical_efficiency*passband*aperture_efficiency*primary_mirror_optical_efficiency


    # single-detector noise
    # this calculates the loading assuming the T_det is the at-detector loading
    # and then calculates NET and NEFD by scaling that noise by the system efficiency and the atmospheric transmission
    P0,net,nefd,nep_at_detector,nerj = calculate_net_and_nep(f_GHz,
                                                        system_efficiency,
                                                        T_det,
                                                        atm_tx,
                                                        telescope_diameter)
    print('Center Frequency = '+str(fcent_nom)+' GHz')
    print('Loading = '+str(P0*1.0e12)+' pW')
    print('NEP at Detector = '+str(nep_at_detector)+' aWrts')
    print('NEFD = '+str(nefd)+' mJrts')

    loading_at_det[i] = P0*1.0e12
    nep_at_det[i] = nep_at_detector 
    nefd_array[i] = nefd
    nerj_array[i] = nerj


# add these rows to the existing data frame, to be saved
df_inst_inputs['Detector Count'] = Ndets_array
df_inst_inputs['Beam in Arcseconds'] = beams
df_inst_inputs['Loading at Detector in pW'] = loading_at_det
df_inst_inputs['NEP at Detectors in aWrtHz'] = nep_at_det
df_inst_inputs['NEFD of a Single Detector'] = nefd_array
df_inst_inputs['NET_RJ of a Single Detector'] = nerj_array
df_inst_inputs['NET_RJ of Entire Frequency Band'] = nerj_array/nm.sqrt(Ndets_array)

# save this out to an excel file
df_inst_inputs.to_csv('results_of_instrument_simulation.csv')
