import numpy as np

from pathlib import Path
from typing import Optional
from collections import namedtuple


ARRAY_CONFIGURATION = namedtuple("ArrayConfiguration", ["stations", "baseline_range"])

AT_SMALL = ARRAY_CONFIGURATION("A0-B2-C1-D0", [11, 34])
AT_MEDIUM = ARRAY_CONFIGURATION("K0-G2-J3-D0", [40, 104])
AT_LARGE = ARRAY_CONFIGURATION("A0-G1-J3-J2", [58, 132])
AT_ASTROMETRIC = ARRAY_CONFIGURATION("A0-G1-K0-J2", [49, 129])
UTS = ARRAY_CONFIGURATION("U1-U2-U3-U4", [47, 130])


DELAY_LINES = {"small": AT_SMALL, "medium": AT_MEDIUM, "large": AT_LARGE,
               "astrometric": AT_ASTROMETRIC, "UTs": UTS}


DATA_DIR = Path(__file__).parents[3] / "data"
DELAY_RESTRICTIONS_DIR = DATA_DIR / "delay_line_restrictions"
CALIBRATOR_CATALOGUES_DIR = DATA_DIR / "calibrator_catalogues"



# NOTE: These only work for the 'paranal' location and VLTI
def get_delay_line_restrictions(array_configuration: str, observatory_name: str):
    """Gets the restrictions for the VLTI/Paranal in azimuth and altitude corresponding to the
    array_configuration

    Parameters
    ----------
    array_configuration: str
        The array configuration. Can be "small", "medium", "large", "astrometric" and
        "UTs2"
    observatory_name: str

    Returns
    -------
    altitude: np.ndarray
    azimuth: np.ndarray
    """
    # NOTE: Maybe make this into a warning
    if observatory_name != "paranal":
        raise IOError("The delay line restrictions are only implemented for paranal!")
    delay_restriction_file_path = DELAY_RESTRICTIONS_DIR / "combined" \
            / f"{DELAY_LINES[array_configuration].stations}.npy"
    return np.load(delay_restriction_file_path, allow_pickle=True)[::-1]


# function insert_zenith,tl,lst,pars
  # tl=[tl[0],tl]
  # tl[0].RA=lst*15d0
  # tl[0].dec=pars.obs_lat
  # tl[0].name='zenith'
  # return,tl
# end

# function fix_Lband_cat,c
  # c.F10*=(3.78/10.)^2 ;; [Jy], assumes rayleigh-jeans regime from 3.78-10 micron
  # return,c
# end

def get_catalogue(update_year: Optional[int] = 2019):
    """Reads the calibrator data from the catalogues"""
    midi_catalogue_path = CALIBRATOR_CATALOGUES_DIR / f"midi_cat_{update_year}.dat"
    cohen_catalogue_path = CALIBRATOR_CATALOGUES_DIR / f"cohen_cat_{update_year}.dat"
    visir_catalogue_path = CALIBRATOR_CATALOGUES_DIR / f"visir_cat_{update_year}.dat"
    msdfcc_catalogue_path = CALIBRATOR_CATALOGUES_DIR / f"msdfcc_cat_{update_year}.dat"
    df = pd.read_csv(midi_catalogue_path)
    print(df)


# function get_catalogue,catalogue
  # ;; the final catalogue is compiled, by first taking the first catalogue,
  # ;; then check which sources in the 2nd catalogue weren't in the first and
  # ;; append these, then for the third.
  # ext='_2019'
  # for i=0,n_elements(catalogue)-1 do begin
     # if (catalogue[i] eq 'midi') then restore,'catalogues/midi_cat'+ext+'.dat'
     # if (catalogue[i] eq 'cohen') then restore,'catalogues/cohen_cat'+ext+'.dat'
     # if (catalogue[i] eq 'visir') then restore,'catalogues/visir_cat'+ext+'.dat'
     # if (catalogue[i] eq 'Lband') then begin
        # restore,'catalogues/msdfcc_cat.dat'
        # nlam=n_elements(calcat[0].spec.lam)
        # calcat=fix_Lband_cat(calcat)
        # if (n_elements(cat[0].spec.lam) gt nlam) then begin
           # spec=calcat.spec
           # calcat=extend_structure(calcat,'spec',empty(cat[0].spec),exclude='SPEC')
           # calcat.spec.lam[0:nlam-1,*]=spec.lam
           # calcat.spec.Fnu[0:nlam-1,*]=spec.Fnu
           # calcat.spec.err[0:nlam-1,*]=spec.err
        # endif
        # if (tag_exist(cat,'VBV_PARS') eq 1) then begin
           # template=empty(cat[0].vbv_pars) & template.match=0 & template.readme=''
           # origin=calcat.origin
           # calcat=extend_structure(calcat,'vbv_pars',template,exclude='ORIGIN')
           # calcat=extend_structure(calcat,'origin',origin[0])
        # endif
     # endif
     # match_dist=5. ;; arcsec
     
     # if (i eq 0) then cat=calcat ;; the first one
     # if (i ge 1) then begin      ;; more than one catalogue, we need to merge them!
        # distances=fltarr(n_elements(calcat))
        # for j=0,n_elements(calcat)-1 do begin
           # dist=3600.*sqrt(((calcat[j].coords.RA-cat.coords.RA)*cos(calcat[j].coords.dec*!pi/180.))^2+$
                           # (calcat[j].coords.dec-cat.coords.dec)^2) ;; [arcsec]
           # distances[j]=min(dist,ix)
        # endfor
        # indx=where(distances gt match_dist)
        # if (i eq 2) then print,catalogue[i]
        # if (i eq 2) then help,/str,cat,calcat
        # if (indx[0] ne -1) then cat=[cat,calcat[indx]]
     # endif
  # endfor
  
  # order=sort(cat.coords.RA) & cat=cat[order]
  # return,cat
# end

# function flux_limit,calcat,pars,ATs=ATs,baseline=baseline,FINITO=FINITO
  # maxdec=90. ;; degrees

  # if keyword_set(ATs) then limit=pars.min_F10_ATs else limit=pars.min_F10_UTs
  # ix=where((calcat.F10 ge limit) and (calcat.coords.dec le maxdec))

  # calcat=calcat[ix]

  # if keyword_set(FINITO) then begin
    # visibilities=fltarr(n_elements(calcat))
    # for i=0,n_elements(calcat)-1 do visibilities[i]= $
      # calibrator_visibility([1.6],calcat[i].diam.theta,baseline)
    # ix=where(visibilities ge pars.finito_Hvis_limit)
    # calcat=calcat[ix]
    # extrastring=', with H-band visibility above '+strtrim(string(pars.finito_Hvis_limit,format='(f10.2)'),2)
    # calcat=calcat[ix]
  # endif else extrastring=''

# ;  ncal=n_elements(calcat)
# ;  print,'there are '+strtrim(string(ncal),2)+' calibrators brighter than '+strtrim(string(limit,format='(i10)'),2)+' Jy at declinations more southern than delta='+strtrim(string(maxdec,format='(i10)'),2)+extrastring
  # return,calcat
# end



# Pro calibrator_find,target,$         ;; name of science target as listed in targetsfile
                    # zenith=zenith,$  ;; if set, the (hypothetical) science target is put at the zenith
                    # start=start,$    ;; [hours], time at which to start the observation set (default: now).
                                     # ;; an observation set consists of a science target and a calibrator
                                     # ;; by default, the science target will be done first and the calibrator second
                    # lst=lst,$        ;; [hours] or ['hours:minutes'] LST at which ot start the evaluation of calibrators (overwrites "start")
                                     # ;; LST can be given in dicimal notation or as a string of 'hours:minutes', so for example
                                     # ;; LST=5.25 and LST='05:15' both work and are equivalent
                    # delay_restrict=delay_restrict,$ ;; optional, one of ['small','medium','large','astrometric','UTs'], draws "pointing" limits
                                                    # ;; for given cofiguration (shadowing, delay line restrictions)
                    # feeling_lucky=feeling_lucky,$   ;; if set, use more optimistic shadowing limits (with the risk that a UT door is in the way!)
                                                    # ;; --> basically don't do this unless you are in despair!!!
                    # noPlot_restrict=noPlot_restrict,$ ;; if set, do not mark the time ranges when science target or calibrator are not observable
                                                      # ;; in the airmass vs. time plot (default: indicate these interval(s) with blue diamonds)
                    # cal=cal,$                       ;; use a specific calibrator (name), useful for detailed planning
                    # duration=duration,$             ;; [min] duration of an observation (default: 30 min)
                    # do_cal_before=do_cal_before,$   ;; if set, we start with the calibrator and then do the science target
                    # offset_cal=offset_cal,$         ;; [min] if set, put a buffer/pause of this length between the science and calibration observation
                    # print_plan=print_plan,$         ;; if set, print lines for integration in observing plan description (cal='...' must be set)
                    # LN=LN,$                         ;; if set, indicate the calibrator is for both L and N band in the observing plan
                    # nmajor=nmajor,$                 ;; use specific start of time for the evaluation of calibrators (in old "timmi2" convention)
                    # ;; input file for science targets, and catalogues to use:
                    # targetsfile=targetsfile,$       ;; file in which science target coordinates are listed (default: sources.txt)
                    # select_sources=select_sources,$ ;; array of source names. If set, only these sources from the targetsfile will be plotted
                                                    # ;; on the sky plot (default: all sources from the targetsfile)
                    # catalogue=catalogue,$           ;; which catalogue to use, one of ['midi','cohen','Lband'], default: midi+cohen
                    # minF10=minF10,$        ;; minimum 10 micron flux of calibrator (defaults taken from "parameters.txt")
                    # minFL=minFL,$          ;; minimumb L-band flux of calibrator (if set, it overwrites minF10)
                    # max_diam=max_diam,$    ;; [mas], if set, limit calibrators to those with diameter <= this value
                    # max_d_az=max_d_az,$    ;; [deg], maximum difference in azimuth (default taken from parameters.txt)
                    # max_d_am=max_d_am,$    ;; [dimensionless], maximum difference in airmass (default taken from parameters.txt)
                    # max_d_ang=max_d_ang,$  ;; [deg], maximum difference in angular distance (default taken from parameters.txt)
                    # zoom=zoom,$            ;; optional, show smaller part of sky (default, entire sky (4*pi sr) , also parts below horizon)
                    # plot_airmasses=plot_airmasses,$ ;; either set /plot_airmasses to plot default grid, or do plot_airmasses=[am1,am2,...,amN]
                                                    # ;; to explicitly specify which airmasses to plot
                    # nolabels=nolabels,$             ;; if set, don't print the names of targets on the sky plot
                    # galactic_plane=galactic_plane,$ ;; if set, draw the galactic plane. If simply /galactic_plane is set, a band of +/-5 deg
                                        # ;; around b=0 is drawn. One cal also set galactic_plane=x, in which case a band of +/-x degrees is drawn
                    # halt=halt,$         ;; end with a "stop" so that variables are still available in meory
                    # baseline=baseline,$  ;; [m] baseline length for which to evaluate calibrator visibilities (if delay_restrict is set, then
                                         # ;; the calibrator visiblities will be calculated for the longest baseline in the particular array,
                                         # ;; ignoring projection effects. Setting the baseline keyword overrides this)
                    # ;; obsolete parameters
                    # ATs=ATs,$     ;; adopt the default flux limits for ATs (default: no flux limit is applied)
                    # FINITO=FINITO,$ ;; if set, calibrator visibility in H-band must be >= minimum value set in parameters.txt
                    # now=now


  # ;; catalogue: choose between 'midi', 'cohen', 'visir', 'msdfcc_cat' (=Pierre's L-band catalogue)
  # ;; combinations are possible, e.g. catalogue=['midi','cohen']

  # if keyword_set(now) then LST=0
  # if not keyword_set(duration) then duration=30. ;; [min]
  # pars=read_parameters('parameters.txt')
  # if keyword_set(lst) then LST=lst else LST=get_LST(pars,start=start) 
  # lst=check_lst_format(lst)
  # if keyword_set(nolabels) then pars.plot_names='no'
  
  # ;; read the catalogue(s)
  # if not keyword_set(catalogue) then catalogue=['midi','cohen']
  # calcat=get_catalogue(catalogue)
  # if ((not keyword_set(baseline)) and (keyword_set(delay_restrict))) then baseline=max_baseline(delay_restrict) ;; [m]
  # if not keyword_set(baseline) then baseline=100. ;; [m]
  # if not keyword_set(feeling_lucky) then aspro_limits=1 else aspro_limits=0
  
  
  # if keyword_set(delay_restrict) then $
    # delrestr=read_delay_line_restrictions(delay_restrict,aspro_limits=aspro_limits)
  
  # if (keyword_set(minF10) or keyword_set(minFL)) then begin
     # if keyword_set(minFL) then minF10=minFL*(3.78/10.)^2 ;; [Jy]
     # pars.min_F10_UTs=minF10 & pars.min_F10_ATs=minF10  ;; [Jy]
  # endif
  # if keyword_set(max_d_az) then pars.max_d_az=max_d_az
  # if keyword_set(max_d_am) then pars.max_d_am=max_d_am
  # if keyword_set(max_d_ang) then pars.max_d_ang=max_d_ang
  
  # if not keyword_set(targetsfile) then file='sources.txt' else file=targetsfile
  # target_list=read_targetlist(file)
  # if keyword_set(select_sources) then target_list=included_stars(target_list)
  # if keyword_set(cal) then begin
     # new=calcat[where(calcat.name eq cal[0])]
     # if (n_elements(cal) gt 1) then for j=1,n_elements(cal)-1 do new=[new,calcat[where(calcat.name eq cal[j])]]
     # calcat=new
     # max_diam=0. & minF10=0. & max_d_az=0. & max_d_am=0.
# ;     get_lun,lun
# ;     openu,lun,'/home/pboley/.xephem/fifos/xephem_in_fifo'
# ;     printf,lun,'Clear'
# ;     printf,lun,'+RA:'+strtrim(calcat.coords.RA*!pi/180.0,2)+' Dec:'+strtrim(calcat.coords.Dec*!pi/180.0,2)
# ;     print,'RA:'+strtrim(calcat.coords.RA*!pi/180.0,2)+' Dec:'+strtrim(calcat.coords.Dec*!pi/180.0,2)
# ;     free_lun,lun
  # endif

  # calcat=flux_limit(calcat,pars,ATs=ATs,baseline=baseline,FINITO=FINITO)
  # if keyword_set(max_diam) then begin
     # ix=where(calcat.diam.theta le max_diam) & if (ix[0] ne -1) then calcat=calcat[ix] else begin
        # print,'No calibrators fullfililng the selection criteria with diameter <='+strtrim(string(max_diam,format='(f20.3)'),2)+' mas'
        # stop
     # endelse
  # endif
  
  # julian=systime(/julian,/utc)                      ;; current julian date in Greenwich mean time
  # julian=julian+pars.clock_correction/(24d0*3600d0) ;; small correction because the laptop's clock is not running exactly on time

  # target_list=insert_zenith(target_list,lst,pars)
  # if keyword_set(zenith) then target='zenith'

  # airmass=get_airmass(target,calcat,target_list,lst,pars,julian,duration=duration,do_cal_before=do_cal_before,offset_cal=offset_cal)
  # find_suitable_calibrators,target,calcat,target_list,airmass,pars,lst,julian,$
                            # cal=cal,baseline=baseline,delrestr=delrestr,plot_airmasses=plot_airmasses,zoom=zoom,print_plan=print_plan,$
                            # minF10=minF10,do_cal_before=do_cal_before,LN=LN,airmassPlotParams=airmassPlotParams,noPlot_restrict=noPlot_restrict,$
                            # duration=duration,offset_cal=offset_cal,galactic_plane=galactic_plane


  # if keyword_set(halt) then stop
# end


if __name__ == "__main__":
    alt, az = get_delay_line_restrictions("small", "paranal")
    print(alt, az, sep="\n\n")
    get_catalogue()

