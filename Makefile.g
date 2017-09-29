O = obj-ifort/
S = src/

DEFINES = -DBRDF_COUPLING -DDISORT_BOUND_SRC -DDISORT_HORIZON -DSEA_BRDF -DDIREM_OPTIMIZATION
LINK       = ifort -g -O -o
COMPILE    = ifort -O0 -g -i4 -fpp $(DEFINES)   -c
EXECUTABLE = Mod5300alpha-ifort.exe

$(O)%.o: $(S)%.f
	$(COMPILE) $< -o $@

$(O)%.o: $(S)%.f90
	$(COMPILE) $< -o $@

OBJS := $(O)Modtrn.o $(O)cirmod.o \
	$(O)cldmod.o $(O)abcdta.o \
	$(O)aerext.o $(O)aermrg.o \
	$(O)aernsm.o $(O)aerprf.o \
	$(O)aertmp.o $(O)aitk.o   \
	$(O)albopn.o $(O)albtrn.o \
	$(O)altrin.o $(O)aprfnu.o \
	$(O)aruexa.o $(O)asymtx.o \
	$(O)atmcon.o $(O)bbfn.o   \
	$(O)betabs.o $(O)bfh2o.o  \
	$(O)binwrt.o $(O)bmcork.o \
	$(O)bmdata.o $(O)bmflux.o \
	$(O)bmload.o $(O)bmod.o   \
	$(O)bmtran.o $(O)bmxld.o  \
	$(O)bo3hh0.o $(O)bo3hh1.o \
	$(O)bo3hh2.o $(O)bs.o     \
	$(O)c6dta.o  $(O)card5.o  \
	$(O)cd1.o    $(O)cd1a.o   \
	$(O)cd2.o    $(O)cd2a.o   \
	$(O)cd2b.o   $(O)cd2c.o   \
	$(O)cd2e.o   $(O)cd3.o    \
	$(O)cd3a.o   $(O)cd3b3c.o \
	$(O)cd4.o    $(O)cd5.o    \
	$(O)check.o  $(O)chekin.o \
	$(O)chkres.o $(O)cirr18.o \
	$(O)cirrus.o $(O)ckd.o    \
	$(O)clddta.o $(O)cldprf.o \
	$(O)cmpint.o $(O)cnvphi.o \
	$(O)compar.o $(O)cool.o   \
	$(O)cool0.o  $(O)cossct.o \
	$(O)cph2o.o  $(O)cpo3.o   \
	$(O)cptrcg.o $(O)cpumix.o \
	$(O)crdeal.o $(O)crdriv.o \
	$(O)crfile.o $(O)crmerg.o \
	$(O)crprof.o $(O)crspec.o \
	$(O)cruprf.o $(O)cruspc.o \
	$(O)cxdta.o  $(O)dbltx.o  \
	$(O)debye.o  $(O)defalt.o \
	$(O)denlay.o $(O)desatt.o \
	$(O)devcbd.o $(O)dfltp.o  \
	$(O)dgrd.o   $(O)disbbf.o \
	$(O)disort.o $(O)dpandx.o \
	$(O)dpfish.o $(O)dprarf.o \
	$(O)dpscht.o $(O)driver.o \
	$(O)dstdta.o $(O)errmsg.o \
	$(O)exabin.o $(O)expint.o \
	$(O)extdt1.o $(O)extdt2.o \
	$(O)facfnc.o $(O)flclos.o \
	$(O)fdbeta.o $(O)fill.o   \
	$(O)filter.o $(O)filtsm.o \
	$(O)findmn.o $(O)flayz.o  \
	$(O)fluxes.o $(O)flxadd.o \
	$(O)flxsum.o $(O)flxtbl.o \
	$(O)fnames.o $(O)fndpth.o \
	$(O)frn296.o $(O)frq5dt.o \
	$(O)ftrang.o $(O)gamfog.o \
	$(O)geo.o    $(O)geodrv.o \
	$(O)geoinp.o $(O)geopth.o \
	$(O)getpf.o  $(O)getvis.o \
	$(O)gmrain.o $(O)gtcldv.o \
	$(O)gtpsap.o $(O)gtstrt.o \
	$(O)h2h2dt.o $(O)h2hedt.o $(O)hengns.o \
	$(O)hertda.o $(O)hno3.o   \
	$(O)indxrf.o $(O)intcor.o \
	$(O)initcd.o \
	$(O)interp.o $(O)irecln.o \
	$(O)irfxn.o  $(O)isamax.o \
	$(O)jou.o    $(O)jubran.o \
	$(O)laycld.o \
	$(O)layvsa.o $(O)lenstr.o \
	$(O)lepoly.o $(O)locate.o \
	$(O)loop.o   $(O)lrdsap.o \
	$(O)m3d.o \
	$(O)m3ddb.o  $(O)m3din.o  \
	$(O)m3drew.o $(O)m3dwrt.o \
	$(O)mapms.o  $(O)mardta.o \
	$(O)marine.o $(O)mccont.o \
	$(O)mcmol.o  $(O)mdta.o   \
	$(O)mlatmb.o $(O)molsct.o \
	$(O)msrad.o  $(O)mssolr.o \
	$(O)n2data.o \
	$(O)newh2.o  $(O)newsrc.o \
	$(O)no2xs.o  $(O)novaer.o \
	$(O)novmrg.o $(O)nunit.o  \
	$(O)o2cont.o $(O)o2inf2.o \
	$(O)o2mate.o $(O)o3chap.o \
	$(O)o3hht0.o $(O)o3hht1.o \
	$(O)o3hht2.o $(O)o3int.o  \
	$(O)o3uv.o   $(O)openbm.o \
	$(O)opnfl.o  $(O)pf.o     \
	$(O)pfexp.o  $(O)phasef.o \
	$(O)phsdta.o $(O)plkavg.o \
	$(O)prfdta.o \
	$(O)prtinp.o $(O)prtint.o \
	$(O)psieca.o $(O)pslct.o  \
	$(O)qgausn.o $(O)rain.o   \
	$(O)ratio.o  $(O)rdcork.o \
	$(O)rdexa.o  $(O)rdnsm.o  \
	$(O)rdsun.o  $(O)reduce.o \
	$(O)rfgeom.o $(O)rflect.o \
	$(O)rfpath.o $(O)rfract.o \
	$(O)right.o  $(O)rmchan.o \
	$(O)rnscat.o $(O)rtbis.o  \
	$(O)s15bd.o  $(O)sasum.o  \
	$(O)saxpy.o  $(O)schrun.o \
	$(O)sclcol.o $(O)scnflx.o \
	$(O)sdot.o   $(O)secsca.o \
	$(O)segmnt.o $(O)setdis.o \
	$(O)setmtx.o $(O)sf260.o  \
	$(O)sf296.o  $(O)sgbco.o  \
	$(O)sgbfa.o  $(O)sgbsl.o  \
	$(O)sgeco.o  $(O)sgefa.o  \
	$(O)sgesl.o  $(O)shade.o  \
	$(O)sinsca.o $(O)sint.o   \
	$(O)slf260.o $(O)slf296.o \
	$(O)slftst.o $(O)smgeo.o  \
	$(O)smprep.o $(O)sncms.o  \
	$(O)so2xs.o  $(O)solazm.o \
	$(O)soleig.o $(O)solve0.o \
	$(O)solve1.o $(O)solzen.o \
	$(O)source.o $(O)spaltr.o \
	$(O)spcflx.o \
	$(O)sscal.o  $(O)sscork.o \
	$(O)ssgeo.o  $(O)ssrad.o  \
	$(O)stdmdl.o $(O)subsol.o \
	$(O)svsola.o \
	$(O)tab.o    $(O)tanht.o  \
	$(O)terpev.o $(O)terpso.o \
	$(O)title.o  $(O)tnrain.o \
	$(O)trans.o  $(O)trlay.o  \
	$(O)tstbad.o $(O)upbeam.o \
	$(O)upcase.o $(O)upisot.o \
	$(O)usrint.o $(O)vsa.o    \
	$(O)vsansm.o $(O)watvap.o \
	$(O)wrtalb.o $(O)wrtbad.o \
	$(O)wrtbuf.o $(O)wrtdim.o \
	$(O)wrtflt.o $(O)wtchan.o \
	$(O)wtcool.o $(O)wtlft.o  \
	$(O)wtrgt.o  $(O)wtsum.o  \
	$(O)xifunc.o $(O)xmlatm.o \
	$(O)zensun.o $(O)zeroal.o \
	$(O)zeroit.o \
	$(O)bright.o $(O)alb_up.o $(O)o2_o2.o $(O)h1_vec.o $(O)ch4_lo.o \
	$(O)ch4_hi.o $(O)d_loop.o $(O)i_loop.o $(O)t_loop.o \
	$(O)seg5.o $(O)seg5ms.o \
	$(O)bmod0.o $(O)loop0.o $(O)doslit.o $(O)wthead.o $(O)wrtplt.o \
	$(O)sbrdf.o

$(EXECUTABLE): $(OBJS)
	$(LINK) $(EXECUTABLE) $(OBJS)

clean:
	rm $(OBJS) *.mod
