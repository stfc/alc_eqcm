!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module containing bibliography information of  the implemented XC 
! approximations and vdW corrections
! corrections
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author    - i.scivetti July 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module references 

  Implicit None

  ! DFT-XC references 
  Character(Len=256), Parameter, Public  :: bib_slater= '[Slater, J.C. Phys. Rev., 81, 385-390 (1951)]'
  Character(Len=256), Parameter, Public  :: bib_ca=     '[D. M. Ceperley et al. Phys. Rev. Lett., 45:566–569, Aug 1980]'
  Character(Len=256), Parameter, Public  :: bib_hl=     '[L Hedin and B I Lundqvist. J. of Phys. C: Solid State Physics,&
                                                        & 4(14):2064–2083, Oct 1971]'
  Character(Len=256), Parameter, Public  :: bib_pz=     '[J. P. Perdew and A. Zunger. Phys. Rev. B, 23:5048–5079, May 1981]' 
  Character(Len=256), Parameter, Public  :: bib_wigner= '[E. Wigner. Phys. Rev., 46:1002–1011, Dec 1934]'  
  Character(Len=256), Parameter, Public  :: bib_vwn=    '[S. H. Vosko, L. Wilk, and M. Nusair. Canadian Journal of Physics,&
                                                        & 58(8):1200–1211 (1980);&
                                                        & S. H. Vosko and L. Wilk. Phys. Rev. B, 22:3812–3815, Oct 1980.]'
  Character(Len=256), Parameter, Public  :: bib_pade=   '[S. Goedecker, M. Teter, and J. Hutter.&
                                                        & Phys. Rev. B, 54:1703–1710, Jul 1996]'
  Character(Len=256), Parameter, Public  :: bib_pw92=   '[J. P. Perdew and Y. Wang. Phys. Rev. B, 45:13244–13249, Jun 1992]'
  Character(Len=256), Parameter, Public  :: bib_pw91=   '[J. P. Perdew et al. Phys. Rev. B, 46:6671–6687, Sep 1992]'
  Character(Len=256), Parameter, Public  :: bib_am05=   '[R. Armiento and A. E. Mattsson. Phys. Rev. B, 72:085108, Aug 2005]'
  Character(Len=256), Parameter, Public  :: bib_pbe=    '[J. P. Perdew, K. Burke, and M. Ernzerhof.&
                                                        & Phys. Rev. Lett., 77:3865–3868, Oct 1996]' 
  Character(Len=256), Parameter, Public  :: bib_rp=     '[B. Hammer, L. B. Hansen, and J. K. Norskov.&
                                                        & Phys. Rev. B, 59:7413–7421, Mar 1999]'
  Character(Len=256), Parameter, Public  :: bib_revpbe= '[Y. Zhang and W. Yang. Phys. Rev. Lett., 80:890–890, Jan 1998]'
  Character(Len=256), Parameter, Public  :: bib_pbesol= '[G. I. Csonka et al. Phys. Rev. B, 79:155107, Apr 2009]'
  Character(Len=256), Parameter, Public  :: bib_blyp=   '[A. D. Becke. Phys. Rev. A, 38:3098–3100, Sep 1988;&
                                                        & C. Lee, W. Yang, and R. G. Parr. Phys. Rev. B, 37:785–789, Jan 1988]'
  Character(Len=256), Parameter, Public  :: bib_xlyp=   '[X. Xu and W. A. Goddard. PNAS, 101(9):2673–2677, 2004]'
  Character(Len=256), Parameter, Public  :: bib_wc=     '[Z. Wu and R. E. Cohen. Rev. B, 73:235116, Jun 2006]'

  ! DFT-vdW references
  Character(Len=256), Parameter, Public  :: bib_g06=     '[S. Grimme. J. Computational Chemistry, 27(15):1787–1799 (2006)]' 
  Character(Len=256), Parameter, Public  :: bib_obs=     '[F. Ortmann, F. Bechstedt, and W. G. Schmidt.&
                                                         & Phys. Rev. B, 73:205101, May 2006.]' 
  Character(Len=256), Parameter, Public  :: bib_jchs=    '[P. Jurecka et al. J Computational Chemistry, 28(2):555–569, (2007)]'
  Character(Len=256), Parameter, Public  :: bib_dftd2=   '[S. Grimme. J. Computational Chemistry, 27(15):1787–1799 (2006)]'
  Character(Len=256), Parameter, Public  :: bib_dftd3=   '[S. Grimme et al. JCP, 132(15):154104 (2010)]'
  Character(Len=256), Parameter, Public  :: bib_dftd3bj= '[S. Grimme et al. J. Computational Chemistry, 32(7):1456–1465 (2011)]'
  Character(Len=256), Parameter, Public  :: bib_ts=      '[A. Tkatchenko and M. Scheffler. Phys. Rev. Lett., 102:073005, Feb 2009]'
  Character(Len=256), Parameter, Public  :: bib_tsh=     '[Tomas Bucko et al. JCP, 141(3):034114 (2014)]'
  Character(Len=256), Parameter, Public  :: bib_mbd=     '[A. Tkatchenko et al. Phys. Rev. Lett., 108:236402, Jun 2012;&
                                                         & A. Ambrosetti et al. JCP, 140(18):18A508 (2014)]'
  Character(Len=256), Parameter, Public  :: bib_ddsc=    '[S. N. Steinmann and C. Corminboeuf. JCTC, 7(11):3567–3577 (2011)]'
  Character(Len=256), Parameter, Public  :: bib_vdwdf=   '[M. Dion et al. Phys. Rev. Lett. (2004)]'
  Character(Len=256), Parameter, Public  :: bib_optpbe=  '[J. Klimes et al. J. Phys. Cond. Mat. 22 (2010)]'
  Character(Len=256), Parameter, Public  :: bib_optb88=  '[J. Klimes et al. J. Phys. Cond. Mat. 22 (2010)]'
  Character(Len=256), Parameter, Public  :: bib_optb86b= '[J. Klimes et al. Phys. Rev. B 83, 195131 (2011)]'
  Character(Len=256), Parameter, Public  :: bib_vdwdf2=  '[K. Lee et al. Phys. Rev. B 82, 081101 (2010)]'
  Character(Len=256), Parameter, Public  :: bib_AVV10S=  '[T. Bjorkman. Phys. Rev. B (2012)]'
  Character(Len=256), Parameter, Public  :: bib_vdwdf2b86r= '[I. Hamada. Phys. Rev. B 89, 121103 (2014)]'   
  Character(Len=256), Parameter, Public  :: bib_SCANrVV10 = '[H. Peng et al. Phys. Rev. X 6, 041005 (2016)]'
  Character(Len=256), Parameter, Public  :: bib_VV10=       '[O. A. Vydrov et al. J. Chem. Phys. (2010);&
                                                            & R. Sabatini et al. Phys. Rev. B (2013)]'

  Character(Len=256), Parameter, Public  :: bib_tunega=   '[D. Tunega et al. J. Chem. Phys. 137, 114105 (2012)]'
  Character(Len=256), Parameter, Public  :: bib_rpw86 =   '[E. D. Murray. JCTC (2009), 5, 10, 2754–2762]'
  Character(Len=256), Parameter, Public  :: bib_fisher=   '[M. Fisher and R. S. Angel. J. Chem. Phys. 146, 174111 (2017)]' 
  Character(Len=256), Parameter, Public  :: bib_scan=     '[J. Sun et al. Phys. Rev. Lett. 115, 036402 (2015)]' 
  Character(Len=256), Parameter, Public  :: bib_rpw86pbe= '[E. D. Murray et al. J. Chem. Theory Comput. (2009), 5, 10,&
                                                          & 2754–2762]'
  Character(Len=256), Parameter, Public  :: web_D3BJ=    '[https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/&
                                                          &dft-d3/functionalsbj]' 
  Character(Len=256), Parameter, Public  :: web_vasp   = 'https://www.vasp.at/'
  Character(Len=256), Parameter, Public  :: web_cp2k   = 'https://manual.cp2k.org/'
  Character(Len=256), Parameter, Public  :: web_castep = 'http://www.castep.org/'
  Character(Len=256), Parameter, Public  :: web_onetep = 'https://www.onetep.org/'
  Character(Len=256), Parameter, Public  :: web_slurm  = 'https://slurm.schedmd.com/documentation.html'

End Module references
