
	 !1/25/05-yyh: for some reason when only TE, program no longer calculate TM fields, change it back to if only TE is launched both TE and TM variables are calculated 
	 !1/11/05-yyh: modify Ey_inc in "Total Field/Scattered field for TE-mode (Hz)"
	 !10/28/04-yyh: add 4 level medium
	 !4/1/03-ssc: up to 16 computers possible
	 !3/25/03-ssc: metal option is added using sigma. 
	 !	Note that the format of "indexprofile.csv" changed.
	 !1/23/02-ssc:nei deallocates after usage:
	 !	major memory: geometry(1), detector(1), (Ez,Hx,Hy)(4) : 14*ie*je Byte
	 !12/20/02-ssc:integrated spectrum detector feature added:
	 !	detector can be used for spectrum detector (specify in "indexprofile.csv") 
	 !12/12/02-ssc:Added several IO options (specify in "userinput.csv": Diagnostic output,
	 ! flux-map output, E(t) ouput (has to be 1 for spectrum analysis), input waveguide position
	 !12/12/02-ssc:spectrum analysis working
	 !12/03/02-ssc
	 !One can specify detector to be Sx or Sy.
	 !need to specify detector:1 for Sx and 2 for Sy in the "indexProfile.csv".
	 !7/xx/02-Ravil
	 !MPICH is used for parallel computing
	 !5/11/02-SSC
	 !If only TE variables (Hz, Ex, Ey) are calculated, high intensity build up from the corner->
	 !Even if only TE is launched both TE and TM variables are calculated, but output is selective 
	 !5/1/02-SSC
	 !Spectrum analysis updated
	 !Flux (Sx) Detector by YY 
	 !UPML 2-D TM & TE case: (Ez, Hx, Hy), (Hz, Ex, Ey) by SSC 07/2001
	 !Slab UPML code --------
	 !4-level model gain material without touching the boundary
	 !fixed dx, independent of wavelength
	 !Atom-field interaction (index:interaction) (1: dE, 2: Ap)
	 !Pauli exclusion included (index:ipauli) (0:no, 1: stop nr-decay, 2: N3(1-N2), N1(1-N0))
	 !Automatically detect waveguide positions (3/15/01)
	 !Three waveguides are possible	 (3/22/01)
	 !Two launching fields are possible for each waveguide  (3/22/01)
	 !TE or TM modes can be launchable (4/2/01)
	 !MPICH multiple processes option added (7/9/02)
 
 !------------------------------------------------------------
 !Required input files:
 !userinput.csv	: input parameters
 !graphic4.pgm	: color-coded waveguide geometry
 !indexProfile.csv : index, gain assignment for each color code
 !observerpoint.csv	 : observation points
 !------------------------------------------------------------

 !--------------------------------------------------------------
 !Output files:
 !Ez##.csv:	snapshot for field
 !ni##.csv:	snapshot for population
 !E-DFT:
 !E-DFT-Normal:
 !---------------------------------------------------------------
 !! 
 ! condition compilation variable USE_MPI_CH used for MPI libraries
 ! use !!DEC.. to build program without MPI; 
 !!DEC$ DEFINE USE_MPI_CH
 !DEC$ DEFINE DBG_FILE_WRITE	
  program fdtd2d
 	USE DFPORT  ! for date function FDATE

!DEC$ IF DEFINED (USE_MPI_CH) 
	use mpi_exhchg
!DEC$ ENDIF ! (USE_MPI_CH) 

   IMPLICIT NONE

!DEC$ IF DEFINED (USE_MPI_CH) 
	include "mpif.h"
!DEC$ ENDIF ! (USE_MPI_CH) 

	character*256 nothing
	real tmax,tdelay0,tdelay,ntspan,epsnclad,xx_old,cent,xlow,xhigh,vallow,valhigh
	real tt0,sinegauss,epsn,delbceta,x1,x2,ffreq,delbc,eta
	real DSbcf,Dsbcr,BHZSbcf,BHZSbcr,DSbcb,BHZSbcb,BHXSbcf,DEXSbcf,BHXSbcb,DEXSbcb,BHYSbcf,DEYSbcf
	real Dsbcl,BHZSbcl,BHYSbcb,DEYSbcb,BHXSbcl,DEXSbcl,BHYSbcl,DEYSbcl,BHXSbcr,DEXSbcr,BHYSbcr,DEYSbcr
	real BHXSbcf1,BHXSbcf2,DEXSbcf1,DEXSbcf2,BHXSbcb1,BHXSbcb2,DEXSbcb1,DEXSbcb2
	real epsrfij,epsrbij,Hy_inc,Hx_inc,Ey_inc,Ex_inc,Ex_inc250,Ez_inc,Hz_inc
	real a1,ez01
	integer il5
	integer s_x_f_d_fid,s_y_f_d_fid,s_f_d_fid
	integer nfile11,nfile12,nfile13,nfile1,nfile2,nfile3,nfile4

	integer nsnapshot,tempdata,nmax,Nwgmax,Nfieldmax,num_index,numfile2
	integer m,inumber,jnumber,ie,je,ib,jb,ibbc,jbbc,iefbc,ibfbc,jefbc,j1_al,j2_al
	integer iff,k,ii,jj,jws1,jwe1,jmax,numfile,if1,if2,il1,ir,jf,ill
	integer num_wg,wg_finder,num_wgx,num_wgy,maxgrid,ilaunchmk ,id,jd
	integer js_np,count
	
	
	real muz,muer,Pi,cc,epsz,epsr,tpi,rmax,rnf,rns,dx,dy,dt, &
          sinusoidal,gaussian,beamshape,dt_over_dx,sigmam, &
		  sigma,bcfactor,kapa,kmax,kfactor,tt, &
		  s_x_f_d_max, s_y_f_d_max, s_f_d_max, geometryij, sigmageoij,  &
		  exmax,exmin,ix,ixmin,ixmax,logix,logixmin,logixmax,&
		  eymax,eymin,iy,iymin,iymax,logiy,logiymin,logiymax, &
		  ezmax,ezmin,iz,izmin,izmax,logiz,logizmin,logizmax

	real temp1,temp2,temp3,temp4	!tempoary storage
	real  wgwidth,wgwidth2,epsnwg,a0,kk,kclad,Fx,Gx,x,xx
	real*8 lamdastart,lamdaend 
	real*8	h_bar,h_bar_inv,mass_e 
	integer ipauli, interaction, TM, TE	, det_num, ftdet_num, metal, gain, gainij, p

	!----gain related
	real,pointer::  lamda_a1(:),lamda_a2(:), wa1(:),wa2(:),dwa1(:),&
		  dwa2(:),ne(:),n0_0(:),n1_0(:),n2_0(:),n3_0(:), current_pumping(:),pumping_rate(:), &
		  Ndensity(:),tau30(:), tau21(:), tau10(:), tau32(:), tau3(:), &
		  dt_over_tau30(:), dt_over_tau21(:), dt_over_tau10(:), dt_over_tau32(:), dt_over_tau3(:)
	integer observer_step, observer_counter, atom_gridsize, atomgroup_number
				
 !---------------------------------------------------------
 !		Fundamental Constants
 !---------------------------------------------------------
    parameter(pi=3.141592654)
	parameter(cc=3.000000000e+8)
	parameter(muz=4.0*pi*1.0e-7)
	parameter(epsz=1./(muz*cc*cc))
	parameter(epsr=1.0/epsz,muer=1.0/muz)
	parameter(tpi=2.0*pi)
	parameter(h_bar=1.05459e-34)
	parameter(h_bar_inv=1./h_bar)
	parameter(mass_e=9.108e-31) 

 !----------------------------------------------------------
 !       UPML parameters
 !----------------------------------------------------------
  	integer iebc,jebc
	parameter(iebc=16,jebc=16)	!UPML grid size
  	integer orderbc,mediabc
    parameter(rmax=1.e-7,orderbc=3,mediabc=2) !Polynomial grading UPML: Rmax, m, 
	 	 
 !---------------------------------------------------------- 
 !		other paramter
 !----------------------------------------------------------
    parameter(rns=1.e0) !rns = free space index of refraction 

!------------------------------------------------------------
!.................Grid Arrays...............................
!------------------------------------------------------------

! fields and fluxes
	real, pointer:: Ez(:,:),Hx(:,:),Hy(:,:),Hz(:,:),Ex(:,:),Ey(:,:)
	real, pointer:: Ex_obs(:,:), Ey_obs(:,:), Ez_obs(:,:) 
	real, pointer:: s_x_f_d(:,:), s_y_f_d(:,:), s_f_d(:,:)
	real, pointer:: s_x_f_d_r(:,:), s_y_f_d_r(:,:), s_f_d_r(:,:)

!  media index	 
  real, pointer:: eps(:),mur(:)		!relative epsilon & mu
  real  ca,da,cb,db,caij,cbij

!------------GAIN RELATED ARRAYS--------------------------  
	!Vector potential
	real, pointer:: Avptz(:,:),Avptz_old(:,:)	
	real, pointer:: Avptx(:,:),Avptx_old(:,:)	
	real, pointer:: Avpty(:,:),Avpty_old(:,:)	

	!polarization and atomic population
	real, Pointer:: pz_1(:,:,:),pz_old1_1(:,:,:),pz_2(:,:,:),pz_old1_2(:,:,:)
	real, Pointer:: px_1(:,:,:),px_old1_1(:,:,:),px_2(:,:,:),px_old1_2(:,:,:)
	real, Pointer:: py_1(:,:,:),py_old1_1(:,:,:),py_2(:,:,:),py_old1_2(:,:,:)
	real, Pointer:: n0(:,:,:),n1(:,:,:),n2(:,:,:),n3(:,:,:)
	real, Pointer:: n0_result(:,:),n1_result(:,:),n2_result(:,:),n3_result(:,:)
	!storage and intermediate array
	real, Pointer:: ez_old(:,:),ex_old(:,:),ey_old(:,:)
	real, Pointer:: n1_old(:,:,:),n2_old(:,:,:),n3_old(:,:,:),n0_old(:,:,:)
	real, Pointer:: n1_obs(:,:),n2_obs(:,:),n3_obs(:,:),n0_obs(:,:)
	real, Pointer::	n0ij(:),n1ij(:),n2ij(:),n3ij(:)
	real, Pointer:: Pz_old2_1(:),Pz_old2_2(:),Px_old2_1(:),Px_old2_2(:),Py_old2_1(:),Py_old2_2(:)
	real, pointer::	Pxg_1ij(:),Pxg_old1_1ij(:),Pxg_2ij(:),Pxg_old1_2ij(:),Pyg_1ij(:),Pyg_old1_1ij(:), &
				Pyg_2ij(:),Pyg_old1_2ij(:),EdP2ij(:),EdP1ij(:),	AP1ij(:), AP2ij(:)
	real exgij, exg_oldij, eygij, eyg_oldij, AVPTXGIJ, AVPTXG_OLDIJ,AVPTyGIJ, AVPTyG_OLDIJ
	!coefficients related to gain 
	real,pointer:: pa1_1(:),pa1_2(:),pa2_1(:),pa2_2(:),pa3_1(:),pa3_2(:),kapa1(:),kapa2(:)
	real*8,pointer:: C_n2_1(:),C_n2_2(:),C_n1_1(:),C_n1_2(:)
	real*8,pointer:: rabico1(:),rabico2(:) 
!-------------END OF GAIN RELATED ARRAYS-------------------

!TM-mode------------------
! UPML grid for Ez, Hx,Hy,epsr (f:front, b:back, l:left, r:right)
  real,pointer:: Ezf(:,:),Ezb(:,:),Ezbcl(:,:),Ezbcr(:,:)
	real,pointer:: epsrf(:,:),epsrb(:,:),epsrbcl(:,:),epsrbcr(:,:) !epsilon
  real,pointer:: Hxf(:,:),Hxb(:,:),Hxbcl(:,:),Hxbcr(:,:),  &
	               Hyf(:,:),Hyb(:,:),Hybcl(:,:),Hybcr(:,:)
! UPML grid for Dz, Bx, and By  
  real,pointer:: DEzf(:,:),DEzb(:,:),DEzbcl(:,:),DEzbcr(:,:)  	    	       
	real,pointer:: BHxf(:,:),BHxb(:,:),BHxbcl(:,:),BHxbcr(:,:),	 &
	               BHyf(:,:),BHyb(:,:),BHybcl(:,:),BHybcr(:,:)
!TE-mode------------------
! UPML grid for Hz, Ex,Ey (f:front, b:back, l:left, r:right)
    real,pointer:: Hzf(:,:),Hzb(:,:),Hzbcl(:,:),Hzbcr(:,:)
    real,pointer:: Exf(:,:),Exb(:,:),Exbcl(:,:),Exbcr(:,:),  &
	               Eyf(:,:),Eyb(:,:),Eybcl(:,:),Eybcr(:,:)
! UPML grid for Hz, Dx, and Dy  
  real,pointer:: BHzf(:,:),BHzb(:,:),BHzbcl(:,:),BHzbcr(:,:)  	    	       
	real,pointer:: DExf(:,:),DExb(:,:),DExbcl(:,:),DExbcr(:,:),	 &
	               DEyf(:,:),DEyb(:,:),DEybcl(:,:),DEybcr(:,:)

! Coefficients for Maxwell's eq. in UPML grid 
	real,pointer:: axef(:) ,ayef(:) ,axhf(:) ,ayhf(:), &
	               axeb(:) ,ayeb(:) ,axhb(:) ,ayhb(:), &
			       axel(:) ,ayel(:) ,axhl(:) ,ayhl(:), &
			       axer(:) ,ayer(:) ,axhr(:) ,ayhr(:), &
	               bxef(:) ,byef(:) ,bxhf(:) ,byhf(:), &
			       bxeb(:) ,byeb(:) ,bxhb(:) ,byhb(:), &
			       bxel(:) ,byel(:) ,bxhl(:) ,byhl(:), &
			       bxer(:) ,byer(:) ,bxhr(:) ,byhr(:), &
	               cxef(:) ,cyef(:) ,cxhf(:) ,cyhf(:), &
			       cxeb(:) ,cyeb(:) ,cxhb(:) ,cyhb(:), &
			       cxel(:) ,cyel(:) ,cxhl(:) ,cyhl(:), &
			       cxer(:) ,cyer(:) ,cxhr(:) ,cyhr(:)
! General bc coefficient Array
	real axebc(iebc,mediabc),ayebc(jebc,mediabc)
	real bxebc(iebc,mediabc),byebc(jebc,mediabc)
	real cxebc(iebc,mediabc),cyebc(jebc,mediabc)
	real axhbc(iebc,mediabc),ayhbc(jebc,mediabc)
	real bxhbc(iebc,mediabc),byhbc(jebc,mediabc)
	real cxhbc(iebc,mediabc),cyhbc(jebc,mediabc)

! waveguide parameters
	!wavguide color-code, index, gain, parameters
	!!
	!Variables below are labeled with index number
	integer,pointer:: i_color(:),det_index(:),spec_index(:),n_skip(:),gain_index(:)
	real,pointer:: n_index(:), sigma_index(:), ni_index(:)
	!below detector related variables are labeled with detector number
	integer, pointer:: det_dir(:),det_xp(:), det_yp(:), det_xlen(:), det_ylen(:)
	integer det_maxlen 
	!index, gain grid (geometry = 1/n^2)
	real,pointer:: geometry(:,:), sigmageo(:,:)
	!!
	integer,pointer:: detector(:,:), gaingeo(:,:)
	real,pointer:: det_result(:) ! sx flux detector  

	!waveguide geometry parameters
	integer wgstart(4),wgend(4),width(4),wgdirection(4)
	real wgindex(4),cladindex(4)
	real, pointer:: rmode(:,:,:)
	character*32 rmodefilename(4,4)

! input field parameters
	real lamda_f(4,4),w_f(4,4),pwidth(4,4),k_f(4,4), kprop(4,4)
	real E_input(4,4),I_input_cgs(4,4),I_input(4,4),  &
		 tdel(4,4)
	integer nfield(4),inputfield(4,4),Ilaunch(4,4),pol(4,4)
		
! observer	
	integer,pointer:: m_obs(:)
	real,pointer:: iobs(:),jobs(:)
	real,pointer:: tt_obs(:)
	integer nobserver, ns_obs, ifmax 
!!
	real det_period,det_time
! snapshot output storage
	character*32 afilename(0:99)
	!medium
	character*32 bfilename(0:99) 
	character*32 cfilename(0:99) 
	character*32 dfilename(0:99) 
	character*32 efilename(0:99)
	!end medium  
	character*32 ffilename(0:99) 
	character*32 gfilename(0:99) 
	character*12 datadir, rootdir	! data directory name

	CHARACTER*24 timestart, timeend

	real snapshotinterval
	integer reso_snapshot_x, reso_snapshot_y, nsnapshotinterval

! Fourier Transform arrays
	real*8, pointer:: ftr_in(:),fti_in(:),fte2_in(:),  &
	ftr_z(:,:),fti_z(:,:),fte2_z(:,:),fte2max_z(:),  &
	ftr_x(:,:),fti_x(:,:),fte2_x(:,:),fte2max_x(:),  &
	ftr_y(:,:),fti_y(:,:),fte2_y(:,:),fte2max_y(:),  &
	ezinc(:), exinc(:), eyinc(:)
	real ftarg,ftsin,ftcos,wstart,wspan, wavelength0,fte2max_in
	integer spectrumanalysis, ndft_input
	real*8, pointer:: ftrd_z(:,:,:),ftid_z(:,:,:),fte2d_z(:,:,:),fte2d_zt(:,:,:),fte2maxd_z(:), &
	ftrd_x(:,:,:),ftid_x(:,:,:),fte2d_x(:,:,:), fte2d_xt(:,:,:),fte2maxd_x(:), &
	ftrd_y(:,:,:),ftid_y(:,:,:),fte2d_y(:,:,:),fte2d_yt(:,:,:),fte2maxd_y(:)

! Other option parameters
	integer outputchoice, period, read_bpm, read_mpi_seg, diag, fluxmap
	integer inputwg, j_inputwg, n_fluxmap
	real t_fluxmap

! pre run spectrum analysis
	character*40 sprfilename(4,4)
	character*32 dftfilename(4,4)
!MPI parameters
	character*32 processor_name
	character*32 fmt_read

!reading bitmap file
 	integer a(75) 
  character*1 CH(75),test,space
	character*3 c3
  integer Number(4)
	integer,pointer:: nei(:,:)
	integer isize,jsize				 
	integer	n, n_st, num_ts_fd, n_fd, j_ffp
	integer u, v, w, rmf_ix1, rmf_ix2
	integer i,j,i1,i2,i1c,i2c,j1c,j2c,js1,js2,is1,is2,j1b,i1b

!DEC$ IF DEFINED (USE_MPI_CH) 
		integer namelen
		integer je_np(0:16), dje_np ! index boundaries arrays for mpich .._np - number of process
		integer pid, numprocs ! this process's id, number of processes
		integer ierr ! error status 
		integer status(MPI_STATUS_SIZE), req ! array of status values
		integer comm_a, comm_b	! communicator id's 
		integer pid_a, pid_b ! process's id in commucators '_a' and '_b' 
		integer pid_g ! n/a 
		integer master ! master process's id
		integer nbrbottom_a, nbrtop_a, nbrbottom_b, nbrtop_b
		real,pointer:: det_result_t(:) ! sx flux detector (_t - total, for mpi)
		integer,pointer:: obs_in_proc(:) ! arrays indicating which process's subgrid has
		integer num_tot									 ! this observer point (mpi) 
		integer maxdje
		real, pointer:: Ezr(:,:),Hzr(:,:),Eyr(:,:), n0r(:,:),n1r(:,:),n2r(:,:),n3r(:,:)
		integer nsend,nsendx,nsendy



!------------------------------------------------------------------------------------
!		end of vatiable declarations
!------------------------------------------------------------------------------------
!
		datadir = '.\..\_dat\'
		rootdir = '.\..\'
		master = 0 ! id of master process, fixed 

! initialize MPICH interface; 
	call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, pid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

	if( numprocs .le. 1 .or. numprocs .gt. 16 )then
		print *, " $223, number of mpich processes N = ", numprocs, ", 1 < N < 16 allowed  !"
		print *, " terminating program ..."
		call MPI_ABORT( MPI_COMM_WORLD, -229, ierr )
		call exit( 223 )
	endif 
	 
	call MPI_GET_PROCESSOR_NAME( processor_name, namelen, ierr )
	write(*,3300) "$229, process#:", pid, "of", numprocs, " is alive on ", processor_name
	
! create communicator 'a' for grid data exchange (1st time half step)
	call MPI_CART_CREATE( MPI_COMM_WORLD, 1, numprocs, .false., .true., comm_a, ierr ) 
! get my position in this communicator 
  call MPI_COMM_RANK( comm_a, pid_a, ierr ) 
  call MPI_Cart_shift( comm_a, 0,  1, nbrbottom_a, nbrtop_a, ierr   ) 

! create communicator 'b' for grid data exchange (1st time half step)
	call MPI_CART_CREATE( MPI_COMM_WORLD, 1, numprocs, .false., .true., comm_b, ierr ) 
! get my position in this communicator 
  call MPI_COMM_RANK( comm_b, pid_b, ierr ) 
  call MPI_Cart_shift( comm_b, 0,  1, nbrbottom_b, nbrtop_b, ierr   ) 
!!
3300 format( a17, i3, a4, i3, a14, a32 )

!DEC$ ELSE ! (USE_MPI_CH) 
		datadir = '.\_dat\'
		rootdir = '.\'
!DEC$ ENDIF ! (USE_MPI_CH) 

!-------------------------------------------------------------------------------
!	
!		 Read 'userinput.csv' 
!
!-------------------------------------------------------------------------------

	open(25, file= (trim(rootdir)//'userinput.csv'), action='read', share='denywr')
	read(25,'(a80)') nothing
	!FDTD parameters
	read(25,'(a80)') nothing
	read(25,*)       dx			   !x-grid spacing
	read(25,'(a80)') nothing
	read(25,*)       dy			   !y-grid spacing
	read(25,'(a80)') nothing
	read(25,*)       tmax		   !total simulation time
	read(25,'(a80)') nothing
	read(25,'(a80)') nothing
	
	!Atom-Field interaction
	read(25,'(a80)') nothing
	read(25,*)       interaction   !Atom-field interaction (1:dE, 2:Ap) 
	read(25,'(a80)') nothing
	read(25,*)       ipauli		   !Pauli-exclusion (0: no, 1: stop nr-decay, 2: n3(1-n2))
	read(25,'(a80)') nothing
	read(25,'(a80)') nothing  

	!Input Field parameters
	read(25,'(a80)') nothing
	read(25,*)       nfield(1)   !number of inputfield in wg1
	read(25,'(a80)') nothing
	read(25,*)       nfield(2)   !number of inputfield in wg2
	read(25,'(a80)') nothing
	read(25,*)       nfield(3)   !number of inputfield in wg3
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(1,1)   !inputfield wavelength(1,1) for wg1, field1
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(1,2)   !inputfield wavelength(1,2)	for wg1, field2
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(1,3)   !inputfield wavelength(1,3) for wg1, field3
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(1,4)   !inputfield wavelength(1,4)	for wg1, field4
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(2,1)   !inputfield wavelength(2,1)	for wg2, field1
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(2,2)   !inputfield wavelength(2,2)	for wg2, field2
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(2,3)   !inputfield wavelength(2,3)	for wg2, field3
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(2,4)   !inputfield wavelength(2,4)	for wg2, field4
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(3,1)   !inputfield wavelength(3,1)	for wg3, field1
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(3,2)   !inputfield wavelength(3,2)	for wg3, field2
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(3,3)   !inputfield wavelength(3,3)	for wg3, field3
	read(25,'(a80)') nothing
	read(25,*)       lamda_f(3,4)   !inputfield wavelength(3,4)	for wg3, field4
	read(25,'(a80)') nothing
	read(25,*)       pol(1,1)   !inputfield polarization(1,1): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(1,2)   !inputfield polarization(1,2): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(1,3)   !inputfield polarization(1,1): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(1,4)   !inputfield polarization(1,2): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(2,1)   !inputfield polarization(2,1): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(2,2)   !inputfield polarization(2,2): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(2,3)   !inputfield polarization(2,3): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(2,4)   !inputfield polarization(2,4): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(3,1)   !inputfield polarization(3,1): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(3,2)   !inputfield polarization(3,2): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(3,3)   !inputfield polarization(3,3): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*)       pol(3,4)   !inputfield polarization(3,4): TM:1, TE:0
	read(25,'(a80)') nothing
	read(25,*) 		 inputfield(1,1)	 ! inputfield(1) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing	 
	read(25,*) 		 inputfield(1,2)	 ! inputfield(2) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing
	read(25,*) 		 inputfield(1,3)	 ! inputfield(3) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing	 
	read(25,*) 		 inputfield(1,4)	 ! inputfield(4) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing
	read(25,*) 		 inputfield(2,1)	 ! inputfield(1) choice 1=pulse, othter=cw	
	read(25,'(a80)') nothing
	read(25,*) 		 inputfield(2,2)	 ! inputfield(2) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing
	read(25,*) 		 inputfield(2,3)	 ! inputfield(3) choice 1=pulse, othter=cw	
	read(25,'(a80)') nothing
	read(25,*) 		 inputfield(2,4)	 ! inputfield(4) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing	 
	read(25,*) 		 inputfield(3,1)	 ! inputfield(1) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing
	read(25,*) 		 inputfield(3,2)	 ! inputfield(2) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing	 
	read(25,*) 		 inputfield(3,3)	 ! inputfield(1) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing
	read(25,*) 		 inputfield(3,4)	 ! inputfield(2) choice 1=pulse, othter=cw
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(1,1)   !input intensity(1)
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(1,2) 
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(1,3)   !input intensity(1)
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(1,4) 
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(2,1)   !input intensity(2)
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(2,2) 
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(2,3)   !input intensity(2)
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(2,4) 
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(3,1)   !input intensity(3)
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(3,2) 
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(3,3)   !input intensity(3)
	read(25,'(a80)') nothing
	read(25,*)       I_input_cgs(3,4) 
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(1,1)		!pulse width (1)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(1,2)		!pulse width (2)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(1,3)		!pulse width (1)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(1,4)		!pulse width (2)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(2,1)		!pulse width (3)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(2,2)		!pulse width (1)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(2,3)		!pulse width (3)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(2,4)		!pulse width (1)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(3,1)		!pulse width (2)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(3,2)		!pulse width (3)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(3,3)		!pulse width (2)
	read(25,'(a80)') nothing
	read(25,*) 		 pwidth(3,4)		!pulse width (3)
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(1,1)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(1,2)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(1,3)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(1,4)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(2,1)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(2,2)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(2,3)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(2,4)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(3,1)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(3,2)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(3,3)		!delay
	read(25,'(a80)') nothing
	read(25,*) 		 tdel(3,4)		!delay
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(1,1)	 ! I(1,1) launch position + from left to right
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(1,2)	 ! I(1,2) launch position - from left to right
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(1,3)	 ! I(1,1) launch position + from left to right
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(1,4)	 ! I(1,2) launch position - from left to right
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(2,1)	 ! I(2,1) launch position 
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(2,2)	 ! I(2,2) launch position + from left to right
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(2,3)	 ! I(2,1) launch position 
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(2,4)	 ! I(2,2) launch position + from left to right
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(3,1)	 ! I(3,1) launch position - from left to right
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(3,2)	 ! I(3,2) launch position 
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(3,3)	 ! I(3,1) launch position - from left to right
	read(25,'(a80)') nothing
	read(25,*)		 Ilaunch(3,4)	 ! I(3,2) launch position 
	read(25,'(a80)') nothing
	read(25,'(a80)') nothing
	
	!gain parameters

	read(25,'(a80)') nothing
	read(25,*)	atomgroup_number
	read(25,'(a80)') nothing
	read(25,*)	atom_gridsize
	
	allocate (lamda_a1(atomgroup_number),lamda_a2(atomgroup_number), wa1(atomgroup_number), &
	wa2(atomgroup_number),dwa1(atomgroup_number),dwa2(atomgroup_number),ne(atomgroup_number), &
	n0_0(atomgroup_number),n1_0(atomgroup_number),current_pumping(atomgroup_number),pumping_rate(atomgroup_number),&
	n2_0(atomgroup_number),n3_0(atomgroup_number),Ndensity(atomgroup_number), tau30(atomgroup_number), &
	 tau21(atomgroup_number), tau10(atomgroup_number), tau32(atomgroup_number), tau3(atomgroup_number),&
	  dt_over_tau30(atomgroup_number), dt_over_tau21(atomgroup_number), dt_over_tau10(atomgroup_number), &
	  dt_over_tau32(atomgroup_number), dt_over_tau3(atomgroup_number))
	allocate (rabico1(atomgroup_number),rabico2(atomgroup_number),pa1_1(atomgroup_number), &
	pa1_2(atomgroup_number),pa2_1(atomgroup_number),pa2_2(atomgroup_number),pa3_1(atomgroup_number), &
	pa3_2(atomgroup_number),kapa1(atomgroup_number),kapa2(atomgroup_number))

	allocate (C_n2_1(atomgroup_number), &
	C_n2_2(atomgroup_number),C_n1_1(atomgroup_number),C_n1_2(atomgroup_number))

	allocate (Pz_old2_1(atomgroup_number),Pz_old2_2(atomgroup_number),Px_old2_1(atomgroup_number), &
	Px_old2_2(atomgroup_number),Py_old2_1(atomgroup_number),Py_old2_2(atomgroup_number), &
	n0ij(atomgroup_number),n1ij(atomgroup_number),n2ij(atomgroup_number), n3ij(atomgroup_number))

	allocate (Pxg_1ij(atomgroup_number),Pxg_old1_1ij(atomgroup_number),Pxg_2ij(atomgroup_number), &
	Pxg_old1_2ij(atomgroup_number),Pyg_1ij(atomgroup_number),Pyg_old1_1ij(atomgroup_number), &
	Pyg_2ij(atomgroup_number),Pyg_old1_2ij(atomgroup_number),EdP2ij(atomgroup_number),EdP1ij(atomgroup_number),&
	AP2ij(atomgroup_number),AP1ij(atomgroup_number))


do 	9876 i=1,atomgroup_number
	read(25,'(a80)') nothing
	read(25,'(a80)') nothing
	read(25,*)       lamda_a1(i)		!atomic resonance wavelegnth between 1&2
	read(25,'(a80)') nothing
	read(25,*)       lamda_a2(i)		!atomic resonance wavelength between 0&3
	read(25,'(a80)') nothing
	read(25,*)       tau30(i)			!radiative decay rate 3 -> 0
	read(25,'(a80)') nothing
	read(25,*)       tau21(i)			!radiative decay rate 2 -> 1
	read(25,'(a80)') nothing
	read(25,*)       tau10(i)			!non-radiative decay rate 1 -> 0
	read(25,'(a80)') nothing
	read(25,*)       tau32(i)			!non-radiative decay rate 3 ->2
	read(25,'(a80)') nothing
	read(25,*)       dwa1(i)
	read(25,'(a80)') nothing
	read(25,*)       dwa2(i)
	read(25,'(a80)') nothing
	read(25,*)       Ndensity(i)		!Numer of atoms per volume in cgs
	Ndensity(i)=Ndensity(i)*atom_gridsize**2
	read(25,'(a80)') nothing
	read(25,*) 		 ne(i)	   			!number of electrons
	read(25,'(a80)') nothing
	read(25,*)  	 n0_0(i) 			!initial value of n0
	read(25,'(a80)') nothing
	read(25,*)  	 n1_0(i) 			!initial value of n1
	read(25,'(a80)') nothing
	read(25,*)  	 n2_0(i) 			!initial value of n2
	read(25,'(a80)') nothing
	read(25,*)  	 n3_0(i) 			!initial value of n3
	read(25,'(a80)') nothing
	read(25,*)  	 current_pumping(i) 
	read(25,'(a80)') nothing
9876 continue
	read(25,'(a80)') nothing
	!--------End of gain parameters--------------
	
	!Spectrum analysis for input pulse 
	read(25,'(a80)') nothing
	read(25,*)       lamdastart
	read(25,'(a80)') nothing
	read(25,*)       lamdaend
	read(25,'(a80)') nothing
	read(25,*)       ifmax			!number of frequency points
	read(25,'(a80)') nothing
	read(25,*)		 tdelay0
	read(25,'(a80)') nothing
	read(25,'(a80)') nothing

	!Spectrum analysis for the output
	read(25,'(a80)') nothing
	read(25,*)       spectrumanalysis  ! 1:yes, 2: no
	read(25,'(a80)') nothing
	read(25,'(a80)') nothing

	!Snapshot parameters
	read(25,'(a80)') nothing
	read(25,*)       snapshotinterval	  ! snapshot interval
	read(25,'(a80)') nothing
	read(25,*)       nsnapshot  ! snapshot interval
	read(25,'(a80)') nothing
	read(25,*)       reso_snapshot_x	  ! snapshot resolution 
	read(25,'(a80)') nothing
	read(25,*)       reso_snapshot_y
	read(25,'(a80)') nothing
	read(25,*)       outputchoice	  ! output choice of snapshot (1:Ez, 2:I, 3:Log(I))
	read(25,'(a80)') nothing
	read(25,*)		period
	read(25,'(a80)') nothing

	! IO option
	read(25,'(a80)') nothing
	read(25,'(a80)') nothing
	read(25,*)  read_bpm		! read waveguide mode from file ( 1: yes, 2: no )
	read(25,'(a80)') nothing
	read(25,*)  read_mpi_seg	  ! read mpich grid segments from file ( 1: yes, 2: no )
	read(25,'(a80)') nothing
	read(25,*)  diag		!(diagnostic output, 1: yes, 0:no)
	read(25,'(a80)') nothing
	read(25,*)  fluxmap		!(fluxmap output, 1: yes, 0:no)
	read(25,'(a80)') nothing
	read(25,*)  t_fluxmap		!(fluxmap output, 1: yes, 0:no)
	read(25,'(a80)') nothing
	read(25,*)  tempdata		!(E(t) output, 1: yes, 0:no)
	read(25,'(a80)') nothing
	read(25,*)  inputwg		!(use specified input waveguide position, 1: yes, 0:no)
	read(25,'(a80)') nothing
	read(25,*)  j_inputwg			!position of input waveguide only works if inputwg=1
	close(25)

! end read user input

	det_period=lamda_f(1,1)/cc*period
!-------------------------------------------------------------
!		derived parameters
!-------------------------------------------------------------
	!FDTD
	dt=dx/cc/2.
	dt_over_dx=dt/dx
	nmax=int(tmax/dt)
!DEC$ IF DEFINED (USE_MPI_CH) 
	if(pid .eq. master)then
		write(*,*) "VERSION:10/28/2004"
		write(*,*) "NMAX: ",nmax
	endif
!DEC$ ELSE ! (USE_MPI_CH) 
		write(*,*) "VERSION:10/28/2004"
		write(*,*) "NMAX: ",nmax
!DEC$ ENDIF ! (USE_MPI_CH) 

	!Field
	TM=0	!number of TM-fields
	TE=0	!number of TE-fields
	Nwgmax = 3
	Nfieldmax=4	!maximum number of field inputs in each waveguide
	
	temp1=tpi*cc
	do 3321 m=1,Nwgmax
		do 3322 n=1,Nfieldmax
			I_input(m,n)=I_input_cgs(m,n)*1.e4
			w_f(m,n)=temp1/lamda_f(m,n)
			k_f(m,n)=tpi/lamda_f(m,n)
3322	continue
		
		do 3323 n=1,nfield(m)			!count the polarization of input fields
			If((I_input(m,n) .ne. 0.) .and. (pol(m,n) .eq. 1))then
				TM=TM+1
			endif
			If((I_input(m,n) .ne. 0.) .and. (pol(m,n) .eq. 0))then
				TE=TE+1
			endif
3323	continue
3321 continue 

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
	write(*,*) "TM: ",TM
	write(*,*) "TE: ",TE
	endif
!DEC$ ELSE
	write(*,*) "TM: ",TM
	write(*,*) "TE: ",TE
!DEC$ ENDIF ! (USE_MPI_CH) 

	!Medium
	do 8765	m=1,atomgroup_number
	Wa1(m)=tpi*cc/lamda_a1(m)
	Wa2(m)=tpi*cc/lamda_a2(m)
	tau3(m) = tau30(m)*tau32(m)/(tau30(m)+tau32(m))
    dt_over_tau30(m)=dt/tau30(m)
	dt_over_tau21(m)=dt/tau21(m)
	dt_over_tau10(m)=dt/tau10(m)
	dt_over_tau32(m)=dt/tau32(m)
	dt_over_tau3(m)=dt/tau3(m)

	kapa1(m)=1./tau21(m)*12.*pi*epsz*cc**3/wa1(m)/wa1(m)*dt*dt/(2.+dWa1(m)*dt)
	rabico1(m)=1./tau21(m)*24.*pi*cc/muz*h_bar_inv/wa1(m)*dt*dt/(2.+dWa1(m)*dt)
	kapa2(m)=1./tau30(m)*12.*pi*epsz*cc**3/wa2(m)/wa2(m)*dt*dt/(2.+dWa2(m)*dt)
	rabico2(m)=1./tau30(m)*12.*pi*cc/muz*h_bar_inv/wa2(m)*dt*dt/(2.+dWa2(m)*dt)
8765 continue

	!snapshot
	nsnapshotinterval = int(snapshotinterval/dt)
	if((tmax/snapshotinterval .ge. 99.).and.(nsnapshot .ge. 99.)) then
		write(*,*) '# of snapshot exceeds 99.'
		pause
		stop
	endif	
						  
!------------------------------------------------------------
! 
! Read Geometry from portable grayscale bitmap
!
!------------------------------------------------------------

	!This is how "graphic4.pgm" looks like----------------------------------------
	!
	!P2
	!# Created by Paint Shop Pro 6
	!500 500
	!255
	!255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 
	!255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 
	!255 255 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
	!0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
	!0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 255 255 
	!255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 
	!........
	!---------------------------------------------------------------

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
	write(*,*) '$475 reading *.pgm file...'   
	endif
!DEC$ ELSE
	write(*,*) '$475 reading *.pgm file...'   
!DEC$ ENDIF ! (USE_MPI_CH) 

	! count number of lines in geometry file
!	linecount=0

!	open(1,file= (trim(rootdir)//'graphic4.pgm'), action='read', share='denywr')
!9997	read (1,*,end=9998) test
!	linecount=linecount+1
!	goto 9997 
!9998 continue
!	close(1)
	   
	! read isize, jsize
	open(1,file= (trim(rootdir)//'graphic4.pgm'), action='read', share='denywr')
	read (1,'(a80)') nothing   !reading header (P2)
    read (1,'(a80)') nothing   !reading header (Created by....)
    read (1,*)  isize,jsize

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
    write(*,*) 'isize= ',isize,', jsize= ', jsize
	write(*,*) ' '
	endif
!DEC$ ELSE
    write(*,*) 'isize= ',isize,', jsize= ', jsize
	write(*,*) ' '
!DEC$ ENDIF ! (USE_MPI_CH) 

	!---------------------------------------------------------------------------------
	! Quick-checking 'observerpoint.csv' to see if they are in range.
	!---------------------------------------------------------------------------------
		open(34,file= (trim(rootdir)//'observerpoint.csv'), action='read', share='denywr')
		read(34,'(a80)') nothing
		read(34,*) nobserver

		do 2340 i=1,nobserver
			read(34,'(a80)') nothing
			read(34,*) inumber
			if(inumber.gt.isize) then
				print *, "$121, observer point out of range"
				call exit (121)
			endif
			read(34,*) jnumber 
			if(jnumber.gt.jsize) then
				print *, "$131, observer point out of range"
				call exit (131)
			endif
2340	 continue
		close(34)
	!-----END OF checking observerpoint------------------

  	allocate (nei(isize,jsize))
	read (1,'(a80)') nothing   !reading the maximum gray scale
	read (1, *) ((nei(i,j),i=1,isize),j=1,jsize) !read color-code and store at nei(i,j)
	close (1)
! adjust for mpi access, if necessary to use with mpich
!DEC$ IF DEFINED (DBG_FILE_WRITE) 
	!--- output portable grayscale file to geonei.csv'

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
!DEC$ ENDIF
	if(diag .eq. 1) then
		open(5, file= (trim(datadir)//'color-code.csv'))
		do 300  j=1,jsize
		do 300  i=1,isize
   			write(5,*) nei(i,j),',',i,',',j
300		continue
		close(5)
	endif
!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF

!DEC$ ENDIF			  
!--------------------------------------------------------
!
!		 allocate the grid matrix
!
!--------------------------------------------------------

	!center grid size (i: x-direction, j: y-direction

	ie = isize ! x size of full grid 
	je = jsize ! y size of full grid
	ib=ie+1
	jb=je+1
	!boundary grid size	(iebc, jebc already defined)
	ibbc=iebc+1
	jbbc=jebc+1
	!center+boundary grid size x - direction
	iefbc=ie+iebc*2
	ibfbc=ib+iebc*2
	!center+boundary grid size y - direction
	jefbc = je + 2*jebc

!DEC$ IF DEFINED (USE_MPI_CH) 

! gs ->
if( read_mpi_seg .eq. 1 ) then	! read grid segmentation from file
!
	call read_segm_mpi( rootdir, je_np, numprocs )
!
	if( je_np(numprocs) .ne. je ) then
		write(*,'(a58)') 'Error 633: incorrect segment info. in seg_info.csv ..'
		call exit
	endif
!
else	! calculate grid segmentation

! ie, je - full grid size; full grid consists from
! the set of all processes subgrids; segmentation (division)
! of full grid is done along x - axis, so subgrids are defined
! by j index array je_np (see below); 
! when allocating size of segments, we consider that number of (arithmetic) operation
! performed on UPML grid is at least twice more than number of operations performed
! on center grid according to UPML equations; thus total number of pixels in vertical
! direction can be considered jt = je + 4*jebc, subgrid size within main grid is 
! dje_np, bottom of master ('0') process grid is at je_np(1); 

! allocate index boundaries for mpich 
	je_np(:) = 0	! je_np(0) = 0, '0' process's subgrid starts from je_np(0)+1 = 1
	je_np(1) = floor((je/real(numprocs)) - (1 - (2/real(numprocs)))*2*jebc)
	dje_np = floor((je + 4*jebc)/real(numprocs)) ! 

	if( numprocs .gt. 2 ) then
	do 304 pid_g = 2, numprocs-1
		je_np(pid_g) = je_np(1) + ( pid_g - 1 )*dje_np
304	continue
	endif

	if( je_np(numprocs-1) .gt. je ) then
		print *, " $590, invalid segmentation configuration, too many processes "
		print *, " should be je > 2*jebc*(N-2) = ", (2*jebc*(numprocs-2)), ", N, &
			(1 < N < 16) is a number of processes !"
		print *, " terminating program ..."
		call exit( 590 )
	endif	

	je_np(numprocs) = je

endif
! <- gs

	if( pid .eq. master ) then !
		308 format(a3,i4,\)
		write(*,'(a12,\)')	' segments:  '
		write(*,308) ((',',je_np(u)),u=0,(numprocs-1))
		write(*,'(a1)') '' ! eof
	endif
 
! (ie, je) grid segmentation is done as follows: 
! process '0' (master): je_np(0) = 1 < j < je_np(1)
! ..
! process 'np' : je_np(np)+1 < j < je_np(np+1)
! ..
! process 'numprocs - 1' : je_np(numprocs - 1)+1 < j < je_np(numprocs) = je
!!
!DEC$ ENDIF ! (USE_MPI_CH) 


	! index
	allocate(epsrf(iefbc,jebc),epsrb(iefbc,jebc), &
  	    	 epsrbcl(iebc,je), epsrbcr(iebc,je))

	!TM mode
	!if(TM .gt. 0)then
!DEC$ IF DEFINED (USE_MPI_CH) 
		j1c = je_np(pid)+1
		j2c = je_np(pid+1)
		if( pid .eq. 0 ) then
			j1_al=j1c
			j2_al=j2c+1
		elseif( pid .eq. (numprocs-1) ) then
			j1_al=j1c-1
			j2_al=j2c
		else
			j1_al=j1c-1
			j2_al=j2c+1
		endif
		allocate(Ez(ie,j1_al:j2_al),Hx(ie,j1_al:j2_al+1),Hy(ib,j1_al:j2_al))

		! medium
		allocate(Ez_old(ie,j1_al:j2_al))
		allocate(Avptz(ie,j1_al:j2_al),Avptz_old(ie,j1_al:j2_al))
		allocate(pz_1(ie,j1_al:j2_al,atomgroup_number),pz_old1_1(ie,j1_al:j2_al,atomgroup_number),pz_2(ie,j1_al:j2_al,atomgroup_number), &
	pz_old1_2(ie,j1_al:j2_al,atomgroup_number))
		allocate(n0(ie,j1_al:j2_al,atomgroup_number),n1(ie,j1_al:j2_al,atomgroup_number),n2(ie,j1_al:j2_al,atomgroup_number),n3(ie,j1_al:j2_al,atomgroup_number),  &
			 n0_old(ie,j1_al:j2_al,atomgroup_number),n1_old(ie,j1_al:j2_al,atomgroup_number),n2_old(ie,j1_al:j2_al,atomgroup_number),n3_old(ie,j1_al:j2_al,atomgroup_number))
		allocate(n0_result(ie,j1_al:j2_al),n1_result(ie,j1_al:j2_al),n2_result(ie,j1_al:j2_al),n3_result(ie,j1_al:j2_al))
		! end of medium

		allocate(Ezf(iefbc,jebc),Ezb(iefbc,jebc),  &
  			   	 Ezbcl(iebc,j1_al:j2_al), Ezbcr(iebc,j1_al:j2_al))
		allocate(Hxf(iefbc,jbbc),Hxb(iefbc,jbbc), &
				Hyf(ibfbc,jebc),Hyb(ibfbc,jebc), &
				Hxbcl(iebc,j1_al:j2_al+1),Hxbcr(iebc,j1_al:j2_al+1), &
				Hybcl(ibbc,j1_al:j2_al),Hybcr(ibbc,j1_al:j2_al))
		allocate(DEzf(iefbc,jebc),DEzb(iefbc,jebc), &
  			     DEzbcl(iebc,j1_al:j2_al), DEzbcr(iebc,j1_al:j2_al))
		allocate(BHxf(iefbc,jbbc),BHxb(iefbc,jbbc), &
				BHyf(ibfbc,jebc),BHyb(ibfbc,jebc), &
			    BHxbcl(iebc,j1_al:j2_al+1),BHxbcr(iebc,j1_al:j2_al+1), &
				BHybcl(ibbc,j1_al:j2_al),BHybcr(ibbc,j1_al:j2_al))
!DEC$ ELSE
		allocate(Ez(ie,je),Hx(ie,jb),Hy(ib,je))
	
		! medium
		allocate(Ez_old(ie,je))
		allocate(Avptz(ie,je),Avptz_old(ie,je))
		allocate(pz_1(ie,je,atomgroup_number),pz_old1_1(ie,je,atomgroup_number),pz_2(ie,je,atomgroup_number), &
			pz_old1_2(ie,je,atomgroup_number))
		allocate(n0(ie,je,atomgroup_number),n1(ie,je,atomgroup_number),n2(ie,je,atomgroup_number),n3(ie,je,atomgroup_number),  &
			 n0_old(ie,je,atomgroup_number),n1_old(ie,je,atomgroup_number),n2_old(ie,je,atomgroup_number),n3_old(ie,je,atomgroup_number))
		allocate(n0_result(ie,je),n1_result(ie,je),n2_result(ie,je),n3_result(ie,je))
		! end of medium
	
		allocate(Ezf(iefbc,jebc),Ezb(iefbc,jebc),  &
  				 Ezbcl(iebc,je), Ezbcr(iebc,je))
		allocate(Hxf(iefbc,jbbc),Hxb(iefbc,jbbc), &
	             Hyf(ibfbc,jebc),Hyb(ibfbc,jebc), &
		         Hxbcl(iebc,jb),Hxbcr(iebc,jb), &
				 Hybcl(ibbc,je),Hybcr(ibbc,je))
		allocate(DEzf(iefbc,jebc),DEzb(iefbc,jebc), &
  				DEzbcl(iebc,je), DEzbcr(iebc,je))
		allocate(BHxf(iefbc,jbbc),BHxb(iefbc,jbbc), &
				BHyf(ibfbc,jebc),BHyb(ibfbc,jebc), &
				BHxbcl(iebc,jb),BHxbcr(iebc,jb), &
				BHybcl(ibbc,je),BHybcr(ibbc,je))
!DEC$ ENDIF ! (USE_MPI_CH)
	
			
	!endif
	
	!TE mode
	if(TE .gt. 0)then
!DEC$ IF DEFINED (USE_MPI_CH) 
		j1c = je_np(pid)+1
		j2c = je_np(pid+1)
		if( pid .eq. 0 ) then
			j1_al=j1c
			j2_al=j2c+1
		elseif( pid .eq. (numprocs-1) ) then
			j1_al=j1c-1
			j2_al=j2c
		else
			j1_al=j1c-1
			j2_al=j2c+1
		endif
		allocate(Hz(ie,j1_al:j2_al),Ex(ie,j1_al:j2_al+1),Ey(ib,j1_al:j2_al))
		
		! medium
		allocate(Ex_old(ie,j1_al:j2_al+1),Ey_old(ib,j1_al:j2_al))
		allocate(Avptx(ie,j1_al:j2_al+1),Avptx_old(ie,j1_al:j2_al+1))
		allocate(Avpty(ib,j1_al:j2_al),Avpty_old(ib,j1_al:j2_al))
 		allocate(px_1(ie,j1_al:j2_al+1,atomgroup_number),px_old1_1(ie,j1_al:j2_al+1,atomgroup_number),px_2(ie,j1_al:j2_al+1,atomgroup_number),px_old1_2(ie,j1_al:j2_al+1,atomgroup_number))
		allocate(py_1(ib,j1_al:j2_al,atomgroup_number),py_old1_1(ib,j1_al:j2_al,atomgroup_number),py_2(ib,j1_al:j2_al,atomgroup_number),py_old1_2(ib,j1_al:j2_al,atomgroup_number))
		allocate(n0(ie,j1_al:j2_al,atomgroup_number),n1(ie,j1_al:j2_al,atomgroup_number),n2(ie,j1_al:j2_al,atomgroup_number),n3(ie,j1_al:j2_al,atomgroup_number),  &
			 n0_old(ie,j1_al:j2_al,atomgroup_number),n1_old(ie,j1_al:j2_al,atomgroup_number),n2_old(ie,j1_al:j2_al,atomgroup_number),n3_old(ie,j1_al:j2_al,atomgroup_number))
		allocate(n0_result(ie,j1_al:j2_al),n1_result(ie,j1_al:j2_al),n2_result(ie,j1_al:j2_al),n3_result(ie,j1_al:j2_al))

		! end medium

		allocate(Hzf(iefbc,jebc),Hzb(iefbc,jebc),  &
  				Hzbcl(iebc,j1_al:j2_al), Hzbcr(iebc,j1_al:j2_al))
		allocate(Exf(iefbc,jbbc),Exb(iefbc,jbbc), &
				Eyf(ibfbc,jebc),Eyb(ibfbc,jebc), &
				Exbcl(iebc,j1_al:j2_al+1),Exbcr(iebc,j1_al:j2_al+1), &
				Eybcl(ibbc,j1_al:j2_al),Eybcr(ibbc,j1_al:j2_al))
		allocate(BHzf(iefbc,jebc),BHzb(iefbc,jebc), &
				BHzbcl(iebc,j1_al:j2_al), BHzbcr(iebc,j1_al:j2_al))
		allocate(DExf(iefbc,jbbc),DExb(iefbc,jbbc), &
				DEyf(ibfbc,jebc),DEyb(ibfbc,jebc), &
				DExbcl(iebc,j1_al:j2_al+1),DExbcr(iebc,j1_al:j2_al+1), &
				DEybcl(ibbc,j1_al:j2_al),DEybcr(ibbc,j1_al:j2_al))
!DEC$ ELSE
		allocate(Hz(ie,je),Ex(ie,jb),Ey(ib,je))
		
		! medium
		allocate(Ex_old(ie,jb),Ey_old(ib,je))
		allocate(Avptx(ie,jb),Avptx_old(ie,jb))
		allocate(Avpty(ib,je),Avpty_old(ib,je))
		allocate(px_1(ie,jb,atomgroup_number),px_old1_1(ie,jb,atomgroup_number),px_2(ie,jb,atomgroup_number),px_old1_2(ie,jb,atomgroup_number))
		allocate(py_1(ib,je,atomgroup_number),py_old1_1(ib,je,atomgroup_number),py_2(ib,je,atomgroup_number),py_old1_2(ib,je,atomgroup_number))
		allocate(n0(ie,je,atomgroup_number),n1(ie,je,atomgroup_number),n2(ie,je,atomgroup_number),n3(ie,je,atomgroup_number),  &
			 n0_old(ie,je,atomgroup_number),n1_old(ie,je,atomgroup_number),n2_old(ie,je,atomgroup_number),n3_old(ie,je,atomgroup_number))
		allocate(n0_result(ie,je),n1_result(ie,je),n2_result(ie,je),n3_result(ie,je))
		! end medium

		allocate(Hzf(iefbc,jebc),Hzb(iefbc,jebc),  &
	  	    	 Hzbcl(iebc,je), Hzbcr(iebc,je))
		allocate(Exf(iefbc,jbbc),Exb(iefbc,jbbc), &
				Eyf(ibfbc,jebc),Eyb(ibfbc,jebc), &
		         Exbcl(iebc,jb),Exbcr(iebc,jb), &
				 Eybcl(ibbc,je),Eybcr(ibbc,je))
		allocate(BHzf(iefbc,jebc),BHzb(iefbc,jebc), &
	  	         BHzbcl(iebc,je), BHzbcr(iebc,je))
		allocate(DExf(iefbc,jbbc),DExb(iefbc,jbbc), &
		         DEyf(ibfbc,jebc),DEyb(ibfbc,jebc), &
		         DExbcl(iebc,jb),DExbcr(iebc,jb), &
				 DEybcl(ibbc,je),DEybcr(ibbc,je))
!DEC$ ENDIF ! (USE_MPI_CH)
	! allocate(Ex_old(ie,jb),Ey_old(ib,je))
	
	endif	!TEmode

																	  
	!UPML coefficients-----------------------
	allocate(axef(iefbc),ayef(jebc),axhf(ibfbc),ayhf(jbbc), &
 			 axeb(iefbc),ayeb(jebc),axhb(ibfbc),ayhb(jbbc), &
			 axel(iebc) ,ayel(je)  ,axhl(ibbc) ,ayhl(jb), &
			 axer(iebc) ,ayer(je)  ,axhr(ibbc) ,ayhr(jb), &
	         bxef(iefbc),byef(jebc),bxhf(ibfbc),byhf(jbbc), &
			 bxeb(iefbc),byeb(jebc),bxhb(ibfbc),byhb(jbbc), &
			 bxel(iebc) ,byel(je)  ,bxhl(ibbc) ,byhl(jb), &
			 bxer(iebc) ,byer(je)  ,bxhr(ibbc) ,byhr(jb), &
	         cxef(iefbc),cyef(jebc),cxhf(ibfbc),cyhf(jbbc), &
			 cxeb(iefbc),cyeb(jebc),cxhb(ibfbc),cyhb(jbbc), &
			 cxel(iebc) ,cyel(je)  ,cxhl(ibbc) ,cyhl(jb), &
			 cxer(iebc) ,cyer(je)  ,cxhr(ibbc) ,cyhr(jb))

    !index grid & gain grid--------------------
	!!
	allocate(geometry(ie,je),detector(ie,je))
	
	!! flux detector, map arrays


!------------------------------------------------------------------- 
!
!				GEOMETRY SPECIFICATION.
!
!-------------------------------------------------------------------

	!---- read index profile
	open(15,file= (trim(rootdir)//'indexProfile.csv'), action='read', share='denywr')

	!-------This is how "indexProfile.csv" looks like----------
	!
	!---index number----- 
	!3
	!--- first color index (black: waveguide without gain)
	!0
	!----refraction  index
	!3.
	!----gain
	!0.
	!----second color index (white: air)
	!255
	!----refraction index
	!1.
	!----gain
	!0.
	!----third color index (gray: waveguide with gain)
	!125
	!----refraction index
	!3.
	!----gain
	!1.
	!
	!--------------------------------------------------------

	read (15,'(a80)') nothing	!reading the header	"---index number-------"
	read (15,*) num_index		!number of color codes

	allocate(i_color(num_index),n_index(num_index),ni_index(num_index),sigma_index(num_index),gain_index(num_index))
	allocate(Mur(num_index),Eps(num_index))
	!!
	allocate(det_index(num_index),det_dir(num_index),det_xp(num_index),det_yp(num_index), &
		det_xlen(num_index),det_ylen(num_index),spec_index(num_index),n_skip(num_index))
	det_num=0
	ftdet_num=0
	metal=0
	gain=0
	sigma_index(:)=0 !gy  added line

	do 270 i=1,num_index
		read (15,'(a80)') nothing
		read(15,*)  i_color(i)
		read (15,'(a80)') nothing
		read(15,*)  n_index(i)
		read (15,'(a80)') nothing
		read(15,*)  ni_index(i)
		read (15,'(a80)') nothing
		read(15,*)  gain_index(i)
		read (15,'(a80)') nothing
		read(15,*)	det_index(i)
		read (15,'(a80)') nothing
		read(15,*)	spec_index(i)
		read (15,'(a80)') nothing
		read(15,*)	n_skip(i)
		if (det_index(i).ge.1) then
			det_num=det_num+1
			det_dir(det_num)=det_index(i)
			det_index(i)=det_num
			if(spec_index(i) .eq. 1)then
				ftdet_num=ftdet_num+1
				spec_index(det_num)=spec_index(i)
				n_skip(det_num)=n_skip(i)
			endif
		endif
		if (ni_index(i).ne.0.) then
			metal=metal+1
			!sigma_index(i)=2*w_f(1,1)*n_index(i)*ni_index(i) 10/04/05
			sigma_index(i)=2*w_f(1,1)*n_index(i)*ni_index(i)*epsz
		endif
		if (gain_index(i).ne.0.) then
			gain=gain+1
		endif

!--------allocate detector results--------------------------
		allocate(det_result(det_num))
!DEC$ IF DEFINED (USE_MPI_CH) 
		allocate(det_result_t(det_num))
!DEC$ ENDIF ! (USE_MPI_CH) 
!-----------------------------------------------------------

270	continue
	close(15)

	if(metal .gt. 0)then
		allocate(sigmageo(ie,je))
	endif
	if(gain .gt. 0)then
		allocate(gaingeo(ie,je))
	endif
 	!!
!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
	write(*,*) 'total detectors:', det_num
	write(*,*) 'total spectrum detectors:', ftdet_num
	write(*,*) 'Number of metallic materials:',metal
	write(*,*) 'Gain material? (0:NO, 1:YES)',gain
	endif
!DEC$ ELSE
	write(*,*) 'total detectors:', det_num
	write(*,*) 'total spectrum detectors:', ftdet_num
	write(*,*) 'Number of metallic materials:',metal
	write(*,*) 'Gain material? (0:NO, 1:YES)',gain
!DEC$ ENDIF ! (USE_MPI_CH) 

	rnf=n_index(1)
	det_time=0.

	!initialize the center grid with free space for index and detector
	do 2100 j=1,je
	do 2100 i=1,ie
		    geometry(i,j)=1.0
			detector(i,j)=0
2100		continue

	!initialize the center grid with free space for conductivity
	if(metal .gt. 0)then
	do 2200 j=1,je
	do 2200 i=1,ie
		    sigmageo(i,j)=0.0
2200		continue
	endif

	!initialize the center grid with free space for conductivity
	if(gain .gt. 0)then
	do 2300 j=1,je
	do 2300 i=1,ie
		    gaingeo(i,j)=0.0
2300		continue
	endif

	!index & gain assignment for the center grid
	! geometry = 1/n^2
	! sigmageo = conductivity
	! gaingeo >= 1 (gain) or 0 (no gain)  
	det_xp(:)=0
	det_yp(:)=0 
	do 210 j=1,je  	  
   	do 210 i=1,ie
		do 2105 k=1,num_index
			if (nei(i,j) .eq. i_color(k)) then
				geometry(i,j)=1./(n_index(k))**2	!geometry=1/n^2
				detector(i,j)=det_index(k)
				if(metal .gt. 0)then
					geometry(i,j)=1./((n_index(k))**2-(ni_index(k))**2) 
					sigmageo(i,j)=sigma_index(k)
				endif
				if(gain .gt. 0)then
					if((mod(i,atom_gridsize).eq.0).and.(mod(j,atom_gridsize).eq.0))then
						gaingeo(i,j)=gain_index(k)
					endif
				endif
				if((det_num .ge. 1).AND.(det_index(k).ge.1)) then  !gy modified, failed when no detectors
					m=det_index(k)
					if((det_xp(m) .eq. 0))then
						det_xp(m)=i
						det_yp(m)=j
					endif
				endif
			endif	 
2105	continue
210 continue
	
	det_xlen(:)=1
	det_ylen(:)=1
	!find out detector length
	do 212 k=1,det_num
		if((det_xp(k).eq.0).OR.(det_yp(k).eq.0)) then !detector is listed in index file but not in bitmap file
			goto 212
		endif
		if(det_dir(k) .eq. 1)then
			ii=det_xp(k)
			do 213 j=det_yp(k)+1,je
				if(nei(ii,j) .eq. nei(ii,j-1)) then
					det_ylen(k)= det_ylen(k)+1
				else
					goto 212
				endif
213			continue
		else
			jj=det_yp(k)
			do 214 i=det_xp(k)+1,ie
				if(nei(i,jj) .eq. nei(i-1,jj)) then
					det_xlen(k)= det_xlen(k)+1
				else
					goto 212
				endif
214			continue
		endif
		
212	continue

	det_maxlen=1
	do 216 k=1,det_num
		if(INT(det_xlen(k)/n_skip(k)) .ge. det_maxlen)then
		det_maxlen= INT(det_xlen(k)/n_skip(k))
		endif
		if(INT(det_ylen(k)/n_skip(k)) .ge. det_maxlen)then
		det_maxlen= INT(det_ylen(k)/n_skip(k))
		endif
216	continue

	deallocate(nei)	

	if((ftdet_num .ge. 1).and.(TM .ge. 1))then
		allocate(ftrd_z(0:ifmax,0:det_maxlen,det_num),ftid_z(0:ifmax,0:det_maxlen,det_num), &
			fte2d_z(0:ifmax,0:det_maxlen,det_num),fte2d_zt(0:ifmax,0:det_maxlen,det_num))
		allocate(fte2maxd_z(det_num))
	endif
	if((ftdet_num .ge. 1).and.(TE .ge. 1))then
		allocate(ftrd_x(0:ifmax,0:det_maxlen,det_num),ftid_x(0:ifmax,0:det_maxlen,det_num), &
			fte2d_x(0:ifmax,0:det_maxlen,det_num),fte2d_xt(0:ifmax,0:det_maxlen,det_num), &
			ftrd_y(0:ifmax,0:det_maxlen,det_num),ftid_y(0:ifmax,0:det_maxlen,det_num),  &
			fte2d_y(0:ifmax,0:det_maxlen,det_num),fte2d_yt(0:ifmax,0:det_maxlen,det_num))
		allocate(fte2maxd_x(det_num),fte2maxd_y(det_num))
	endif
	

	! -------- output geometry file -----------------
!DEC$ IF DEFINED (DBG_FILE_WRITE) 

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
!DEC$ ENDIF ! (USE_MPI_CH) 

	if(diag .eq. 1)then
		open(2102,file= (trim(datadir)//'geometry.csv'))
		open(2104,file= (trim(datadir)//'det_profile.csv'))
		if(metal .gt. 0)then
		open(2106,file= (trim(datadir)//'sigmageo.csv'))
		endif
		if(gain .gt. 0)then
		open(2103,file= (trim(datadir)//'gaingeo.csv'))
		endif
		do 2101 j=1,je
		do 2101 i=1,ie
			!write(2102,*) geometry(i,j) 10/4/05
			caij=(1-sigmageo(i,j)*dt/2./epsz*geometry(i,j))/(1+sigmageo(i,j)*dt/2./epsz*geometry(i,j))
			write(2102,*) caij
			write(2104,*) detector(i,j)
			if(metal.gt.0)then
			write(2106,*) sigmageo(i,j)
			endif
			if(gain.gt.0)then
			write(2103,*) gaingeo(i,j)
			endif
2101	continue
		close(2102)
		close(2104)
		if(metal.gt.0)then
		close(2106)
		endif
		if(gain.gt.0)then
		close(2103)
		endif
	endif

!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF ! (USE_MPI_CH) 

!DEC$ ENDIF			

	! ---------Optical-Intensity to Electric-field conversion
	temp1=2./(cc*epsz*rnf)
	do 979 m=1,Nwgmax
		do 980 n=1,Nfieldmax
		E_input(m,n)=sqrt(I_input(m,n)*temp1) 
980		continue
979	continue


	!----------UPML index assignment

	!initialize the UPML index "epsr" = epsz
		! front & back BC
		do 2050 j=1,jebc
		do 2050 i=1,iefbc
			epsrf(i,j)=epsz
			epsrb(i,j)=epsz
2050	continue

		! left & right BC
		do 2051 j=1,je
		do 2051 i=1,iebc
			epsrbcl(i,j)=epsz
			epsrbcr(i,j)=epsz
2051	continue

	!assign the UPML index "epsr" = epsz/geometry = epsz*n^2
		! front & back BC
		do 2052 j=1,jebc
		do 2052 i=1,ie
			epsrf(i+iebc,j)=epsz/geometry(i,1)
			epsrb(i+iebc,j)=epsz/geometry(i,je)
2052	continue

		do 2054 j=1,je
		do 2054 i=1,iebc
			epsrbcl(i,j)=epsz/geometry(1,j)
			epsrbcr(i,j)=epsz/geometry(ie,j)
2054	continue

!DEC$ IF DEFINED (DBG_FILE_WRITE) 

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
!DEC$ ENDIF ! (USE_MPI_CH) 

		if(diag .eq. 1) then
			open(222,file= (trim(datadir)//'epsrbcl.csv'))
			open(223,file= (trim(datadir)//'epsrbcr.csv'))
			! left & right BC
			do 2055 j=1,je
			do 2055 i=1,iebc
				write(222,*) epsrbcl(i,j)
				write(223,*) epsrbcr(i,j)
2055		continue
			close(222)
			close(223)
		endif

!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF ! (USE_MPI_CH) 

!DEC$ ENDIF			
	
!----------------------------------------------------------------------
!
!		waveguide spec & userinputs
!
!----------------------------------------------------------------------

	!-------------------------------------------------------------------------- 
	!	Detection of waveguide geometry from the bitmap input file
	!	
	!	1. # of waveguides: num_wg = num_wgx (in x-direction) + num_wgy (in y-direction)
	!	2. wgstart() : starting position of waveguides
	!	3. wgend() : ending position of waveguides
	!	4. width() : waveguide width
	!	5. wgindex() : waveguide index
	!	6. wgdirection() : waveguide direction (1:x,2:y)
	!
	!	Output the waveguide position:'x_wgposition.csv' & "y_wgposition.csv'
	!	generates two columns in the format of {j, geometry(1,j)}
	!	geometry number is not one when for waveguide region
	!---------------------------------------------------------------------------

	if ( read_bpm .eq. 2 ) then
		
		if(inputwg .eq. 1)then !if input wageguide position is specified
		
			num_wg = 1
			wgdirection(1)=1
			do 2107 j=j_inputwg,1,-1	
				if(geometry(1,j) .gt. geometry(1,j+1))then
					wgstart(num_wg)=j+1
					if((geometry(1,j_inputwg).lt.0).or.(geometry(1,j).lt.0)) then
						print *, "$972, real index is negative at the boundary."
						print *, "Do not put metal at the boundary"
						call exit (972)
					endif
					wgindex(num_wg)=sqrt(1./geometry(1,j_inputwg))
					cladindex(num_wg)=sqrt(1./geometry(1,j))
					goto 21077
				endif
2107		continue
21077		continue
			do 2108 j=j_inputwg,je	
				if(geometry(1,j) .gt. geometry(1,j-1))then
					wgend(num_wg)=j-1
					width(num_wg)=wgend(num_wg)-wgstart(num_wg)+1
					goto 21088
				endif
2108		continue
21088		continue		
		else !input wave position condition
		
			num_wg = 0
			wg_finder = 0
			do 2109 j=1,je
				if(j .ge. 2)then
					if(geometry(1,j) .lt. geometry(1,j-1))then
						wg_finder = 1
						if (num_wg.ge.4) then
							num_wg = num_wg
						else
							num_wg = num_wg + 1
						endif
						wgdirection(num_wg)=1
						wgstart(num_wg)=j
						if((geometry(1,j).lt.0).or.(geometry(1,j-1).lt.0)) then
							print *, "$973, real index is negative at the boundary."
							print *, "Do not put metal at the boundary"
							call exit (973)
						endif
						wgindex(num_wg)=sqrt(1./geometry(1,j))
						cladindex(num_wg)=sqrt(1./geometry(1,j-1))
					elseif((geometry(1,j) .gt. geometry(1,j-1)).and.(wg_finder.eq.1))then
						wgend(num_wg)=j-1
						width(num_wg)=wgend(num_wg)-wgstart(num_wg)+1
					endif
				endif
2109		continue
		endif !end of input wg position condition

	else	! read mode from BPM file: 1 waveguide in x direction; 
		num_wg = 1
		wgdirection(1)=1
	endif

	num_wgx=num_wg
	wg_finder = 0

	do 2110 i=1,ie
		if(i .ge. 2)then
			if(geometry(i,1) .lt. geometry(i-1,1))then
				if (num_wg.ge.4) then
					 num_wg = num_wg
				else
					num_wg = num_wg + 1
				endif
				wg_finder = 1
				wgdirection(num_wg)=2
				wgstart(num_wg)=i
				if((geometry(i,1).lt.0).or.(geometry(i-1,1).lt.0)) then
						print *, "$974, real index is negative at the boundary."
						print *, "Do not put metal at the boundary"
						call exit (974)
				endif
				wgindex(num_wg)=sqrt(1./geometry(i,1))
				cladindex(num_wg)=sqrt(1./geometry(i-1,1))
			elseif((geometry(i,1) .gt. geometry(i-1,1)).and.(wg_finder.eq.1))then
				wgend(num_wg)=i-1
				width(num_wg)=wgend(num_wg)-wgstart(num_wg)+1
			endif
2110 endif

!DEC$ IF DEFINED (DBG_FILE_WRITE) 

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
!DEC$ ENDIF ! (USE_MPI_CH) 
	
	if(diag .eq. 1)then
		open(333,file= (trim(datadir)//'x_wgposition.csv'))
		do 2112 j=1,je
			write(333,*) j,',',geometry(1,j)
2112	continue
		close(333)

		open(334,file= (trim(datadir)//'y_wgposition.csv'))
		do 2114 i=1,ie
			write(334,*) i,',',geometry(i,1)
2114	continue
		close(334)
	endif

!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF ! (USE_MPI_CH) 

!DEC$ ENDIF			

	num_wgy=num_wg-num_wgx
	maxgrid=max(ie,je)
	allocate(rmode(num_wg,Nfieldmax,maxgrid))

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
!DEC$ ENDIF ! (USE_MPI_CH) 

	write(*,*) 'total number of waveguide: ',num_wg
	write(*,*) 'number of waveguide in x-direction: ',num_wgx
	write(*,*) 'number of waveguide in y-direction: ',num_wgy
	do 2111 m=1,num_wg
		write(*,*) 'waveguide #',m
		if(wgdirection(m) .eq. 1)then
			if ( read_bpm .eq. 2 ) then
			write(*,*) 'x-waveguide start',wgstart(m)
			write(*,*) 'x-waveguide end',wgend(m)
 			write(*,*) 'x-waveguide width',width(m),'=',width(m)*dx*1.e6,'um'
			write(*,*) 'x-waveguide index',wgindex(m)
			write(*,*) 'x-waveguide cladding index',cladindex(m)
			endif
		else
			write(*,*) 'y-waveguide start',wgstart(m)
			write(*,*) 'y-waveguide end',wgend(m)
 			write(*,*) 'y-waveguide width',width(m),'=',width(m)*dy*1.e6,'um'
			write(*,*) 'y-waveguide index',wgindex(m)
			write(*,*) 'y-waveguide cladding index',cladindex(m)
		endif
2111 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF ! (USE_MPI_CH) 

!--------------------------------------------------------------------------------
!
!			Fundamental waveguide mode solver (assumed symmetric slab waveguide)
!				Both TM and TE
!
!--------------------------------------------------------------------------------
rmode(:,:,:)=0.
if ( read_bpm .eq. 2 ) then
    do 2155 m=1,num_wg
		do 2156 k=1,nfield(m)
			wgwidth=width(m)*dx	! = 2d
 			wgwidth2=wgwidth/2.	 ! = d
			epsnwg=epsz*wgindex(m)*wgindex(m)  !eps*n^2
			epsnclad=epsz*cladindex(m)*cladindex(m)
			a0=(w_f(m,k)*wgwidth2)**2*muz*(epsnwg-epsnclad)	  ! (k0*d)^2*(n^2-1)
			xx=0.0000001  ! initial value for Newton's method
			xx_old=xx
			If(pol(m,k) .eq. 1)then
				fx=tan(xx)-sqrt(a0-xx**2)/xx
			else
				fx=(epsnclad/epsnwg)*tan(xx)-sqrt(a0-xx**2)/xx
			endif
			do 1920 n=1,10000
				if(pol(m,k) .eq. 1)then
					Gx=1/(cos(xx)**2)-(-xx**2/sqrt(a0-xx**2)-sqrt(a0-xx**2))/xx**2
				else
					Gx=(epsnclad/epsnwg)/(cos(xx)**2)-(-xx**2/sqrt(a0-xx**2)-sqrt(a0-xx**2))/xx**2
				endif
				xx_old=xx			
				xx=xx-Fx/Gx*0.1
				if(abs(xx-xx_old) .le. 1.e-7*abs(xx_old))then
				goto 1930
				endif
				If(pol(m,k) .eq. 1)then	 !TM
					fx=tan(xx)-sqrt(a0-xx**2)/xx
				else					 !TE
					fx=(epsnclad/epsnwg)*tan(xx)-sqrt(a0-xx**2)/xx
				endif
1920		continue
1930		continue
			kk=xx/(wgwidth2)   !transversal wave vector in the waveguide (h in Yariv's notation)
			kclad=(a0/wgwidth2**2-kk**2)**0.5   ! exponent outside the waveguide (p in Yariv's notation)
			kprop(m,k)=	sqrt((k_f(m,k)*wgindex(m))**2-kk**2)  !propagation wave vector, (beta in Yariv's notation)
  
			jws1=wgstart(m)
			jwe1=wgend(m)
			cent=(float(jws1)+ float(jwe1-jws1)/2.0)*dx		
			!Cosine in core region
			do 141 j=jws1,jwe1
				x=dx*float(j)-cent
				rmode(m,k,j)=cos(kk*x)
141			continue
			! Exponential decay in lower (front) air region
			xlow=dx*float(jws1)-cent
			vallow=cos(kk*xlow)
			do 143 j=4,jws1-1
				x=dx*float(j-jws1)
				rmode(m,k,j)=exp(kclad*x)*vallow
143			continue		
			!Exponential decay in upper (back) air region
			xhigh=dx*float(jwe1)-cent
			valhigh=cos(kk*xhigh)
			if(wgdirection(m) .eq. 1)then
				jmax = je
			else
				jmax = ie
			endif
			do 145 j=jwe1+1,jmax-4
				x=dx*float(jwe1-j)
				rmode(m,k,j)=exp(kclad*x)*valhigh
145			continue  

2156	continue 
2155 continue

!DEC$ IF DEFINED (DBG_FILE_WRITE) 

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
!DEC$ ENDIF ! (USE_MPI_CH) 

			if(diag .eq. 1) then
				open(10,file= (trim(datadir)//'wg_mode.csv'))   !write converging solution via Newton's method
				do 2157 m=1,num_wg
				do 2157 k=1,nfield(m)
					write(10,*) 'k0',',',k_f(m,k)
					write(10,*) 'kk (h)',',',kk
					write(10,*) 'kair (p)',',',kclad
					write(10,*) 'kprop (beta)',',',kprop(m,k)
2157			continue
				close(10)
			endif

!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF ! (USE_MPI_CH) 

!DEC$ ENDIF			

endif

! get rmode(1,1,j) data from 's_tx_0.m00'; 
! works for a single waveguide, x-direction;  
if ( read_bpm .eq. 1 ) then

	if(pol(1,1) .eq. 1)then	 !TM
	open(16, file= (trim(rootdir)//'s_tm_0.m00'), action='read', share='denywr')
	else					 !TE
	open(16, file= (trim(rootdir)//'s_te_0.m00'), action='read', share='denywr')
	endif

!	read(16,'(a80)') nothing
!	read(16,'(a80)') nothing
!	read(16,'(a80)') nothing ! shift array by 1 

	do 302  j=1,jsize
   		read(16,*) temp1,rmode(1,1,j)
302		continue

endif


!DEC$ IF DEFINED (DBG_FILE_WRITE) 

!DEC$ IF DEFINED (USE_MPI_CH) 
if( pid .eq. master ) then
!DEC$ ENDIF ! (USE_MPI_CH) 
if(diag .eq. 1)then
! create filenames, like: rmodefilename(1,1)='rmode1_1.csv', etc.
	do 144 u = 1,Nwgmax
	do 144 v = 1,Nfieldmax
		rmf_ix1 = 48 + u
		rmf_ix2 = 48 + v
		rmodefilename(u,v) = trim(datadir)//'rmode'//char(rmf_ix1)//'_'//char(rmf_ix2)//'.csv' 
144 continue
!!
	do 147 m=1,num_wg
		do 148 k=1,nfield(m)
			open(1000+m*100+k,file=trim(rmodefilename(m,k)))
			if(wgdirection(m) .eq. 1)then
				jmax = je
			else
				jmax = ie
			endif
			do 146 j=1,jmax             
146				write(1000+m*100+k,*)j,',',rmode(m,k,j)
			close(1000+m*100+k)
148		continue
147 continue
endif

!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF ! (USE_MPI_CH) 

!DEC$ ENDIF			

!------------------------------------------------------------------------
!
!		 filenames for snapshot output
!
!------------------------------------------------------------------------
! afilename(1)='Ez01.csv' .. afilename(99)='Ez99.csv'
! ffilename(1)='Ex01.csv' .. ffilename(99)='Ex99.csv'
! gfilename(1)='Ey01.csv' .. ffilename(99)='Ey99.csv'
!!
	do 149 u = 0,9
	do 149 v = 0,9
		rmf_ix1 = 48 + u
		rmf_ix2 = 48 + v
		w = u*10 + v
		afilename(w) = trim(datadir)//'Ez'//char(rmf_ix1)//char(rmf_ix2)//'.csv' 
!medium
		bfilename(w) = trim(datadir)//'N0'//char(rmf_ix1)//char(rmf_ix2)//'.csv' 
		cfilename(w) = trim(datadir)//'N1'//char(rmf_ix1)//char(rmf_ix2)//'.csv' 
		dfilename(w) = trim(datadir)//'N2'//char(rmf_ix1)//char(rmf_ix2)//'.csv' 
		efilename(w) = trim(datadir)//'N3'//char(rmf_ix1)//char(rmf_ix2)//'.csv' 
!end medium
		ffilename(w) = trim(datadir)//'Hz'//char(rmf_ix1)//char(rmf_ix2)//'.csv' 
		gfilename(w) = trim(datadir)//'Ey'//char(rmf_ix1)//char(rmf_ix2)//'.csv' 
149 continue
!!
! set size, allocate observer arrays
	ns_obs = 512 ! number of time steps saved in observer arrays
	allocate(m_obs(ns_obs))
	allocate(tt_obs(ns_obs))
	m_obs(:) = 0
	tt_obs(:) = 0.
!

!---------------------------------------------------------------------------------
!
!		Specify Observers: read from 'observerpoint.csv'
!
!---------------------------------------------------------------------------------
	open(35,file= (trim(rootdir)//'observerpoint.csv'), action='read', share='denywr')
	read(35,'(a80)') nothing
	read(35,*) nobserver
	allocate(iobs(nobserver),jobs(nobserver))

!DEC$ IF DEFINED (USE_MPI_CH) 
		allocate(obs_in_proc(nobserver))
!DEC$ ENDIF ! (USE_MPI_CH) 

	do 2350 i=1,nobserver
		read(35,'(a80)') nothing
		read(35,*) iobs(i)
		read(35,*) jobs(i) 
2350 continue
	iobs(3)=atom_gridsize*int(iobs(3)/atom_gridsize)
	jobs(3)=atom_gridsize*int(jobs(3)/atom_gridsize)

	close(35)

	if(spectrumanalysis .eq. 1)then
		allocate(ftr_in(0:ifmax),fti_in(0:ifmax),fte2_in(0:ifmax))
		if(TM .ge.1)then
			allocate(ftr_z(nobserver,0:ifmax),fti_z(nobserver,0:ifmax),fte2_z(nobserver,0:ifmax),fte2max_z(nobserver),  &
			ezinc(nobserver))
		endif
		if(TE .ge.1)then
			allocate(ftr_x(nobserver,0:ifmax),fti_x(nobserver,0:ifmax),fte2_x(nobserver,0:ifmax),fte2max_x(nobserver),  &
			ftr_y(nobserver,0:ifmax),fti_y(nobserver,0:ifmax),fte2_y(nobserver,0:ifmax),fte2max_y(nobserver),  &
			exinc(nobserver),eyinc(nobserver))
		endif
	endif
	
!----------------------------------------------------------------------------
!
!	Source pre-run: input pulse in time domain : 'sourceprerun.csv'
!								in frequency domain	: 'srcE-DFT.csv'
!
!----------------------------------------------------------------------------
if(spectrumanalysis .eq. 1)then

	!create source prerun filenames, like: sprfilename(1,1)='sourceprerun1_1.csv', etc.
	do 2352 u = 1,Nwgmax
	do 2352 v = 1,Nfieldmax
		rmf_ix1 = 48 + u
		rmf_ix2 = 48 + v
		sprfilename(u,v) = trim(datadir)//'sourceprerun'//char(rmf_ix1)//'_'//char(rmf_ix2)//'.csv' 
		dftfilename(u,v) = trim(datadir)//'srcE-DFT'//char(rmf_ix1)//'_'//char(rmf_ix2)//'.csv' 
2352 continue

	wstart=tpi*cc/lamdastart
	wspan=tpi*cc*(1./lamdaend-1./lamdastart)
	
!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
!DEC$ ENDIF ! (USE_MPI_CH) 

	do 2401 m=1,num_wg
		do 2402 k=1,nfield(m)
			numfile=3000+m*100+k
			open(numfile,file=trim(sprfilename(m,k)))
			
			tdelay=tdelay0+1.5*pwidth(m,k)
			ntspan=int(2.*tdelay/dt)

			! write temporal pulse shape
			do 20500 n=1, ntspan
				tt0=float(n)*dt - tdelay
				gaussian=exp(-tt0**2/(pwidth(m,k)/2)**2)
				sinusoidal=sin(w_f(m,k)*tt0)
				sinegauss=sinusoidal*gaussian
				write(numfile,*) 	tt0,',',sinegauss

			! Fourier transform of input pulse (real & imaginary amplitude)
			do 20501 iff=0,ifmax
				ffreq=wstart + float(iff)*wspan/float(ifmax)
				ftarg=ffreq*float(n)*dt
				ftcos=cos(ftarg)
				ftsin=sin(ftarg)
				ftr_in(iff) = ftr_in(iff) + E_input(m,k)*sinegauss*ftcos
				fti_in(iff) = fti_in(iff) + E_input(m,k)*sinegauss*ftsin
20501		continue
20500		continue
			close (numfile)

			! output the prerun spectrum of pulse and gain spectrum
			numfile2=4000+m*100+k
			open(numfile2,file=trim(dftfilename(m,k)))

			fte2max_in = 0.
				
			! Fourier transform (amplitude^2 & max-finder)
			do 30308 iff=0,ifmax				
				fte2_in(iff) = ftr_in(iff)**2 + fti_in(iff)**2
				if (fte2_in(iff) .gt. fte2max_in) then
				fte2max_in=fte2_in(iff)
				endif
30308		continue
				
			
			!Output Fourier spectrum into 'src-DFT(m).csv'
			do 30301 iff=0,ifmax
				ffreq=wstart + float(iff)*wspan/float(ifmax)
				wavelength0=tpi*cc/ffreq
				write(numfile2,8701)  wavelength0,',',ffreq/tpi,',',fte2_in(iff)/fte2max_in, &
				',',fte2_in(iff)
30301		continue

			close(numfile2)
2402	continue
2401 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF ! (USE_MPI_CH)

endif !spectrum analysis option

8701 format(e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3)
8702 format(e13.6e3,a3,e13.6e3,a3,e13.6e3,a3)

!.......Specify Electrical Properties......
	eps(1)=rns*rns
	eps(2)=rnf*rnf
	epsn=eps(2)*epsz
  mur(1)=1.0
  mur(2)=1.0   
 
!---------------------------------------------------------------------------
!
!	Output parameters into "src-INFO.csv"
!
!---------------------------------------------------------------------------
!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
!DEC$ ENDIF ! (USE_MPI_CH) 
!
	open(102,file= (trim(datadir)//'src-INFO.csv'),ACCESS='SEQUENTIAL')
	!FDTD parameters
	CALL FDATE(timestart)
	write(102,*) 'start time:', timestart
	write(102,*) 'FDTD parameters--------------------'
	write(102,*) 'dx',',',dx
	write(102,*) 'dy',',',dy
	write(102,*) 'dt',',',dt
	write(102,*) 'tmax',',',tmax
	write(102,*) 'nmax',',',nmax

	!Grid information
	write(102,*) 'Grid parameters-------------------'
	write(102,*) 'isize',',',isize
	write(102,*) 'jsize',',',jsize
	write(102,*) 'isize-actual length',',',isize*dx
	write(102,*) 'jsize_actual',',',jsize*dy
	write(102,*) 'iebc',',',iebc
	write(102,*) 'jebc',',',jebc
	write(102,*) 'rmax',',',rmax
	write(102,*) 'orderbc',',',orderbc
    write(102,*) 'mediabc',',',mediabc

	!Waveguide information
	write(102,*) 'Waveguide information-------------'
	write(102,*) 'total number of waveguide',',',num_wg
	write(102,*) 'number of waveguide in x-direction',',',num_wgx
	write(102,*) 'number of waveguide in y-direction',',',num_wgy
	do 244 m=1,num_wg
		write(102,*) 'waveguide #',m
		if(wgdirection(m) .eq. 1)then
			write(102,*) 'x-waveguide start',',',wgstart(m)
			write(102,*) 'x-waveguide end',',',wgend(m)
 			write(102,*) 'x-waveguide width',',',width(m),'=',width(m)*dx*1.e6,'um'
			write(102,*) 'x-waveguide index',',',wgindex(m)
		else
			write(102,*) 'y-waveguide start',',',wgstart(m)
			write(102,*) 'y-waveguide end',',',wgend(m)
 			write(102,*) 'y-waveguide width',',',width(m),'=',width(m)*dy*1.e6,'um'
			write(102,*) 'y-waveguide index',',',wgindex(m)
		endif
244	continue

	!Atom-field interaction
	write(102,*) 'Atom-Field interaction information------------'
	write(102,*) 'interaction',',',interaction   
	write(102,*) 'ipauli',',',ipauli

	!Inputfield information
	write(102,*) 'Inputfield information------------'
	write(102,*) 'nfield(1)',',',nfield(1)
	write(102,*) 'nfield(2)',',',nfield(2)
	write(102,*) 'nfield(3)',',',nfield(3)
	write(102,*) 'lamda_f(1,1)',',',lamda_f(1,1)   
	write(102,*) 'lamda_f(1,2)',',',lamda_f(1,2) 
	write(102,*) 'lamda_f(1,3)',',',lamda_f(1,3)   
	write(102,*) 'lamda_f(1,4)',',',lamda_f(1,4)    
	write(102,*) 'lamda_f(2,1)',',',lamda_f(2,1) 
	write(102,*) 'lamda_f(2,2)',',',lamda_f(2,2) 
	write(102,*) 'lamda_f(2,3)',',',lamda_f(2,3) 
	write(102,*) 'lamda_f(2,4)',',',lamda_f(2,4)   
	write(102,*) 'lamda_f(3,1)',',',lamda_f(3,1)   
	write(102,*) 'lamda_f(3,2)',',',lamda_f(3,2) 
	write(102,*) 'lamda_f(3,3)',',',lamda_f(3,3)   
	write(102,*) 'lamda_f(3,4)',',',lamda_f(3,4) 
	write(102,*) 'pol(1,1)',',',pol(1,1)
	write(102,*) 'pol(1,2)',',',pol(1,2)
	write(102,*) 'pol(1,3)',',',pol(1,3)
	write(102,*) 'pol(1,4)',',',pol(1,4)
	write(102,*) 'pol(2,1)',',',pol(2,1)
	write(102,*) 'pol(2,2)',',',pol(2,2)
	write(102,*) 'pol(2,3)',',',pol(2,3)
	write(102,*) 'pol(2,4)',',',pol(2,4)
	write(102,*) 'pol(3,1)',',',pol(3,1)
	write(102,*) 'pol(3,2)',',',pol(3,2)
	write(102,*) 'pol(3,3)',',',pol(3,3)
	write(102,*) 'pol(3,4)',',',pol(3,4)
	write(102,*) 'TE',',',TE
	write(102,*) 'TM',',',TM
	write(102,*) 'w_f(1,1)',',',w_f(1,1) 
	write(102,*) 'w_f(1,2)',',',w_f(1,2)
	write(102,*) 'w_f(1,3)',',',w_f(1,3) 
	write(102,*) 'w_f(1,4)',',',w_f(1,4)  
	write(102,*) 'w_f(2,1)',',',w_f(2,1) 
	write(102,*) 'w_f(2,2)',',',w_f(2,2) 
	write(102,*) 'w_f(2,3)',',',w_f(2,3) 
	write(102,*) 'w_f(2,4)',',',w_f(2,4) 
	write(102,*) 'w_f(3,1)',',',w_f(3,1)  
	write(102,*) 'w_f(3,2)',',',w_f(3,2)
	write(102,*) 'w_f(3,3)',',',w_f(3,3)  
	write(102,*) 'w_f(3,4)',',',w_f(3,4) 
	write(102,*) 'inputfield(1,1) (0:CW, 1:pulse)',',',inputfield(1,1)	 
	write(102,*) 'inputfield(1,2) (0:CW, 1:pulse)',',',inputfield(1,2)
	write(102,*) 'inputfield(1,3) (0:CW, 1:pulse)',',',inputfield(1,3)	 
	write(102,*) 'inputfield(1,4) (0:CW, 1:pulse)',',',inputfield(1,4)		 
	write(102,*) 'inputfield(2,1) (0:CW, 1:pulse)',',',inputfield(2,1)	 
	write(102,*) 'inputfield(2,2) (0:CW, 1:pulse)',',',inputfield(2,2)
	write(102,*) 'inputfield(2,3) (0:CW, 1:pulse)',',',inputfield(2,3)	 
	write(102,*) 'inputfield(2,4) (0:CW, 1:pulse)',',',inputfield(2,4)	 
	write(102,*) 'inputfield(3,1) (0:CW, 1:pulse)',',',inputfield(3,1)	 
	write(102,*) 'inputfield(3,2) (0:CW, 1:pulse)',',',inputfield(3,2)
	write(102,*) 'inputfield(3,3) (0:CW, 1:pulse)',',',inputfield(3,3)	 
	write(102,*) 'inputfield(3,4) (0:CW, 1:pulse)',',',inputfield(3,4)	 

	write(102,*) 'I_input_cgs(1,1)',',',I_input_cgs(1,1)   
	write(102,*) 'I_input_cgs(1,2)',',',I_input_cgs(1,2) 
	write(102,*) 'I_input_cgs(1,3)',',',I_input_cgs(1,3)   
	write(102,*) 'I_input_cgs(1,4)',',',I_input_cgs(1,4)
	write(102,*) 'I_input_cgs(2,1)',',',I_input_cgs(2,1)
	write(102,*) 'I_input_cgs(2,2)',',',I_input_cgs(2,2)
	write(102,*) 'I_input_cgs(2,3)',',',I_input_cgs(2,3)
	write(102,*) 'I_input_cgs(2,4)',',',I_input_cgs(2,4)    
	write(102,*) 'I_input_cgs(3,1)',',',I_input_cgs(3,1) 
	write(102,*) 'I_input_cgs(3,2)',',',I_input_cgs(3,2)
	write(102,*) 'I_input_cgs(3,3)',',',I_input_cgs(3,3) 
	write(102,*) 'I_input_cgs(3,4)',',',I_input_cgs(3,4)
	write(102,*) 'I_input(1,1): mks',',',I_input(1,1)   
	write(102,*) 'I_input(1,2): mks',',',I_input(1,2) 
	write(102,*) 'I_input(1,3): mks',',',I_input(1,3)   
	write(102,*) 'I_input(1,4): mks',',',I_input(1,4) 
	write(102,*) 'I_input(2,1): mks',',',I_input(2,1)
	write(102,*) 'I_input(2,2): mks',',',I_input(2,2)
	write(102,*) 'I_input(2,3): mks',',',I_input(2,3)
	write(102,*) 'I_input(2,4): mks',',',I_input(2,4)    
	write(102,*) 'I_input(3,1): mks',',',I_input(3,1) 
	write(102,*) 'I_input(3,2): mks',',',I_input(3,2)
	write(102,*) 'I_input(3,3): mks',',',I_input(3,3) 
	write(102,*) 'I_input(3,4): mks',',',I_input(3,4)
	write(102,*) 'E_input(1,1) (V/m)',',',E_input(1,1)   
	write(102,*) 'E_input(1,2) (V/m)',',',E_input(1,2)
	write(102,*) 'E_input(1,3) (V/m)',',',E_input(1,3)   
	write(102,*) 'E_input(1,4) (V/m)',',',E_input(1,4)  
	write(102,*) 'E_input(2,1) (V/m)',',',E_input(2,1) 
	write(102,*) 'E_input(2,2) (V/m)',',',E_input(2,2) 
	write(102,*) 'E_input(2,3) (V/m)',',',E_input(2,3) 
	write(102,*) 'E_input(2,4) (V/m)',',',E_input(2,4)  
	write(102,*) 'E_input(3,1) (V/m)',',',E_input(3,1) 
	write(102,*) 'E_input(3,2) (V/m)',',',E_input(3,2) 
	write(102,*) 'E_input(3,3) (V/m)',',',E_input(3,3) 
	write(102,*) 'E_input(3,4) (V/m)',',',E_input(3,4) 
	write(102,*) 'pwidth(1,1)',',',pwidth(1,1)	
	write(102,*) 'pwidth(1,2)',',',pwidth(1,2)
	write(102,*) 'pwidth(1,3)',',',pwidth(1,3)	
	write(102,*) 'pwidth(1,4)',',',pwidth(1,4)
	write(102,*) 'pwidth(2,1)',',',pwidth(2,1)
	write(102,*) 'pwidth(2,2)',',',pwidth(2,2)
	write(102,*) 'pwidth(2,3)',',',pwidth(2,3)
	write(102,*) 'pwidth(2,4)',',',pwidth(2,4)	
	write(102,*) 'pwidth(3,1)',',',pwidth(3,1)
	write(102,*) 'pwidth(3,2)',',',pwidth(3,2)
	write(102,*) 'pwidth(3,3)',',',pwidth(3,3)
	write(102,*) 'pwidth(3,4)',',',pwidth(3,4)
	write(102,*) 'Ilaunch(1,1)',',',Ilaunch(1,1)	 
	write(102,*) 'Ilaunch(1,2)',',',Ilaunch(1,2)
	write(102,*) 'Ilaunch(1,3)',',',Ilaunch(1,3)	 
	write(102,*) 'Ilaunch(1,4)',',',Ilaunch(1,4)
	write(102,*) 'Ilaunch(2,1)',',',Ilaunch(2,1)
	write(102,*) 'Ilaunch(2,2)',',',Ilaunch(2,2)
	write(102,*) 'Ilaunch(2,3)',',',Ilaunch(2,3)
	write(102,*) 'Ilaunch(2,4)',',',Ilaunch(2,4)	 
	write(102,*) 'Ilaunch(3,1)',',',Ilaunch(3,1)
	write(102,*) 'Ilaunch(3,2)',',',Ilaunch(3,2)
	write(102,*) 'Ilaunch(3,3)',',',Ilaunch(3,3)
	write(102,*) 'Ilaunch(3,4)',',',Ilaunch(3,4)
	
	!Spectrum analysis for input pulse 
	write(102,*) 'Specrtrum Analysis parameters------------'
	write(102,*) 'lamdastart',',', lamdastart
	write(102,*) 'lamdaend',',', lamdaend
	write(102,*) 'wstart',',',wstart
	write(102,*) 'wspan',',',wspan
	write(102,*) 'ifmax',',',ifmax			
	write(102,*) 'lamda resolution',',', (tpi*cc/wstart-tpi*cc/(wstart+wspan))/ifmax
	write(102,*) 'tdelay0',',',tdelay0

	!Spectrum analysis for the output
	write(102,*) 'spectrumanalysis',',',spectrumanalysis  ! 1:yes, 2: no

	!Snapshot parameters
	write(102,*) 'Snapshot information------------'
	write(102,*) 'snapshotinterval',',',snapshotinterval	  
	write(102,*) 'reso_snapshot_x',',',reso_snapshot_x	  
	write(102,*) 'reso_snapshot_x',',',reso_snapshot_x	  
	write(102,*) 'outputchoice',',', outputchoice

	!Detector related
	write(102,*) 'Detector parameters------------'
	write(102,*) 'Number of detectors',',',det_num
	write(102,*) 'Number of spectrum detectors',',',ftdet_num
	do 215 k=1,det_num
		write(102,*) 'Detector length',',',det_xlen(k),det_ylen(k)
215	continue
	write(102,*) 'maximum detector length',',',det_maxlen
	write(102,*) 'detector period',',',period
	
	!Metal related
	write(102,*) 'number of metallic materials',',',metal

	!Gain related
	write(102,*) 'Gain parameters-----------------'
	write(102,*) 'gain?',',',gain
	if(gain .gt. 0)then
		write(102,*) 'number of atomic groups:',',',atomgroup_number
		write(102,*) 'skipping factor in atomic distribution:',',',atom_gridsize

		do 	9877 i=1,atomgroup_number
			write(102,*) 'lamda_a1(',i,')',',', lamda_a1(i)		!atomic resonance wavelegnth between 1&2
			write(102,*) 'lamda_a2(',i,')',',', lamda_a2(i)		!atomic resonance wavelength between 0&3
			write(102,*) 'tau30(',i,')',',', tau30(i)			!radiative decay rate 3 -> 0
			write(102,*) 'tau21(',i,')',',', tau21(i)			!radiative decay rate 2 -> 1
			write(102,*) 'tau10(',i,')',',', tau10(i)			!non-radiative decay rate 1 -> 0
			write(102,*) 'tau32(',i,')',',', tau32(i)			!non-radiative decay rate 3 ->2
			write(102,*) 'dwa1(',i,')',',', dwa1(i)
			write(102,*) 'dwa2(',i,')',',', dwa2(i)
			write(102,*) 'Ndensity(',i,')',',', Ndensity(i)		!Numer of atoms per volume in cgs
			write(102,*) 'ne(',i,')',',', ne(i)	   			!number of electrons
			write(102,*) 'n0_0(',i,')',',',	n0_0(i) 			!initial value of n0
			write(102,*) 'n1_0(',i,')',',', n1_0(i) 			!initial value of n1
			write(102,*) 'n2_0(',i,')',',', n2_0(i) 			!initial value of n2
			write(102,*) 'n3_0(',i,')',',',	n3_0(i) 			!initial value of n3
			write(102,*) 'current_pumiping(',i,')',',',	current_pumping(i) 
9877	continue
	endif !gain.gt.0

	!options
	write(102,*) 'Other options------------'
	write(102,*) 'waveguide mode read from file',',',read_bpm	  
	write(102,*) 'mpi segment read from file',',',read_mpi_seg	  
	write(102,*) 'diagnostic output',',',diag	
	write(102,*) 'flux snapshot',',',fluxmap
	write(102,*) 'flux snapshot taken at t = ',',',t_fluxmap
	write(102,*) 'E(t) output',',',tempdata
	write(102,*) 'Specify input waveguide position?',',',inputwg
	write(102,*) 'Input waveguide position',',',j_inputwg

	close(102)  ! information append to this file later


!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF ! (USE_MPI_CH) 
!
!---------------------------------------------------------------------
!
!					Fill UPML Grid
!
!---------------------------------------------------------------------
      
	  !TM-initialization
	  do 8200 i=1,iefbc
			axef(i)=1.0
			axeb(i)=1.0
			bxef(i)=1.0
			bxeb(i)=1.0
			cxef(i)=1.0
8200  cxeb(i)=1.0

    do 8201 j=1,jebc
			ayef(j)=1.0
			ayeb(j)=1.0
			byef(j)=1.0
			byeb(j)=1.0
			cyef(j)=1.0
8201  cyeb(j)=1.0

    do 8202 i=1,iebc
      axel(i)=1.0
			axer(i)=1.0
			bxel(i)=1.0
			bxer(i)=1.0
			cxel(i)=1.0
8202  cxer(i)=1.0

    do 8207 j=1,je
			ayel(j)=1.0
			ayer(j)=1.0
			byel(j)=1.0
			byer(j)=1.0
			cyel(j)=1.0
8207  cyer(j)=1.0

	  do 8203 i=1,ibfbc
      axhf(i)=1.0
			axhb(i)=1.0
			bxhf(i)=1.0
			bxhb(i)=1.0
			cxhf(i)=1.0
8203  cxhb(i)=1.0

	  do 8204 j=1,jbbc
      ayhf(j)=1.0
			ayhb(j)=1.0
			byhf(j)=1.0
			byhb(j)=1.0
			cyhf(j)=1.0
8204  cyhb(j)=1.0

    do 8205 i=1,ibbc
      axhl(i)=1.0
			axhr(i)=1.0
			bxhl(i)=1.0
			bxhr(i)=1.0
			cxhl(i)=1.0
8205  cxhr(i)=1.0

	  do 8206 j=1,jb
			ayhl(j)=1.0
			ayhr(j)=1.0
			byhl(j)=1.0
			byhr(j)=1.0
			cyhl(j)=1.0
8206  cyhr(j)=1.0

	!-----------------------------------------------------------------------
	!		UMPL coefficients definition
	!-----------------------------------------------------------------------
	  kmax=1.
      do 8230 m=1,mediabc
      delbc=float(iebc)*dx			!actual thickness of absorbing boundary

	  eta=sqrt(muz*mur(m)/epsz/eps(m))	  
      sigmam=-log(rmax)*(orderbc+1.0)/(2.0*delbc*eta*eps(m))
		bcfactor=sigmam/(dx*(delbc**orderbc)*(orderbc+1.0))
		kfactor=(Kmax-1.)/(dx*(delbc**orderbc)*(orderbc+1.0)) 

        do 824 i=1,iebc
          x1=dx*float(i)
          x2=dx*(float(i)-1.0)
          sigma=bcfactor*(x1**(orderbc+1.0)-x2**(orderbc+1.0))
					kapa=1.+kfactor*(x1**(orderbc+1.0)-x2**(orderbc+1.0))

          axebc(i,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz)
          ayebc(i,m)= axebc(i,m)
					bxebc(i,m)=(kapa+sigma*dt/2.0/epsz)
					byebc(i,m)= bxebc(i,m)
					cxebc(i,m)=(kapa-sigma*dt/2.0/epsz)
					cyebc(i,m)= cxebc(i,m)
824		continue

		x1=dx*0.50
		sigma=bcfactor*x1**(orderbc+1.0)
		kapa=1.+kfactor*x1**(orderbc+1.0)

		axhbc(1,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz)
		ayhbc(1,m)= axhbc(1,m)
		bxhbc(1,m)=(kapa+sigma*dt/2.0/epsz)
		byhbc(1,m)= bxhbc(1,m)
		cxhbc(1,m)=(kapa-sigma*dt/2.0/epsz)
		cyhbc(1,m)= cxhbc(1,m)
        
        do 825 i=2,iebc
          x1=dx*(float(i)-0.50)
          x2=dx*(float(i)-1.50)
          sigma=bcfactor*(x1**(orderbc+1.0)-x2**(orderbc+1.0))
		  kapa=1.+kfactor*(x1**(orderbc+1.0)-x2**(orderbc+1.0))

			axhbc(i,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz)
			ayhbc(i,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz)
			bxhbc(i,m)=(kapa+sigma*dt/2.0/epsz)
			byhbc(i,m)=(kapa+sigma*dt/2.0/epsz)
			cxhbc(i,m)=(kapa-sigma*dt/2.0/epsz)
			cyhbc(i,m)=(kapa-sigma*dt/2.0/epsz)
825		continue 

8230	continue


	!------------------------------------------------------------------------
	!		Assignment of UPML coefficients
	!------------------------------------------------------------------------
	  do 8260 i=1,iebc
		 if1=iebc+1-i
	   if2=ie+iebc+i
		 ill=iebc+1-i
		 ir=i
	   axef(if1) =axebc(i,1)
		 axef(if2) =axebc(i,1)
		 axeb(if1) =axebc(i,1)
		 axeb(if2) =axebc(i,1)
		 axel(ill) =axebc(i,1)
		 axer(ir)  =axebc(i,1)

		 bxef(if1) =bxebc(i,1)
		 bxef(if2) =bxebc(i,1)
		 bxeb(if1) =bxebc(i,1)
		 bxeb(if2) =bxebc(i,1)
		 bxel(ill) =bxebc(i,1)
		 bxer(ir)  =bxebc(i,1)

		 cxef(if1) =cxebc(i,1)
		 cxef(if2) =cxebc(i,1)
		 cxeb(if1) =cxebc(i,1)
		 cxeb(if2) =cxebc(i,1)
		 cxel(ill) =cxebc(i,1)
		 cxer(ir)  =cxebc(i,1)

		 axhf(if1+1) =axhbc(i,1)
		 axhf(if2)   =axhbc(i,1)
		 axhb(if1+1) =axhbc(i,1)
		 axhb(if2)   =axhbc(i,1)
		 axhl(ill+1) =axhbc(i,1)
		 axhr(ir)    =axhbc(i,1)
		 
		 bxhf(if1+1) =bxhbc(i,1)
		 bxhf(if2)   =bxhbc(i,1)
		 bxhb(if1+1) =bxhbc(i,1)
		 bxhb(if2)   =bxhbc(i,1)
		 bxhl(ill+1) =bxhbc(i,1)
		 bxhr(ir)    =bxhbc(i,1)

		 cxhf(if1+1) =cxhbc(i,1)
		 cxhf(if2)   =cxhbc(i,1)
		 cxhb(if1+1) =cxhbc(i,1)
		 cxhb(if2)   =cxhbc(i,1)
		 cxhl(ill+1) =cxhbc(i,1)
		 cxhr(ir)    =cxhbc(i,1)

   8260 continue

	   do 8261 j=1,jebc
	   jf=jebc+1-j
		 
		 ayef(jf)  =ayebc(j,1)
		 ayeb(j)   =ayebc(j,1)
		 
		 byef(jf)  =byebc(j,1)
		 byeb(j)   =byebc(j,1)
		 
		 cyef(jf)  =cyebc(j,1)
		 cyeb(j)   =cyebc(j,1)
		 	 
		 ayhf(jf+1) =ayhbc(j,1)
		 ayhb(j)    =ayhbc(j,1)
		 
		 byhf(jf+1) =byhbc(j,1)
		 byhb(j)    =byhbc(j,1)
		 
		 cyhf(jf+1) =cyhbc(j,1)
		 cyhb(j)    =cyhbc(j,1)
		 
  8261 continue

!----------------------------------------------------------------------------
!			Initialize (Ez,Hx,Hy) in the main grid
!---------------------------------------------------------------------------- 
!write(*,*) 'Initialize (Ez,Hx,Hy) in the main grid'
		!if(TM .gt. 0)then !TM-mode
		Ez = 0.
		Hx = 0.
		Hy = 0.
			!medium
			if(gain .gt. 0) then
			pz_1=0.
			pz_2=0.
			pz_old1_1=0.
			pz_old1_2=0.

			endif
			!end medium
		!endif

		if(TE .gt. 0)then !TE-mode
		Ex = 0.	
		Ey = 0.
		Hz = 0.
			!medium
			if(gain .gt. 0) then
			px_1=0.
			px_2=0.
			px_old1_1=0.
			px_old1_2=0.
			py_1=0.
			py_2=0.
			py_old1_1=0.
			py_old1_2=0.
			endif
			!end medium
		endif

		!medium
		if(gain .gt. 0) then
			n0=0.
			n1=0.
			n2=0.
			n3=0.
!DEC$ IF DEFINED (USE_MPI_CH)
			do 8050 i=1,ie
			do 8050 j=j1_al,j2_al
				if (gaingeo(i,j) .gt. 0) then
				p=gaingeo(i,j)
				n0(i,j,p)=n0_0(p)	 ! the population inversion number set to be 1e20 to 1e23
				n1(i,j,p)=n1_0(p)
				n2(i,j,p)=n2_0(p)
				n3(i,j,p)=n3_0(p)
				endif
8050		continue
!DEC$ ELSE
			do 8051 i=1,ie
			do 8051 j=1,je
				if (gaingeo(i,j) .gt. 0) then
				p=gaingeo(i,j)
				n0(i,j,p)=n0_0(p)	 ! the population inversion number set to be 1e20 to 1e23
				n1(i,j,p)=n1_0(p)
				n2(i,j,p)=n2_0(p)
				n3(i,j,p)=n3_0(p)
				endif
8051		continue
!DEC$ ENDIF ! (USE_MPI_CH)
		endif
		!end medium
!------------------------------------------------------------------------------
!			Initialize (Ez,Dz,Hx,Hy,Bx,By in the UPML grid)
!------------------------------------------------------------------------------
!write(*,*) 'Initialize (Ez,Dz,Hx,Hy,Bx,By in the UPML grid)'
	! front & back BC
		!if(TM .gt. 0)then !TM-mode
		ezf=0.0
		ezb=0.0
		dezf=0.0
		dezb=0.0
		!endif

		if(TE .gt. 0)then !TE-mode
		hzf=0.0
		hzb=0.0
		bhzf=0.0
		bhzb=0.0
		endif 

		!if(TM .gt. 0)then 
		hxf = 0.0
		hxb = 0.0
		Bhxf=0.0
		Bhxb=0.0
		!endif
		
		if(TE .gt. 0)then !TE-mode
		exf = 0.0
		exb = 0.0
		dexf = 0.0
		dexb = 0.0
		endif	

		!if(TM .gt. 0)then !TM-mode
		hyf = 0.0
		hyb = 0.0
		Bhyf=0.0
		Bhyb=0.0
		!endif
		
		if(TE .gt. 0)then !TE-mode
		eyf = 0.0
		eyb = 0.0
		deyf=0.0
		deyb=0.0
		endif 

	! left & right BC
		!if(TM .gt. 0)then !TM-mode
		ezbcl=0.0
		ezbcr=0.0
		dezbcl=0.0
		dezbcr=0.0
		!endif
		
		if(TE .gt. 0)then !TE-mode
		hzbcl=0.0
		hzbcr=0.0
		bhzbcl=0.0
		bhzbcr=0.0
		endif	

		!if(TM .gt. 0)then !TM-mode
		hxbcl = 0.0
		hxbcr = 0.0
		Bhxbcl=0.0
		Bhxbcr=0.0
		!endif
		
		if(TE .gt. 0)then !TE-mode
		exbcl = 0.0
		exbcr = 0.0
		dexbcl=0.0
		dexbcr=0.0
		endif	

		!if(TM .gt. 0)then !TM-mode
		hybcl = 0.0
		hybcr = 0.0
		Bhybcl=0.0
		Bhybcr=0.0
		!endif
		
		if(TE .gt. 0)then !TE-mode
		eybcl = 0.0
		eybcr = 0.0
		deybcl=0.0
		deybcr=0.0
		endif	


!-------- file open for time step output ---------------
!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. master ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 
!
	if (tempdata .eq. 1) then 
	!medium
	if (gain .gt. 0) then
 		open(123,access='sequential',file= (trim(datadir)//'observerN0.csv'))
		open(124,access='sequential',file= (trim(datadir)//'observerN1.csv'))
		open(125,access='sequential',file= (trim(datadir)//'observerN2.csv'))
		open(126,access='sequential',file= (trim(datadir)//'observerN3.csv'))
	endif
	!end medium
 	if(TE .eq. 0)then
		open(122,access='sequential',file= (trim(datadir)//'observerTM.csv'))	!file for Ez(t) 
	elseif(TM .eq. 0)then
		open(122,access='sequential',file= (trim(datadir)//'observerTE.csv'))
	else
		open(122,access='sequential',file= (trim(datadir)//'observerTM_TE.csv'))
	endif
	!
	endif !tempdata option
!
!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

! allocate and initialize observer arrays.. 

	!medium
	if (gain .gt. 0) then
 	allocate(n0_obs(ns_obs,nobserver),n1_obs(ns_obs,nobserver),n2_obs(ns_obs,nobserver),n3_obs(ns_obs,nobserver))
	endif
	!end medium

 	if(TE .eq. 0)then
			allocate(Ez_obs(ns_obs,nobserver))
			Ez_obs(:,:) = 0.
	elseif(TM .eq. 0)then
			allocate(Ex_obs(ns_obs,nobserver), Ey_obs(ns_obs,nobserver))
			Ex_obs(:,:) = 0.
			Ey_obs(:,:) = 0.
	else
			allocate(Ez_obs(ns_obs,nobserver))
			allocate(Ex_obs(ns_obs,nobserver), Ey_obs(ns_obs,nobserver))
			Ez_obs(:,:) = 0.
			Ex_obs(:,:) = 0.
			Ey_obs(:,:) = 0.
	endif
!
!DEC$ IF DEFINED (USE_MPI_CH) 
! identify observer points that fall on the grid of this process:
	do 306 pid_g = 0, (numprocs-1)
	do 306 u = 1,nobserver
		if( je_np(pid_g) .lt. jobs(u) .and. jobs(u) .le. je_np(pid_g+1) ) then
			obs_in_proc(u) = pid_g						
		endif
306 continue
!DEC$ ENDIF ! (USE_MPI_CH) 
!!
	  ix=0
	  iy=0
	  iz=0
	  logix=0
	  logiy=0
	  logiz=0
	  exmax=0
	  exmin=0
	  ixmin=0
	  ixmax=0
	  logixmin=0
	  logixmax=0
	  eymax=0
	  eymin=0
	  iymin=0
	  iymax=0
	  logiymin=0
	  logiymax=0
	  ezmax=0
	  ezmin=0
	  izmin=0
	  izmax=0
	  logiz=0
	  logizmin=0
	  logizmax=0

	det_result(:)=0.
!DEC$ IF DEFINED (USE_MPI_CH) 
	det_result_t(:)=0.
!DEC$ ENDIF ! (USE_MPI_CH) 

!---------------------------------------------
!        Coefficient definition
!---------------------------------------------
	cb=dt/epsz/dx
	db=dt/muz/dx

	do 7654	m=1,atomgroup_number
	pa1_1(m)=2.*dt*dt/(2.+dWa1(m)*dt)
	pa2_1(m)=(2-(wa1(m))**2*dt*dt)*2./(2.+dWa1(m)*dt)  !(2.+dwa1*dt) !(2./(dt**2)-Wa1**2)
	pa3_1(m)=(dwa1(m)*dt-2.)/(2.+dWa1(m)*dt) !(dWa1/2./dt-1./(dt**2))

	pa1_2(m)=2.*dt*dt/(2.+dWa2(m)*dt)
	pa2_2(m)=(2./dt/dt-(wa2(m))**2)*pa1_2(m)  !(2.+dwa1*dt) !(2./(dt**2)-Wa1**2)
	pa3_2(m)=(dwa2(m)/dt/2.-1./dt/dt)*pa1_2(m) !(dWa1/2./dt-1./(dt**2))

 	
    C_n2_1(m)=dt*wa2(m)*h_bar_inv
    C_n1_1(m)=dt*wa1(m)*h_bar_inv
	C_n2_2(m)=1./(wa2(m)*h_bar)
	C_n1_2(m)=1./(wa1(m)*h_bar)
	pumping_rate(m)=dt*current_pumping(m)
7654 continue

!---------------------------------------------
!        adjusting boundaries for mpich
!---------------------------------------------
!write(*,*) 'adjusting boundaries for mpich'
! no mpi (default):
	j1c = 1
	j2c = je
	i1c = 1		! center grid left
	i2c = ie	! center grid right
	js1 = 4
	js2 = je-4
	is1 = 4
	is2 = ie-4
	j1b = 2
	i1b = 2
!
!DEC$ IF DEFINED (USE_MPI_CH) 
	j1c = je_np(pid)+1
	j2c = je_np(pid+1)
	js1 = j1c ! j1c for '1' .. 'np-1'; 4, for '0'; np - number of processes; 
	js2 = j2c ! j2c for '0' .. 'np-2'; je-4, for 'np-1';
	j1b = j1c ! j1c for '1' .. 'np-1'; 2, for '0' 

	if( pid .eq. 0 ) then
		js1 = 4
		j1b = 2 
	elseif( pid .eq. (numprocs-1) ) then
		js2 = je-4
	endif

    call MPI_BARRIER( comm_a, ierr ) ! synchronize processes;

		if( pid .eq. master ) then
			open(1345,file= (trim(datadir)//'detector.csv'))		
		endif

!DEC$ ELSE ! (USE_MPI_CH) 
	open(1345,file= (trim(datadir)//'detector.csv'))		
!DEC$ ENDIF ! (USE_MPI_CH) 

	n_st = 1 ! initialize observer index

!Initialize fft matrix for spectrum detector
	if(ftdet_num .ge. 1)then
				if(TM .gt. 0)then
				ftrd_z(0:ifmax,0:det_maxlen,1:det_num)=0.
				ftid_z(0:ifmax,0:det_maxlen,1:det_num)=0.
				fte2d_z(0:ifmax,0:det_maxlen,1:det_num)=0.
				fte2d_zt(0:ifmax,0:det_maxlen,1:det_num)=0.
				fte2maxd_z(1:det_num)=0.
				endif
				if(TE .gt. 0)then
				ftrd_x(0:ifmax,0:det_maxlen,1:det_num)=0.
				ftid_x(0:ifmax,0:det_maxlen,1:det_num)=0.
				fte2d_x(0:ifmax,0:det_maxlen,1:det_num)=0.
				fte2d_xt(0:ifmax,0:det_maxlen,1:det_num)=0.
				ftrd_y(0:ifmax,0:det_maxlen,1:det_num)=0.
				ftid_y(0:ifmax,0:det_maxlen,1:det_num)=0.
				fte2d_y(0:ifmax,0:det_maxlen,1:det_num)=0.
				fte2d_yt(0:ifmax,0:det_maxlen,1:det_num)=0.
				fte2maxd_x(1:det_num)=0.
				fte2maxd_y(1:det_num)=0.
				endif
	endif

!fluxmap related
	if(fluxmap .eq. 1)then
		allocate(s_x_f_d(ie,j1c:j2c), s_y_f_d(ie,j1c:j2c), s_f_d(ie,j1c:j2c))
		s_x_f_d(:,:) = 0.
		s_y_f_d(:,:) = 0.
		s_f_d(:,:) = 0.

		num_ts_fd = int(det_period/dt) + 1 ! = lamda_f(1,1)/dx = (f.e.) 1.55/0.1 = 15.5
		n_fluxmap = int(t_fluxmap/dt)
	endif

!-------------------------------------------------------------------
!
!					TIME-STEPPING LOOP
!
!-------------------------------------------------------------------
do 500 n=1,nmax

!write(*,*) "inside the time loop"
	tt=float(n)*dt

!DEC$ IF DEFINED (USE_MPI_CH) 
    call MPI_BARRIER( comm_a, ierr ) ! synchronize entrance;
		if( pid .eq. master ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 
!
	if( mod(n,10) .eq. 0 ) then
	write(*,'(I8,$)') n ! if pid .eq. master 
	endif
!
!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 
  
	!--------------------------------------------------------------------
	!		Ez UPDATE in the MAIN grid      
	!-------------------------------------------------------------------- 

! Ez (Hz) update on the main (center) grid; 
	do 101 j=j1c,j2c
	do 101 i=i1c,i2c
		geometryij=geometry(i,j)
	
	if (gain.gt.0) then	
		ez_old(i,j)=ez(i,j)
		gainij=gaingeo(i,j)
	else
		gainij=0
	endif
				

    if (gainij.eq.0)then
		!if(TM .gt. 0)then
			if(metal .eq. 0)then
				Ez(i,j) =Ez(i,j)+cb*geometryij*(Hy(i+1,j)-Hy(i,j)+Hx(i,j)-Hx(i,j+1))
			else
				caij=(1-sigmageo(i,j)*dt/2./epsz*geometryij)/(1+sigmageo(i,j)*dt/2./epsz*geometryij)
				cbij=cb/(1+sigmageo(i,j)*dt/2./epsz*geometryij)
				Ez(i,j)=caij*Ez(i,j)+cbij*geometryij*(Hy(i+1,j)-Hy(i,j)+Hx(i,j)-Hx(i,j+1))
			endif
		!endif

		if(TE .gt. 0)then
		Hz(i,j) =Hz(i,j)-db*(Ey(i+1,j)-Ey(i,j)+Ex(i,j)-Ex(i,j+1))
		endif
   	else
	!medium
		p=gainij
		n3_old(i,j,p)=n3(i,j,p)
		n2_old(i,j,p)=n2(i,j,p)
		n1_old(i,j,p)=n1(i,j,p)
		n0_old(i,j,p)=n0(i,j,p)

		!if (TM .gt.0) then
			Pz_old2_1(p)=Pz_old1_1(i,j,p)
			Pz_old1_1(i,j,p)=Pz_1(i,j,p)

			Pz_old2_2(p)  =Pz_old1_2(i,j,p)
			Pz_old1_2(i,j,p)=Pz_2(i,j,p)

	  		if(interaction .eq.1) then
				Pz_1(i,j,p) = kapa1(p)*(n1(i,j,p)-n2(i,j,p))*Ez(i,j) &
			 +pa2_1(p)*Pz_old1_1(i,j,p)  &
			 +pa3_1(p)*pz_old2_1(p)
				Pz_2(i,j,p) = kapa2(p)*(n0(i,j,p)-n3(i,j,p))*Ez(i,j) &
			 +pa2_2(p)*Pz_old1_2(i,j,p)  &
			 +pa3_2(p)*pz_old2_2(p)
			else
				Pz_1(i,j,p) = kapa1(p)*(n1(i,j,p)-n2(i,j,p))*Ez(i,j) &
			 +(pa2_1(p)-rabico1(p)*(Avptz(i,j))**2)*Pz_old1_1(i,j,p) &
			 +pa3_1(p)*pz_old2_1(p)
				Pz_2(i,j,p) = kapa2(p)*(n0(i,j,p)-n3(i,j,p))*Ez(i,j) &
			 +(pa2_2(p)-rabico2(p)*(Avptz(i,j))**2)*Pz_old1_2(i,j,p) &
			 +pa3_2(p)*pz_old2_2(p)
			endif
		!endif
 		!Electric field  
		
		!if (TM .gt. 0) then		  
			ez(i,j) =ez(i,j)+cb*geometryij*(Hy(i+1,j)-Hy(i,j)+Hx(i,j)-Hx(i,j+1))
			ez(i,j) =ez(i,j)-geometryij/epsz*Ndensity(p)*(Pz_1(i,j,p) &
				+Pz_2(i,j,p)-Pz_old1_1(i,j,p)-Pz_old1_2(i,j,p))
		!endif

		if(TE .gt. 0)then
			hz(i,j) =hz(i,j)-db*(ey(i+1,j)-ey(i,j)+ex(i,j)-ex(i,j+1))
		endif

		!Vector potential
		if(interaction .eq. 2)then
			Avptz_old(i,j)=Avptz(i,j)
			Avptz(i,j)=Avptz(i,j)-dt/2.*(Ez(i,j)+Ez_old(i,j))
		endif
			
		!Population
		!TM-mode
		if(TE .eq. 0)then
			if(interaction .eq. 1)then
				if (ipauli .eq. 2)then
			EdP1ij(p)=(Ez(i,j)+Ez_old(i,j))/2*(Pz_1(i,j,p)-Pz_old1_1(i,j,p))*C_n1_2(p)
			EdP2ij(p)=(Ez(i,j)+Ez_old(i,j))/2*(Pz_2(i,j,p)-Pz_old1_2(i,j,p))*C_n2_2(p)

			n3(i,j,p)=n3(i,j,p)*(1.-dt_over_tau30(p)*(1.-n0(i,j,p)))-n3(i,j,p)*(1.-n2(i,j,p))*dt_over_tau32(p) & 
			  +EdP2ij(p)+pumping_rate(p)*(1-n3(i,j,p))
			n2(i,j,p)=n2(i,j,p)*(1.-dt_over_tau21(p)*(1.-n1(i,j,p)))+n3(i,j,p)*(1.-n2(i,j,p))*dt_over_tau32(p) & 
			  +EdP1ij(p)
			n1(i,j,p)=n1(i,j,p)*(1.-dt_over_tau10(p)*(1.-n0(i,j,p)))+n2(i,j,p)*(1.-n1(i,j,p))*dt_over_tau21(p) & 
			  -EdP1ij(p)
			n0(i,j,p)=ne(p)-n1(i,j,p)-n2(i,j,p)-n3(i,j,p)-pumping_rate(p)*(1-n3(i,j,p))
				endif
			else      !A.p
				if(ipauli .eq. 2)then
			AP2ij(p)=C_n2_1(p)*Avptz(i,j)*Pz_2(i,j,p)
			AP1ij(p)=C_n1_1(p)*Avptz(i,j)*Pz_1(i,j,p)
			n3(i,j,p)=n3(i,j,p)*(1.-dt_over_tau30(p)*(1.-n0(i,j,p)))-n3(i,j,p)*(1.-n2(i,j,p))*dt_over_tau32(p) &
				-AP2ij(p)+pumping_rate(p)*(1-n3(i,j,p))
			n2(i,j,p)=n2(i,j,p)*(1.-dt_over_tau21(p)*(1.-n1(i,j,p)))+n3(i,j,p)*(1.-n2(i,j,p))*dt_over_tau32(p) &
				-AP1ij(p)
			n1(i,j,p)=n1(i,j,p)*(1.-dt_over_tau10(p)*(1.-n0(i,j,p)))+n2(i,j,p)*(1.-n1(i,j,p))*dt_over_tau21(p) & 
				+AP1ij(p)
			n0(i,j,p)=ne(p)-n1(i,j,p)-n2(i,j,p)-n3(i,j,p)-pumping_rate(p)*(1-n3(i,j,p))
				endif
			endif
		endif	!TM-mode
	!end medium
	endif
101 continue
 

	!-----------------------------------------------------------------------------------
	!		Ez (Hz) UPDATE in the UPML grid
	!-----------------------------------------------------------------------------------

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! front for '0' process in mpi
		do 1010 j=1,jebc
		do 1010 i=1,iefbc
	  	  
			!TM-mode
			!if(TM .gt. 0)then
			DSbcf=DEZf(i,j)
			DEZf(i,j)=axef(i)*DEZf(i,j)			  &
				+(HYf(i+1,j)-HYf(i,j)-HXf(i,j+1)+HXf(i,j))*dt_over_dx/bxef(i)
			EZf(i,j)=ayef(j)*EZf(i,j)			     &
				+(DEZf(i,j)-DSbcf)/epsrf(i,j)/byef(j)		   ! epsrfb(i,j)
			!endif
			!TE-mode
			If(TE .gt. 0)then
			BHZSbcf=BHZf(i,j)
			BHZf(i,j)=axef(i)*BHZf(i,j)			  &
				-(EYf(i+1,j)-EYf(i,j)-EXf(i,j+1)+EXf(i,j))*dt_over_dx/bxef(i)
			Hzf(i,j)=ayef(j)*Hzf(i,j)			     &
				+(BHZf(i,j)-BHZSbcf)*muer/byef(j)		   ! epsrfb(i,j)
			endif	

1010   continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! back region - for 'numprocs-1' process in mpi
		do 1011 j=1,jebc
		do 1011 i=1,iefbc
	  	  
			!TM-mode
			!if(TM .gt. 0)then
			DSbcb=Dezb(i,j)
			DEZb(i,j)=axeb(i)*DEZb(i,j)			  &
				+(HYb(i+1,j)-HYb(i,j)-HXb(i,j+1)+HXb(i,j))*dt_over_dx/bxeb(i)
			EZb(i,j)=ayeb(j)*EZb(i,j)			     &
				+(DEZb(i,j)-DSbcb)/epsrb(i,j)/byeb(j)		   ! epsrfb(i,j)
			!endif
			!TE-mode
			If(TE .gt. 0)then
			BHZSbcb=BHzb(i,j)
			BHZb(i,j)=axeb(i)*BHZb(i,j)			  &
				-(EYb(i+1,j)-EYb(i,j)-EXb(i,j+1)+EXb(i,j))*dt_over_dx/bxeb(i)
			Hzb(i,j)=ayeb(j)*Hzb(i,j)			     &
				+(BHZb(i,j)-BHZSbcb)*muer/byeb(j)		   ! epsrfb(i,j)
			endif	

1011   continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

! left and right regions are segmented according to je_np array
! same for mpi or without it

		! left region
		do 1012 j=j1c,j2c
		do 1012 i=1,iebc

			!TM-mode
			!if(TM .gt. 0)then
			Dsbcl=Dezbcl(i,j)
			DEZbcl(i,j)=axel(i)*DEZbcl(i,j)			  &
				+(HYbcl(i+1,j)-HYbcl(i,j)-HXbcl(i,j+1)+HXbcl(i,j))*dt_over_dx/bxel(i)
			EZbcl(i,j)=ayel(j)*EZbcl(i,j)			     &
				+(DEZbcl(i,j)-DSbcl)/epsrbcl(i,j)/byel(j)	
			!endif
			!TE-mode
			If(TE .gt. 0)then
			BHZSbcl=BHzbcl(i,j)
			BHZbcl(i,j)=axel(i)*BHZbcl(i,j)			  &
				-(EYbcl(i+1,j)-EYbcl(i,j)-EXbcl(i,j+1)+EXbcl(i,j))*dt_over_dx/bxel(i)
			Hzbcl(i,j)=ayel(j)*Hzbcl(i,j)			     &
				+(BHZbcl(i,j)-BHZSbcl)*muer/byel(j)		   ! epsrfb(i,j)
			endif	
			   
1012   continue

		! right region
		do 1013 j=j1c,j2c
		do 1013 i=1,iebc

			!TM-mode
			!if(TM .gt. 0)then
			Dsbcr=Dezbcr(i,j)
			DEZbcr(i,j)=axer(i)*DEZbcr(i,j)			  &
				+(HYbcr(i+1,j)-HYbcr(i,j)-HXbcr(i,j+1)+HXbcr(i,j))*dt_over_dx/bxer(i)
			EZbcr(i,j)=ayer(j)*EZbcr(i,j)			     &
				+(DEZbcr(i,j)-DSbcr)/epsrbcr(i,j)/byer(j)				  
			!endif
			!TE-mode
			If(TE .gt. 0)then
			BHZSbcr=BHZbcr(i,j)
			BHZbcr(i,j)=axer(i)*BHZbcr(i,j)			  &
				-(EYbcr(i+1,j)-EYbcr(i,j)-EXbcr(i,j+1)+EXbcr(i,j))*dt_over_dx/bxer(i)
			Hzbcr(i,j)=ayer(j)*Hzbcr(i,j)			     &
				+(BHZbcr(i,j)-BHZSbcr)*muer/byer(j)		   ! epsrfb(i,j)
			endif	
			   
1013   continue

	!-------------------------------------------------------------------
	!		Source for Ez or Hz
	!-------------------------------------------------------------------
	do 1014 m=1,num_wg
		do 1015 k=1,nfield(m)
			if ((inputfield(m,k) .eq. 1) .and.(k .eq. 1) ) then	! First Pulse
				tdelay=tdelay0+1.5*pwidth(m,k)+tdel(m,k)
				beamshape=exp(-(tt-tdelay)**2/(pwidth(m,k)/2)**2)
			elseif (inputfield(m,k) .eq. 1) then	!second pulse
				tdelay=tdelay0+1.5*pwidth(m,k)+tdel(m,1)+tdel(m,2)
				beamshape=exp(-(tt-tdelay)**2/(pwidth(m,k)/2)**2)
			elseif(tt .lt. tdel(m,k))then	  !CW
				beamshape=0.
			else
				beamshape=1.
			endif

		!-----------------------------------------------------------
		!		Total Field/Scattered field for TM-mode (Ez)
		!-----------------------------------------------------------
			
			if ((pol(m,k) .eq. 1) .and. (wgdirection(m) .eq. 1).and.(ilaunch(m,k) .gt.0)) then
				do 1020 j=js1,js2
					Hy_inc=-beamshape*(kprop(m,k)/w_f(m,k)/muz)*E_input(m,k)*		&
							sin(w_f(m,k)*(tt-dt/2.) + kprop(m,k)*dx/2.)*rmode(m,k,j)
					Ez(ilaunch(m,k),j)=Ez(ilaunch(m,k),j) - cb*geometry(ilaunch(m,k),j)*Hy_inc
1020			continue

			elseif ((pol(m,k) .eq. 1) .and. (wgdirection(m) .eq. 1).and.(ilaunch(m,k) .lt.0)) then
				do 1021 j=js1,js2
					Hy_inc=-beamshape*(kprop(m,k)/w_f(m,k)/muz)*E_input(m,k)*		&
							sin(w_f(m,k)*(tt-dt/2.) - kprop(m,k)*dx/2.)*rmode(m,k,j)
					ilaunchmk = ie + ilaunch(m,k)
					Ez(ilaunchmk,j)=Ez(ilaunchmk,j) + cb*geometry(ilaunchmk,j)*Hy_inc
1021			continue

			elseif((pol(m,k) .eq. 1) .and. (wgdirection(m) .eq. 2).and.(ilaunch(m,k).gt.0))then
			! front only

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

				do 1022 i=is1,is2
					Hx_inc=beamshape*(kprop(m,k)/w_f(m,k)/muz)*E_input(m,k)*		&
							sin(w_f(m,k)*(tt-dt/2.) + kprop(m,k)*dx/2.)*rmode(m,k,i)
					Ez(i,ilaunch(m,k))=Ez(i,ilaunch(m,k)) + cb*geometry(i,ilaunch(m,k))*Hx_inc
1022			continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

			elseif((pol(m,k) .eq. 1) .and. (wgdirection(m) .eq. 2).and.(ilaunch(m,k).lt.0))then
			! back only

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

				do 1023 i=is1,is2
					Hx_inc=beamshape*(kprop(m,k)/w_f(m,k)/muz)*E_input(m,k)*		&
							sin(w_f(m,k)*(tt-dt/2.) - kprop(m,k)*dx/2.)*rmode(m,k,i)
					ilaunchmk = je + ilaunch(m,k)
					Ez(i,ilaunchmk)=Ez(i,ilaunchmk) - cb*geometry(i,ilaunchmk)*Hx_inc
1023			continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

		!-----------------------------------------------------------
		!		Total Field/Scattered field for TE-mode (Hz)
		!-----------------------------------------------------------

			elseif((pol(m,k) .eq. 0) .and. (wgdirection(m) .eq. 1) .and. (ilaunch(m,k) .gt. 0))then
				do 1024 j=js1,js2
					Ey_inc=cc*cc*geometry(ilaunch(m,k),j)*kprop(m,k)*kprop(m,k)/w_f(m,k)/w_f(m,k)*beamshape*E_input(m,k)*		&
							sin(w_f(m,k)*(tt-dt/2.)+kprop(m,k)*dx/2.)*rmode(m,k,j)

!! old version		Ey_inc=beamshape*E_input(m,k)*sin(w_f(m,k)*(tt-dt/2.)+kprop(m,k)*dx/2.)*rmode(m,k,j)
					Hz(ilaunch(m,k),j)=Hz(ilaunch(m,k),j)+db*Ey_inc
1024			continue

			elseif((pol(m,k) .eq. 0) .and. (wgdirection(m) .eq. 1) .and. (ilaunch(m,k) .lt. 0))then
				do 1025 j=js1,js2
					ilaunchmk = ie + ilaunch(m,k)

					Ey_inc=cc*cc*geometry(ilaunchmk,j)*kprop(m,k)*kprop(m,k)/w_f(m,k)/w_f(m,k)*beamshape*E_input(m,k)*		&
							sin(w_f(m,k)*(tt-dt/2.)-kprop(m,k)*dx/2.)*rmode(m,k,j)
!! old version		Ey_inc=beamshape*E_input(m,k)*sin(w_f(m,k)*(tt-dt/2.)-kprop(m,k)*dx/2.)*rmode(m,k,j)

					Hz(ilaunchmk,j)=Hz(ilaunchmk,j)-db*Ey_inc
1025			continue

 			elseif((pol(m,k) .eq. 0) .and. (wgdirection(m) .eq. 2) .and. (ilaunch(m,k) .gt. 0))then
			! front only
!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 
				do 1026 i=is1,is2

					Ex_inc=cc*cc*geometry(ilaunch(m,k),j)*kprop(m,k)*kprop(m,k)/w_f(m,k)/w_f(m,k)*beamshape*E_input(m,k)*		&
							sin(w_f(m,k)*(tt-dt/2.)+kprop(m,k)*dx/2.)*rmode(m,k,j)
!! old version		Ex_inc=beamshape*E_input(m,k)*sin(w_f(m,k)*(tt-dt/2.)+kprop(m,k)*dx/2.)*rmode(m,k,i)
					if(i.eq. 250) then
						Ex_inc250=Ex_inc
					endif
					Hz(i,ilaunch(m,k))=Hz(i,ilaunch(m,k))-db*Ex_inc
1026			continue
!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

			else ! back only

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

				do 1027 i=is1,is2
					ilaunchmk = je + ilaunch(m,k)

					Ex_inc=cc*cc*geometry(ilaunchmk,j)*kprop(m,k)*kprop(m,k)/w_f(m,k)/w_f(m,k)*beamshape*E_input(m,k)*		&
							sin(w_f(m,k)*(tt-dt/2.)+kprop(m,k)*dx/2.)*rmode(m,k,j)
!!					Ex_inc=beamshape*E_input(m,k)*sin(w_f(m,k)*(tt-dt/2.)-kprop(m,k)*dx/2.)*rmode(m,k,i)
					Hz(i,ilaunchmk)=Hz(i,ilaunchmk)+db*Ex_inc
1027			continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

			endif

1015	continue
1014	continue
!write(*,*) "Ez source update"  

!DEC$ IF DEFINED (USE_MPI_CH) 
		! (m) main grid: exchange lines (j=..) between processes before getting Hx, etc. ..
		! send bottom line of top process to bottom process 
		! no data exchange needed for calculating Hy, Ey ..

    call MPI_BARRIER( comm_a, ierr ) ! synchronize entrance;
!if( pid .eq. 0 ) then 
!Ez(1,1)=10
!write(*,*) Ez(1:10,1)
!endif
		!if(TM .gt. 0)then
			!(m) for Hx, send Ez 0 -> 1, 1 -> 2, etc., 
			call exchng3( Ez, ie, j1_al, j2_al, je_np, pid_a, nbrtop_a, &
										nbrbottom_a, numprocs, n, comm_a )
!write(*,*) 'exchange' 
!if( pid .eq. 0 ) then
!write(*,*) Ez(1:10,1)
!endif
			! (upml)	EZbcl, EZbcr 0 -> 1, 1 -> 2, etc.
			call exchng3( Ezbcl, iebc, j1_al, j2_al, je_np, pid_a, nbrtop_a, &
										nbrbottom_a, numprocs, n, comm_a )
			call exchng3( Ezbcr, iebc, j1_al, j2_al, je_np, pid_a, nbrtop_a, &
										nbrbottom_a, numprocs, n, comm_a )

		!endif

		if(TE .gt. 0)then
			!(m) for Ex, send Hz 0 -> 1, 1 -> 2, etc.,
			call exchng3( Hz, ie, j1_al, j2_al, je_np, pid_a, nbrtop_a, &
										nbrbottom_a, numprocs, n, comm_a )
			! (upml)	EZbcl, EZbcr 0 -> 1, 1 -> 2, etc.
			call exchng3( Hzbcl, iebc, j1_al, j2_al, je_np, pid_a, nbrtop_a, &
										nbrbottom_a, numprocs, n, comm_a )
			call exchng3( Hzbcr, iebc, j1_al, j2_al, je_np, pid_a, nbrtop_a, &
										nbrbottom_a, numprocs, n, comm_a )
		endif

!DEC$ ENDIF ! (USE_MPI_CH) 
!write(*,*) 'Ez exchanged' 

	!------------------------------------------------------------------------ 
	!			HX (Ex) and HY (Ey) UPDATE in the MAIN grid
	!------------------------------------------------------------------------												  
	do 4501 j=j1b,j2c
	do 4501 i=i1c,i2c

		!TM-mode
		!if(TM .gt. 0)then
		Hx(i,j)=Hx(i,j)-db*(Ez(i,j)-Ez(i,j-1))
		!endif

		!TE-mode
		if(TE .gt. 0)then
			geometryij=max(geometry(i,j),geometry(i,j-1))
			
			if (gain .gt. 0) then
			gainij=min(gaingeo(i,j),gaingeo(i,j-1))
			else
			gainij=0
			endif
			
			if  (gainij.eq.0)	then
				if(metal .eq. 0)then
				Ex(i,j)=Ex(i,j)+cb*geometryij*(Hz(i,j)-Hz(i,j-1))
				else
					if (geometry(i,j) .ge. geometry(i,j-1))then
					sigmageoij=sigmageo(i,j)
					else
					sigmageoij=sigmageo(i,j-1)
					endif
				caij=(1-sigmageoij*dt/2./epsz*geometryij)/(1+sigmageoij*dt/2./epsz*geometryij)
				cbij=cb/(1+sigmageoij*dt/2./epsz*geometryij)
				Ex(i,j)=caij*Ex(i,j)+cbij*geometryij*(Hz(i,j)-Hz(i,j-1))
				endif
			else

			!Polarizaion
			Px_old2_1(:)=px_old1_1(i,j,:)
			px_old1_1(i,j,:)=px_1(i,j,:)

			Px_old2_2(:)=px_old1_2(i,j,:)
			px_old1_2(i,j,:)=px_2(i,j,:)
	
			n0ij(:)=(n0(i,j,:)+n0(i,j-1,:))/2.
			n1ij(:)=(n1(i,j,:)+n1(i,j-1,:))/2.
			n2ij(:)=(n2(i,j,:)+n2(i,j-1,:))/2.
			n3ij(:)=(n3(i,j,:)+n3(i,j-1,:))/2.

			p=gainij
				if(interaction .eq.1) then
			px_1(i,j,p) = kapa1(p)*(n1ij(p)-n2ij(p))*Ex(i,j) &
			 +pa2_1(p)*px_old1_1(i,j,p)  &
			 +pa3_1(p)*px_old2_1(p)
			px_2(i,j,p) = kapa2(p)*(n0ij(p)-n3ij(p))*Ex(i,j) &
			 +pa2_2(p)*px_old1_2(i,j,p)  &
			 +pa3_2(p)*px_old2_2(p)
				else
			px_1(i,j,p) = kapa1(p)*(n1ij(p)-n2ij(p))*Ex(i,j) &
			 +(pa2_1(p)-rabico1(p)*(Avptx(i,j))**2)*px_old1_1(i,j,p) &
			 +pa3_1(p)*px_old2_1(p)
			px_2(i,j,p) = kapa2(p)*(n0ij(p)-n3ij(p))*Ex(i,j) &
			 +(pa2_2(p)-rabico2(p)*(Avptx(i,j))**2)*px_old1_2(i,j,p) &
			 +pa3_2(p)*px_old2_2(p)
				endif
	 	
			!Electric field  
			ex_old(i,j)=ex(i,j)
			ex(i,j)=ex(i,j)+cb*geometryij*(Hz(i,j)-Hz(i,j-1))
			ex(i,j)=ex(i,j)-1/epsz*geometryij*Ndensity(p)*(px_1(i,j,p)+px_2(i,j,p)   &
			-px_old1_1(i,j,p)-px_old1_2(i,j,p))
		
			!Vector potential
				if(interaction .eq. 2)then
			Avptx_old(i,j)=Avptx(i,j)
			Avptx(i,j)=Avptx(i,j)-dt/2.*(Ex(i,j)+Ex_old(i,j))
				endif
			endif
		endif	!TE-mode

4501 continue
       
	do 4502 j=j1c,j2c
	do 4502 i=i1b,i2c

		!TM-mode
		!if(TM .gt. 0)then
		Hy(i,j)=Hy(i,j)+db*(Ez(i,j)-Ez(i-1,j))
		!endif

		!TE-mode
		if(TE .gt. 0)then
		geometryij=min(geometry(i,j),geometry(i-1,j))
		  if (gain .gt. 0) then
		  gainij=max(gaingeo(i,j),gaingeo(i-1,j))
		  else
		  gainij=0
		  endif
			if (gainij.eq.0)then
				if(metal .eq. 0)then
 				Ey(i,j) =Ey(i,j)+cb*geometryij*(Hz(i-1,j)-Hz(i,j))
				else
					if (geometry(i,j) .ge. geometry(i-1,j))then
					sigmageoij=sigmageo(i-1,j)
					else
					sigmageoij=sigmageo(i,j)
					endif
				caij=(1-sigmageoij*dt/2./epsz*geometryij)/(1+sigmageoij*dt/2./epsz*geometryij)
				cbij=cb/(1+sigmageoij*dt/2./epsz*geometryij)
				Ey(i,j) =caij*Ey(i,j)+cbij*geometryij*(Hz(i-1,j)-Hz(i,j))
				endif
			else
			!Polarizaion
			Py_old2_1(:) =py_old1_1(i,j,:)
			py_old1_1(i,j,:)=py_1(i,j,:)

			Py_old2_2(:)  =py_old1_2(i,j,:)
			py_old1_2(i,j,:)=py_2(i,j,:)

			n0ij(:)=(n0(i,j,:)+n0(i-1,j,:))/2.
			n1ij(:)=(n1(i,j,:)+n1(i-1,j,:))/2.
			n2ij(:)=(n2(i,j,:)+n2(i-1,j,:))/2.
			n3ij(:)=(n3(i,j,:)+n3(i-1,j,:))/2.

			p=gainij
				if(interaction .eq.1) then
			py_1(i,j,p) = kapa1(p)*(n1ij(p)-n2ij(p))*Ey(i,j) &
			 +pa2_1(p)*py_old1_1(i,j,p)  &
			 +pa3_1(p)*py_old2_1(p)
			py_2(i,j,p) =kapa2(p)*(n0ij(p)-n3ij(p))*Ey(i,j) &
			 +pa2_2(p)*py_old1_2(i,j,p)  &
			 +pa3_2(p)*py_old2_2(p)
				else
			py_1(i,j,p) = kapa1(p)*(n1ij(p)-n2ij(p))*Ey(i,j) &
			 +(pa2_1(p)-rabico1(p)*(Avpty(i,j))**2)*py_old1_1(i,j,P) &
			 +pa3_1(p)*py_old2_1(p)
			py_2(i,j,p) = kapa2(p)*(n0ij(p)-n3ij(p))*Ey(i,j) &
			 +(pa2_2(p)-rabico2(p)*(Avpty(i,j))**2)*py_old1_2(i,j,p) &
			 +pa3_2(p)*py_old2_2(p)
				endif
			!Electric field  
			ey_old(i,j)=ey(i,j)
			ey(i,j) =ey(i,j)+cb*geometryij*(Hz(i-1,j)-Hz(i,j))
			ey(i,j) =ey(i,j)-1/epsz*geometryij*Ndensity(p)*(py_1(i,j,p)+py_2(i,j,p)   &
				-py_old1_1(i,j,p)-py_old1_2(i,j,p))
			!Vector potential
			Avpty_old(i,j)=Avpty(i,j)
			Avpty(i,j)=Avpty(i,j)-dt/2.*(Ey(i,j)+Ey_old(i,j))
			endif
		endif	!TE-mode
4502 continue
	 !-----------------------------------------------------------------------
	 !			Hx (Ex) and Hy (Ey) UPDATE in the UPML grid
	 !----------------------------------------------------------------------- 

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! front Hx, Hy update
		do 4503 j=2,jebc
		do 4503 i=1,iefbc

			!TM-mode
			!if(TM .gt. 0)then
			BHXSbcf=BHXf(i,j)
			BHXf(i,j)=ayhf(j)*BHXf(i,j)-(EZf(i,j)-EZf(i,j-1))*dt_over_dx/byhf(j)
			HXf(i,j)=HXf(i,j)+(bxef(i)*BHXf(i,j)-cxef(i)*BHXSbcf)*muer			   
			!endif
			!TE-mode
			if(TE .gt. 0)then
			DEXSbcf=DEXf(i,j)
			DEXf(i,j)=ayhf(j)*DEXf(i,j)+(HZf(i,j)-HZf(i,j-1))*dt_over_dx/byhf(j)
			EXf(i,j)=EXf(i,j)+(bxef(i)*DEXf(i,j)-cxef(i)*DEXSbcf)/epsrf(i,j)
			endif
			  
4503 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		elseif( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! back Hx, Hy update
		do 4513 j=2,jebc
		do 4513 i=1,iefbc

			!TM-mode
			!if(TM .gt. 0)then
			BHXSbcb=BHXb(i,j)
			BHXb(i,j)=ayhb(j)*BHXb(i,j)-(EZb(i,j)-EZb(i,j-1))*dt_over_dx/byhb(j)
			HXb(i,j)=HXb(i,j)+(bxeb(i)*BHXb(i,j)-cxeb(i)*BHXSbcb)*muer
			!endif
			!TE-mode
			if(TE .gt. 0)then
			DEXSbcb=DEXb(i,j)
			DEXb(i,j)=ayhb(j)*DEXb(i,j)+(HZb(i,j)-HZb(i,j-1))*dt_over_dx/byhb(j)
			EXb(i,j)=EXb(i,j)+(bxeb(i)*DEXb(i,j)-cxeb(i)*DEXSbcb)/epsrb(i,j)
			endif  

4513 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! front
		do 4504 j=1,jebc 
		do 4504 i=2,iefbc

			!TM-mode
			!if(TM .gt. 0)then
			BHYSbcf=BHYf(i,j)
			BHYf(i,j)=BHYf(i,j)+(EZf(i,j)-EZf(i-1,j))*dt_over_dx
			HYf(i,j)=axhf(i)*HYf(i,j)+(byef(j)*BHYf(i,j)-cyef(j)*BHYSbcf)*muer/bxhf(i)
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEYSbcf=DEYf(i,j)
			epsrfij=min(epsrf(i,j),epsrf(i-1,j))
			DEYf(i,j)=DEYf(i,j)-(Hzf(i,j)-HZf(i-1,j))*dt_over_dx
			EYf(i,j)=axhf(i)*EYf(i,j)+(byef(j)*DEYf(i,j)-cyef(j)*DEYSbcf)/epsrfij/bxhf(i)
			endif

4504 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		elseif( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! back
		do 4514 j=1,jebc 
		do 4514 i=2,iefbc

			!TM-mode
			!if(TM .gt. 0)then
			BHYSbcb=BHYb(i,j)
			BHYb(i,j)=BHYb(i,j)+(EZb(i,j)-EZb(i-1,j))*dt_over_dx
			HYb(i,j)=axhb(i)*HYb(i,j)+(byeb(j)*BHYb(i,j)-cyeb(j)*BHYSbcb)*muer/bxhb(i)
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEYSbcb=DEYb(i,j)
			epsrbij=min(epsrb(i,j),epsrb(i-1,j))
			DEYb(i,j)=DEYb(i,j)-(HZb(i,j)-HZb(i-1,j))*dt_over_dx
			EYb(i,j)=axhb(i)*EYb(i,j)+(byeb(j)*DEYb(i,j)-cyeb(j)*DEYSbcb)/epsrbij/bxhb(i)
			endif

4514 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

		! left Hx (Ex), Hy (Ey) update
		do 4505 j=j1b,j2c
		do 4505 i=1,iebc

			!TM-mode
			!if(TM .gt. 0)then
			BHXSbcl=BHXbcl(i,j)
			BHXbcl(i,j)=ayhl(j)*BHXbcl(i,j)-(EZbcl(i,j)-EZbcl(i,j-1))*dt_over_dx/byhl(j)    
			HXbcl(i,j)=HXbcl(i,j)+(bxel(i)*BHXbcl(i,j)-cxel(i)*BHXSbcl)*muer			   
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEXSbcl=DEXbcl(i,j)
			DEXbcl(i,j)=ayhl(j)*DEXbcl(i,j)+(HZbcl(i,j)-HZbcl(i,j-1))*dt_over_dx/byhl(j)    
			EXbcl(i,j)=EXbcl(i,j)+(bxel(i)*DEXbcl(i,j)-cxel(i)*DEXSbcl)/epsrbcl(i,j)			   
			endif  

4505 continue

    do 4506 j=j1c,j2c
		do 4506 i=i1b,iebc

			!TM-mode
			!if(TM .gt. 0)then
			BHYSbcl=BHYbcl(i,j)
			BHYbcl(i,j)=BHYbcl(i,j)+(EZbcl(i,j)-EZbcl(i-1,j))*dt_over_dx
			HYbcl(i,j)=axhl(i)*HYbcl(i,j)+(byel(j)*BHYbcl(i,j)-cyel(j)*BHYSbcl)*muer/bxhl(i)
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEYSbcl=DEYbcl(i,j)
			DEYbcl(i,j)=DEYbcl(i,j)-(HZbcl(i,j)-HZbcl(i-1,j))*dt_over_dx
			EYbcl(i,j)=axhl(i)*EYbcl(i,j)+(byel(j)*DEYbcl(i,j)- &
						cyel(j)*DEYSbcl)/epsrbcl(i,j)/bxhl(i)
			endif  
4506 continue
	  
		! right Hx (Ex), Hy (Ey) update
		do 4515 j=j1b,j2c 
		do 4515 i=1,iebc

			!TM-mode
			!if(TM .gt. 0)then
			BHXSbcr=BHXbcr(i,j)
			BHXbcr(i,j)=ayhr(j)*BHXbcr(i,j)-(EZbcr(i,j)-EZbcr(i,j-1))*dt_over_dx/byhr(j)
			HXbcr(i,j)=HXbcr(i,j)+(bxer(i)*BHXbcr(i,j)-cxer(i)*BHXSbcr)*muer
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEXSbcr=DEXbcr(i,j)
			DEXbcr(i,j)=ayhr(j)*DEXbcr(i,j)+(HZbcr(i,j)-HZbcr(i,j-1))*dt_over_dx/byhr(j)
			EXbcr(i,j)=EXbcr(i,j)+(bxer(i)*DEXbcr(i,j)-cxer(i)*DEXSbcr)/epsrbcr(i,j)
			endif  

4515 continue
		
    do 4516 j=j1c,j2c
		do 4516 i=i1b,iebc

			!TM-mode
			!if(TM .gt. 0)then
			BHYSbcr=BHYbcr(i,j)
			BHYbcr(i,j)=BHYbcr(i,j)+(EZbcr(i,j)-EZbcr(i-1,j))*dt_over_dx
			HYbcr(i,j)=axhr(i)*HYbcr(i,j)+(byer(j)*BHYbcr(i,j)-cyer(j)*BHYSbcr)*muer/bxhr(i)
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEYSbcr=DEYbcr(i,j)
			DEYbcr(i,j)=DEYbcr(i,j)-(HZbcr(i,j)-HZbcr(i-1,j))*dt_over_dx
			EYbcr(i,j)=axhr(i)*EYbcr(i,j)+ &
						(byer(j)*DEYbcr(i,j)-cyer(j)*DEYSbcr)/epsrbcr(i,j)/bxhr(i)
			endif  
4516 continue


	!-------------------------------------------------------------------------------
	!		Hx (Ex),Hy (Ey) UPDATE in the Grid/PML interface 
	!----------------------------------------------------------------------------- 
!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! front
		do 4601 i=i1c,i2c
			i1=iebc+i

			!TM-mode
			!if(TM .gt. 0)then
			BHXSbcf=BHXf(i1,jbbc)
			BHXf(i1,jbbc)=ayhf(jbbc)*BHXf(i1,jbbc)-(EZ(i,1)-EZf(i1,jebc))*dt_over_dx/byhf(jbbc)
			HXf(i1,jbbc)=HXf(i1,jbbc)+(bxef(i1)*BHXf(i1,jbbc)-cxef(i1)*BHXSbcf)*muer			   
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEXSbcf=DEXf(i1,jbbc)
			DEXf(i1,jbbc)=ayhf(jbbc)*DEXf(i1,jbbc)+(HZ(i,1)-HZf(i1,jebc))*dt_over_dx/byhf(jbbc)
			EXf(i1,jbbc)=EXf(i1,jbbc)+(bxef(i1)*DEXf(i1,jbbc)-cxef(i1)*DEXSbcf)/epsrf(i1,jebc)			   
			endif  
4601 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		elseif( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! back
		do 4611 i=i1c,i2c
			i1=iebc+i

			!TM-mode
			!if(TM .gt. 0)then
			BHXSbcb=BHXb(i1,1)
			BHXb(i1,1)=ayhb(1)*BHXb(i1,1)-(EZb(i1,1)-EZ(i,je))*dt_over_dx/byhb(1)
			HXb(i1,1)=HXb(i1,1)+(bxeb(i1)*BHXb(i1,1)-cxeb(i1)*BHXSbcb)*muer
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEXSbcb=DEXb(i1,1)
			DEXb(i1,1)=ayhb(1)*DEXb(i1,1)+(HZb(i1,1)-HZ(i,je))*dt_over_dx/byhb(1)
			EXb(i1,1)=EXb(i1,1)+(bxeb(i1)*DEXb(i1,1)-cxeb(i1)*DEXSbcb)/epsrb(i1,1)
			endif  
4611 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

		! left
		do 4602 j=j1c,j2c

			!TM-mode
			!if(TM .gt. 0)then
			BHYSbcl=BHYbcl(ibbc,j)
			BHYbcl(ibbc,j)=BHYbcl(ibbc,j)+(EZ(1,j)-EZbcl(iebc,j))*dt_over_dx
			HYbcl(ibbc,j)=axhl(ibbc)*HYbcl(ibbc,j)+        &
							(byel(j)*BHYbcl(ibbc,j)-cyel(j)*BHYSbcl)*muer/bxhl(ibbc)
			!endif
			!TE-mode
			if(TE .gt. 0)then
			DEYSbcl=DEYbcl(ibbc,j)
			DEYbcl(ibbc,j)=DEYbcl(ibbc,j)-(HZ(1,j)-HZbcl(iebc,j))*dt_over_dx
			EYbcl(ibbc,j)=axhl(ibbc)*EYbcl(ibbc,j)+        &
					(byel(j)*DEYbcl(ibbc,j)-cyel(j)*DEYSbcl)/epsrbcl(iebc,j)/bxhl(ibbc)
			endif

4602 continue

		! right
		do 4612 j=j1c,j2c

			!TM-mode
			!if(TM .gt. 0)then
			BHYSbcr=BHYbcr(1,j)
			BHYbcr(1,j)=BHYbcr(1,j)+(EZbcr(1,j)-EZ(ie,j))*dt_over_dx
			HYbcr(1,j)=axhr(1)*HYbcr(1,j)+        &
										(byer(j)*BHYbcr(1,j)-cyer(j)*BHYSbcr)*muer/bxhr(1)
			!endif
			!TE-mode
			if(TE .gt. 0)then
			DEYSbcr=DEYbcr(1,j)
			DEYbcr(1,j)=DEYbcr(1,j)-(HZbcr(1,j)-HZ(ie,j))*dt_over_dx
			EYbcr(1,j)=axhr(1)*EYbcr(1,j)+        &
									(byer(j)*DEYbcr(1,j)-cyer(j)*DEYSbcr)/epsrbcr(1,j)/bxhr(1)
			endif  
4612 continue


	!-------------------------------------------------------------------
	!		Hx (Ex),Hy (Ey) UPDATE in the UPML/UPML interface 
	!-------------------------------------------------------------------

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! front, (left, right) UPML/UPML interface
		do 4701 i1=1,iebc
			i2=iebc+ie+i1

			!TM-mode
			!if(TM .gt. 0)then
			BHXSbcf1=BHXf(i1,jbbc)
			BHXSbcf2=BHXf(i2,jbbc)
			BHXf(i1,jbbc)=ayhf(jbbc)*BHXf(i1,jbbc)-				&
					(EZbcl(i1,1)-EZf(i1,jebc))*dt_over_dx/byhf(jbbc)
			BHXf(i2,jbbc)=ayhf(jbbc)*BHXf(i2,jbbc)-				&
					(EZbcr(i1,1)-EZf(i2,jebc))*dt_over_dx/byhf(jbbc)
			HXf(i1,jbbc)=HXf(i1,jbbc)+(bxef(i1)*BHXf(i1,jbbc)-cxef(i1)*BHXSbcf1)*muer
			HXf(i2,jbbc)=HXf(i2,jbbc)+(bxef(i2)*BHXf(i2,jbbc)-cxef(i2)*BHXSbcf2)*muer			   
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEXSbcf1=DEXf(i1,jbbc)
			DEXSbcf2=DEXf(i2,jbbc)
			DEXf(i1,jbbc)=ayhf(jbbc)*DEXf(i1,jbbc)+	&
															(Hzbcl(i1,1)-Hzf(i1,jebc))*dt_over_dx/byhf(jbbc)
			DEXf(i2,jbbc)=ayhf(jbbc)*DEXf(i2,jbbc)+	&
															(Hzbcr(i1,1)-Hzf(i2,jebc))*dt_over_dx/byhf(jbbc)
			EXf(i1,jbbc)=EXf(i1,jbbc)+	&
		            (bxef(i1)*DEXf(i1,jbbc)-cxef(i1)*DEXSbcf1)/epsz
			EXf(i2,jbbc)=EXf(i2,jbbc)+	&
		            (bxef(i2)*DEXf(i2,jbbc)-cxef(i2)*DEXSbcf2)/epsz			   
			endif  !TE-mode

4701 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		elseif( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! back, (left, right) UPML/UPML interface 
		do 4711 i1=1,iebc
			i2=iebc+ie+i1

			!TM-mode
			!if(TM .gt. 0)then
			BHXSbcb1=BHXb(i1,1)
	   	BHXSbcb2=BHXb(i2,1)
			BHXb(i1,1)=ayhb(1)*BHXb(i1,1)-(EZb(i1,1)-EZbcl(i1,je))*dt_over_dx/byhb(1)
			BHXb(i2,1)=ayhb(1)*BHXb(i2,1)-(EZb(i2,1)-EZbcr(i1,je))*dt_over_dx/byhb(1)
			HXb(i1,1)=HXb(i1,1)+(bxeb(i1)*BHXb(i1,1)-cxeb(i1)*BHXSbcb1)*muer
			HXb(i2,1)=HXb(i2,1)+(bxeb(i2)*BHXb(i2,1)-cxeb(i2)*BHXSbcb2)*muer
			!endif

			!TE-mode
			if(TE .gt. 0)then
			DEXSbcb1=DEXb(i1,1)
   		DEXSbcb2=DEXb(i2,1)
			DEXb(i1,1)=ayhb(1)*DEXb(i1,1)+(Hzb(i1,1)-Hzbcl(i1,je))*dt_over_dx/byhb(1)
			DEXb(i2,1)=ayhb(1)*DEXb(i2,1)+(Hzb(i2,1)-Hzbcr(i1,je))*dt_over_dx/byhb(1)
			EXb(i1,1)=EXb(i1,1)+(bxeb(i1)*DEXb(i1,1)-cxeb(i1)*DEXSbcb1)/epsz
			EXb(i2,1)=EXb(i2,1)+(bxeb(i2)*DEXb(i2,1)-cxeb(i2)*DEXSbcb2)/epsz
			endif  !TE-mode

4711 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! put front values to left & right at PML/PMl interface 	
		do 4702 i1=1,iebc

			!TM-mode
			!if(TM .gt. 0)then
			i2=iebc+ie+i1
			hxbcl(i1,1)=hxf(i1,jbbc)
			hxbcr(i1,1)=hxf(i2,jbbc)
			Bhxbcl(i1,1)=Bhxf(i1,jbbc)	
			Bhxbcr(i1,1)=Bhxf(i2,jbbc)
			!endif

			!TE-mode
			if(TE .gt. 0)then
			exbcl(i1,1)=exf(i1,jbbc)
			exbcr(i1,1)=exf(i2,jbbc)
			DExbcl(i1,1)=DExf(i1,jbbc)	
			DExbcr(i1,1)=DExf(i2,jbbc)
			endif
4702 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		elseif( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! put back values to left & right at PML/PMl interface 	
		do 4712 i1=1,iebc

			!TM-mode
			!if(TM .gt. 0)then
			i2=iebc+ie+i1
			hxbcl(i1,jb)=hxb(i1,1)
			hxbcr(i1,jb)=hxb(i2,1)
			Bhxbcl(i1,jb)=Bhxb(i1,1)
			Bhxbcr(i1,jb)=Bhxb(i2,1)
			!endif

			!TE-mode
			if(TE .gt. 0)then
			exbcl(i1,jb)=exb(i1,1)
			exbcr(i1,jb)=exb(i2,1)
			DExbcl(i1,jb)=DExb(i1,1)
			DExbcr(i1,jb)=DExb(i2,1)
			endif
4712 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 


!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

    ! put BC Hx (Ex),Hy (Ey) upml/grid interface to grid Hx, Hy

		! front
		do  4704 i=i1c,i2c
			i2=iebc+i
			!TM-mode
			!if(TM .gt. 0)then
			Hx(i,1)=Hxf(i2,jbbc)
			!endif
			!TE-mode
			if(TE .gt. 0)then
			Ex(i,1)=Exf(i2,jbbc)
			endif
4704 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		elseif( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

		! back
		do  4714 i=i1c,i2c
			i2=iebc+i
			!TM-mode
			!if(TM .gt. 0)then
			Hx(i,jb)=Hxb(i2,1)
			!endif
			!TE-mode
			if(TE .gt. 0)then
			Ex(i,jb)=Exb(i2,1)
			endif
4714 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

		! left
		do  4706 j=j1c,j2c
			!TM-mode
			!if(TM .gt. 0)then
			Hy(1,j) =Hybcl(ibbc,j)
			!endif
			!TE-mode
			if(TE .gt. 0)then
			Ey(1,j) =Eybcl(ibbc,j)
			endif  !TE-mode
4706 continue

		! right
		do  4716 j=j1c,j2c
			!TM-mode
			!if(TM .gt. 0)then
			Hy(ib,j)=Hybcr(1,j)
			!endif
			!TE-mode
			if(TE .gt. 0)then
			Ey(ib,j)=Eybcr(1,j)
			endif  !TE-mode
4716 continue

		!-------------------------------------------------------------------
		!		Source for Hx,Hy (Ex,Ey)
		!-------------------------------------------------------------------
	do 2014 m=1,num_wg
		do 2015 k=1,nfield(m)
			if ((inputfield(m,k) .eq. 1) .and.(k .eq. 1) ) then	! First Pulse
				tdelay=tdelay0+1.5*pwidth(m,k)+tdel(m,k)
				beamshape=exp(-(tt-tdelay)**2/(pwidth(m,k)/2)**2)
			elseif (inputfield(m,k) .eq. 1) then	!second pulse
				tdelay=tdelay0+1.5*pwidth(m,k)+tdel(m,1)+tdel(m,2)
				beamshape=exp(-(tt-tdelay)**2/(pwidth(m,k)/2)**2)
			elseif(tt .lt. tdel(m,k))then	  !CW
				beamshape=0.
			else
				beamshape=1.
			endif

		!-----------------------------------------------------------
		!		Total Field/Scattered field for TM-mode (Hx, Hy)
		!-----------------------------------------------------------
				
			if ((pol(m,k) .eq. 1) .and. (wgdirection(m) .eq. 1) .and. (ilaunch(m,k) .gt. 0)) then
				do 2020 j=js1,js2
					Ez_inc=beamshape*E_input(m,k)*sin(w_f(m,k)*tt)*rmode(m,k,j)
					Hy(ilaunch(m,k),j)=Hy(ilaunch(m,k),j)-db*Ez_inc
2020			continue
			
			elseif ((pol(m,k) .eq. 1) .and. (wgdirection(m) .eq. 1) .and. (ilaunch(m,k) .lt. 0)) then
				do 2021 j=js1,js2
					Ez_inc=beamshape*E_input(m,k)*sin(w_f(m,k)*tt)*rmode(m,k,j)
					ilaunchmk = ie + ilaunch(m,k)
					Hy(ilaunchmk,j)=Hy(ilaunchmk,j)+db*Ez_inc
2021			continue

			! front only
			elseif ((pol(m,k) .eq. 1) .and. (wgdirection(m) .eq. 2) .and. (ilaunch(m,k) .gt. 0)) then
!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

				do 2022 i=is1,is2
					Ez_inc=beamshape*E_input(m,k)*sin(w_f(m,k)*tt)*rmode(m,k,i)
					Hx(i,ilaunch(m,k))=Hx(i,ilaunch(m,k))+db*Ez_inc
2022		continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

			! back only
			elseif ((pol(m,k) .eq. 1) .and. (wgdirection(m) .eq. 2) .and. (ilaunch(m,k) .lt. 0)) then

!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

				do 2023 i=is1,is2
					Ez_inc=beamshape*E_input(m,k)*sin(w_f(m,k)*tt)*rmode(m,k,i)
					ilaunchmk = je + ilaunch(m,k)
					Hx(i,ilaunchmk)=Hx(i,ilaunchmk)-db*Ez_inc
2023		continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

		!-----------------------------------------------------------
		!		Total Field/Scattered field for TE-mode (Ex, Ey)
		!-----------------------------------------------------------

			elseif ((pol(m,k) .eq. 0) .and. (wgdirection(m) .eq. 1) .and. (ilaunch(m,k) .gt. 0)) then
				do 2024	j=js1,js2
					Hz_inc=beamshape*(kprop(m,k)/w_f(m,k)/muz)*E_input(m,k)*sin(w_f(m,k)*tt)*rmode(m,k,j)
					Ey(ilaunch(m,k),j)=Ey(ilaunch(m,k),j)+cb*geometry(ilaunch(m,k),j)*Hz_inc
2024			continue

			elseif ((pol(m,k) .eq. 0) .and. (wgdirection(m) .eq. 1) .and. (ilaunch(m,k) .lt. 0)) then
				do 2025	j=js1,js2
					ilaunchmk = ie + ilaunch(m,k)
					Hz_inc=beamshape*(kprop(m,k)/w_f(m,k)/muz)*E_input(m,k)*sin(w_f(m,k)*tt)*rmode(m,k,j)
					Ey(ilaunchmk,j)=Ey(ilaunchmk,j)-cb*geometry(ilaunchmk,j)*Hz_inc
2025			continue

			! front only
			elseif ((pol(m,k) .eq. 0) .and. (wgdirection(m) .eq. 2) .and. (ilaunch(m,k) .gt. 0)) then
!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. 0 ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

				do 2026	i=is1,is2
					Hz_inc=-beamshape*(kprop(m,k)/w_f(m,k)/muz)*E_input(m,k)*sin(w_f(m,k)*tt)*rmode(m,k,i)
					Ex(i,ilaunch(m,k))=Ex(i,ilaunch(m,k))-cb*geometry(i,ilaunch(m,k))*Hz_inc
2026			continue
!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

			! back only
			else
!DEC$ IF DEFINED (USE_MPI_CH) 
		if( pid .eq. (numprocs-1) ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

				do 2027	i=is1,is2
					Hz_inc=-beamshape*(kprop(m,k)/w_f(m,k)/muz)*E_input(m,k)*sin(w_f(m,k)*tt)*rmode(m,k,i)
					ilaunchmk = je + ilaunch(m,k)
					Ex(i,ilaunchmk)=Ex(i,ilaunchmk)+cb*geometry(i,ilaunchmk)*Hz_inc
2027			continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

			endif

2015	continue
2014	continue

if(TE .gt. 0)then
	do 103 j=j1c,j2c
	do 103 i=i1c,i2c
	
		geometryij=geometry(i,j)
		
		if (gain .gt. 0) then
		gainij=gaingeo(i,j)
		else
		gainij=0
		endif

      if (gainij.eq.0)then
	  else
		Exgij= (Ex(i,j)+Ex(i,j+1))/2.
		Exg_oldij= (Ex_old(i,j)+Ex_old(i,j+1))/2.
		Eygij= (Ey(i,j)+Ey(i+1,j))/2.
		Eyg_oldij= (Ey_old(i,j)+Ey_old(i+1,j))/2.
		Pxg_1ij(:)=(px_1(i,j,:)+Px_1(i,j+1,:))/2.
		Pxg_old1_1ij(:)=(px_old1_1(i,j,:)+Px_old1_1(i,j+1,:))/2.
		Pxg_2ij(:)=(px_2(i,j,:)+Px_2(i,j+1,:))/2.
		Pxg_old1_2ij(:)=(px_old1_2(i,j,:)+Px_old1_2(i,j+1,:))/2.
		Pyg_1ij(:)=(py_1(i,j,:)+Py_1(i+1,j,:))/2.
		Pyg_old1_1ij(:)=(py_old1_1(i,j,:)+Py_old1_1(i+1,j,:))/2.
		Pyg_2ij(:)=(Py_2(i,j,:)+Py_2(i+1,j,:))/2.
		Pyg_old1_2ij(:)=(Py_old1_2(i,j,:)+Py_old1_2(i+1,j,:))/2.
		Avptxgij= (Avptx(i,j)+Avptx(i,j+1))/2.
		Avptxg_oldij= (Avptx_old(i,j)+Avptx_old(i,j+1))/2.
		Avptygij= (Avpty(i,j)+Avpty(i+1,j))/2.
		Avptyg_oldij= (Avpty_old(i,j)+Avpty_old(i+1,j))/2.

		p=gainij
	  	 if(interaction .eq. 1)then
			EdP2ij(p)=((Exgij+Exg_oldij)*(Pxg_2ij(p)-Pxg_old1_2ij(p))+		&
				(Eygij+Eyg_oldij)*(Pyg_2ij(p)-Pyg_old1_2ij(p))+		&
				(Ez(i,j)+Ez_old(i,j))*(Pz_2(i,j,p)-Pz_old1_2(i,j,p)))*C_n2_2(p)
 			EdP1ij=((Exgij+Exg_oldij)*(Pxg_1ij(p)-Pxg_old1_1ij(p))+		&
				(Eygij+Eyg_oldij)*(Pyg_1ij(p)-Pyg_old1_1ij(p))+		&
				(Ez(i,j)+Ez_old(i,j))*(Pz_1(i,j,p)-Pz_old1_1(i,j,p)))*C_n1_2(p)
  
			if (ipauli .eq. 2)then
			n3(i,j,p)=n3(i,j,p)*(1.-dt_over_tau30(p)*(1.-n0(i,j,p)))-n3(i,j,p)*(1.-n2(i,j,p))*dt_over_tau32(p) &
			  +EdP2ij(p)+pumping_rate(p)*(1-n3(i,j,p))
			n2(i,j,p)=n2(i,j,p)*(1.-dt_over_tau21(p)*(1.-n1(i,j,p)))+n3(i,j,p)*(1.-n2(i,j,p))*dt_over_tau32(p) & 
			  +EdP1ij(p)
			n1(i,j,p)=n1(i,j,p)*(1.-dt_over_tau10(p)*(1.-n0(i,j,p)))+n2(i,j,p)*(1.-n1(i,j,p))*dt_over_tau21(p) & 
			  -EdP1ij(p)
			n0(i,j,p)=ne(p)-n1(i,j,p)-n2(i,j,p)-n3(i,j,p)-pumping_rate(p)*(1-n3(i,j,p))
			endif
		 else      !A.p

			AP2ij(p)=C_n2_1(P)*(Avptxgij*Pxg_2ij(p)+Avptygij*Pyg_2ij(p)+Avptz(i,j)*Pz_2(i,j,p))
 			AP1ij(p)=C_n1_1(P)*(Avptxgij*Pxg_1ij(p)+Avptygij*Pyg_1ij(p)+Avptz(i,j)*Pz_1(i,j,p))

			if(ipauli .eq. 2)then
			n3(i,j,p)=n3(i,j,p)*(1.-dt_over_tau30(p)*(1.-n0(i,j,p)))-n3(i,j,p)*(1.-n2(i,j,p))*dt_over_tau32(p) & 
				-AP2ij(p)+pumping_rate(p)*(1-n3(i,j,p))
			n2(i,j,p)=n2(i,j,p)*(1.-dt_over_tau21(p)*(1.-n1(i,j,p)))+n3(i,j,p)*(1.-n2(i,j,p))*dt_over_tau32(p) &
				-AP1ij(p)
			n1(i,j,p)=n1(i,j,p)*(1.-dt_over_tau10(p)*(1.-n0(i,j,p)))+n2(i,j,p)*(1.-n1(i,j,p))*dt_over_tau21(p) & 
				+AP1ij(p)
			n0(i,j,p)=ne(p)-n1(i,j,p)-n2(i,j,p)-n3(i,j,p)-pumping_rate(p)*(1-n3(i,j,p))
			endif
		 endif
		endif
103	continue
endif	!TE-mode

!write(*,*) "H source update" 
!DEC$ IF DEFINED (USE_MPI_CH) 
		! (m) main grid: exchange lines (j=..) between processes before getting Hx, etc. ..
		! send bottom line of top process to bottom process 
		! no data exchange needed for calculating Hy, Ey ..

    call MPI_BARRIER( comm_b, ierr ) ! synchronize entrance;

		!if(TM .gt. 0)then
			! (m) for Ez, Hx 2 -> 1, 1 -> 0, etc.	
			call exchng4( Hx, ie, j1_al, j2_al+1, je_np, pid_b, nbrtop_b, &
										nbrbottom_b, numprocs, n, comm_b )

			! (upml)			HXbcl, HXbcr 2 -> 1, 1 -> 0, etc.	
			call exchng4( HXbcl, iebc, j1_al, j2_al+1, je_np, pid_b, nbrtop_b, &
										nbrbottom_b, numprocs, n, comm_b )
			call exchng4( HXbcr, iebc, j1_al, j2_al+1, je_np, pid_b, nbrtop_b, &
										nbrbottom_b, numprocs, n, comm_b )
		!endif

		if(TE .gt. 0)then
			! (m) for Hz, Ex 2 -> 1, 1 -> 0, etc.	
			call exchng4( Ex, ie, j1_al, j2_al+1, je_np, pid_b, nbrtop_b, &
										nbrbottom_b, numprocs, n, comm_b )
			! (upml)		  EXbcl, EXbcr 2 -> 1, 1 -> 0, etc.	
			call exchng4( EXbcl, iebc, j1_al, j2_al+1, je_np, pid_b, nbrtop_b, &
										nbrbottom_b, numprocs, n, comm_b )
			call exchng4( EXbcr, iebc, j1_al, j2_al+1, je_np, pid_b, nbrtop_b, &
										nbrbottom_b, numprocs, n, comm_b )
		endif

!DEC$ ENDIF ! (USE_MPI_CH) 

!write(*,*) "H exchanged" 
	! ======================== file output ===================================

	! collect observer data only when tempdata output option is 1.  
	if (tempdata .eq. 1) then

! read observer point values to array(s)
	do 901 u = 1,nobserver
!DEC$ IF DEFINED (USE_MPI_CH) 
	if((jobs(u).ge.j1c .and. jobs(u).le.j2c).and.(iobs(u).ge.i1c .and. iobs(u).le.i2c))then
!DEC$ ENDIF ! (USE_MPI_CH) 

		if(TE .eq. 0)then
			Ez_obs(n_st, u) = Ez(iobs(u), jobs(u))
		elseif(TM .eq. 0)then
			Ex_obs(n_st, u) = (Ex(iobs(u), jobs(u))+Ex(iobs(u), jobs(u)+1))/2.
			Ey_obs(n_st, u) = (Ey(iobs(u), jobs(u))+Ey(iobs(u)+1, jobs(u)))/2.
		else
			Ez_obs(n_st, u) = Ez(iobs(u), jobs(u))
			Ex_obs(n_st, u) = (Ex(iobs(u), jobs(u))+Ex(iobs(u), jobs(u)+1))/2.
			Ey_obs(n_st, u) = (Ey(iobs(u), jobs(u))+Ey(iobs(u)+1, jobs(u)))/2.
		endif
		
		if (gain .gt. 0) then
		n0_obs(n_st, u) = sum(n0(iobs(u), jobs(u),:))
		n1_obs(n_st, u) = sum(n1(iobs(u), jobs(u),:))
		n2_obs(n_st, u) = sum(n2(iobs(u), jobs(u),:))
		n3_obs(n_st, u) = sum(n3(iobs(u), jobs(u),:))
		endif
!DEC$ IF DEFINED (USE_MPI_CH) 
	endif
!DEC$ ENDIF ! (USE_MPI_CH) 

901 continue
!!
	m_obs(n_st) = n ! record time step
	tt_obs(n_st) = tt ! record time
	n_st = n_st+1 ! increment 

	if( n_st .gt. ns_obs .or. n .eq. nmax ) then 	! (W) -> write to file

!DEC$ IF DEFINED (USE_MPI_CH) 
	! for mpich - transfer values to master

    call MPI_BARRIER( comm_b, ierr ) ! synchronize entrance;

		if(TE .eq. 0)then
			call obs_dat_to_master(Ez_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
		elseif(TM .eq. 0)then
			call obs_dat_to_master(Ex_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
			call obs_dat_to_master(Ey_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
		else
			call obs_dat_to_master(Ez_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
			call obs_dat_to_master(Ex_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
			call obs_dat_to_master(Ey_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
		endif

		if (gain .gt. 0) then
			call obs_dat_to_master(n0_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
			call obs_dat_to_master(n1_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
			call obs_dat_to_master(n2_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
			call obs_dat_to_master(n3_obs, ns_obs, nobserver, obs_in_proc, pid_b, comm_b )
		endif
!DEC$ ENDIF ! (USE_MPI_CH) 

	! write values .. 
904 format(i8,a3,e13.6e3,\)
905 format(a3,e13.6e3,\)
 

!DEC$ IF DEFINED (USE_MPI_CH) 
! for mpich write only if you are master ..
		if( pid .eq. master ) then 
!DEC$ ENDIF ! (USE_MPI_CH) 

do 902 v = 1, (n_st-1)

	write(122,904) m_obs(v),',',tt_obs(v) ! step number and time

	if (gain .gt. 0) then
		write(123,904) m_obs(v),',',tt_obs(v) ! step number and time
		write(124,904) m_obs(v),',',tt_obs(v) ! step number and time
		write(125,904) m_obs(v),',',tt_obs(v) ! step number and time
		write(126,904) m_obs(v),',',tt_obs(v) ! step number and time

!
		write(123,905) ( ( ',',n0_obs(v, u) ), u = 1,nobserver ) 
		write(124,905) ( ( ',',n1_obs(v, u) ), u = 1,nobserver ) 
		write(125,905) ( ( ',',n2_obs(v, u) ), u = 1,nobserver ) 
		write(126,905) ( ( ',',n3_obs(v, u) ), u = 1,nobserver ) 
	endif

	if(TE .eq. 0)then
		write(122,905) ( ( ',', Ez_obs(v, u) ), u = 1,nobserver ) 
	elseif(TM .eq. 0)then
		write(122,905) ( ( ',', Ex_obs(v, u) ), u = 1,nobserver ) 
		write(122,905) ( ( ',', Ey_obs(v, u) ), u = 1,nobserver ) 
	else
		write(122,905) ( ( ',', Ez_obs(v, u) ), u = 1,nobserver ) 
		write(122,905) ( ( ',', Ex_obs(v, u) ), u = 1,nobserver ) 
		write(122,905) ( ( ',', Ey_obs(v, u) ), u = 1,nobserver ) 
	endif
!
	write(122,'(a1)') ''! end of line
	
	if (gain .gt. 0) then
		write(123,'(a1)') ''! end of line
		write(124,'(a1)') ''! end of line
		write(125,'(a1)') ''! end of line
		write(126,'(a1)') ''! end of line
	endif

902 continue

!DEC$ IF DEFINED (USE_MPI_CH) 
		endif 
!DEC$ ENDIF ! (USE_MPI_CH) 

	! reset observer arrays, counter ..
	if (gain .gt. 0) then
		n0_obs(:,:)=0
		n1_obs(:,:)=0
		n2_obs(:,:)=0
		n3_obs(:,:)=0
	endif

	if(TE .eq. 0)then
			Ez_obs(:,:) = 0.
	elseif(TM .eq. 0)then
			Ex_obs(:,:) = 0.
			Ey_obs(:,:) = 0.
	else
			Ez_obs(:,:) = 0.
			Ex_obs(:,:) = 0.
			Ey_obs(:,:) = 0.
	endif
	!
	n_st = 1
	tt_obs(:) = 0.
	m_obs(:) = 0

	endif ! <- (W) 

	endif	! tempdata output .eq. 1
	 
!write(*,*) "Enter detector part"	
!!!----------------------------------------------------------------------
	!		 detector part
	!----------------------------------------------------------------------
do 1234 j=j1c,j2c
do 1234 i=i1c,i2c
  if (detector(i,j).gt.0) then
	do 1235 k=1,det_num
		if (detector(i,j).eq.k) then
		if(det_dir(k).eq.1)then
			if (TE .eq. 0)then
			det_result(k)=det_result(k)+(-Ez(i,j)*Hy(i,j))
			elseif(TM .eq. 0) then
			det_result(k)=det_result(k)+(Ey(i,j)*Hz(i,j))
			else
			det_result(k)=det_result(k)+(-Ez(i,j)*Hy(i,j)+Ey(i,j)*Hz(i,j))
	   		endif
		else
			if (TE .eq. 0)then
			det_result(k)=det_result(k)+(Ez(i,j)*Hx(i,j))
			elseif(TM .eq. 0) then
			det_result(k)=det_result(k)+(-Ex(i,j)*Hz(i,j))
			else
			det_result(k)=det_result(k)+(Ez(i,j)*Hx(i,j)-Ex(i,j)*Hz(i,j))
	   		endif
		endif
	endif								 																
1235 continue
  endif
 1234 continue

det_time=det_time+dt

if (det_time .ge. det_period) then ! write detector data

!DEC$ IF DEFINED (USE_MPI_CH)

  call MPI_BARRIER( comm_a, ierr ) ! synchronize entrance;

	call MPI_Allreduce( det_result, det_result_t, det_num, MPI_REAL,	&
                           MPI_SUM, MPI_COMM_WORLD, ierr ) 
 
	if( pid .eq. master ) then
		write(1345,904) n,',',tt ! 
		write(1345,905) ( ( ',', det_result_t(u) ), u = 1,det_num )
		write(1345,'(a1)') '' ! end if line
	endif
	det_result_t(:)=0.

!DEC$ ELSE ! (USE_MPI_CH) 
	write(1345,904) n,',',tt ! 
	write(1345,905) ( ( ',', det_result(u) ), u = 1,det_num )
	write(1345,'(a1)') '' ! end if line 

!DEC$ ENDIF ! (USE_MPI_CH) 
	
	det_result(:)=0.
	det_time=0
endif


!write(*,*) "Enter integrated spectrum detector"
!------------------------------------------------------------------------------------
!			integrated spectrum detector
!-------------------------------------------------------------------------------------
if(ftdet_num .ge. 1)then
	!special time limit for DFT of the input pulse
	ndft_input=int(2*(tdelay0+1.5*pwidth(1,1)+tdel(1,1))/dt)

	do 1236 k=1,det_num
		
		if((spec_index(k) .eq. 1) .and. (det_dir(k) .eq. 1))then
			jj=0
			do 1137 j=det_yp(k),det_yp(k)+det_ylen(k)-1,n_skip(k)
				id=det_xp(k)
				jj=jj+1
				if((j.ge.j1c).and.(j.le.j2c).and.(id.ge.i1c).and.(id.le.i2c))then
				if((k .ne. 1) .or. (n .le. ndft_input))then
					do 1237 iff=0,ifmax
						ffreq=wstart + float(iff)*wspan/float(ifmax)
						ftarg=ffreq*float(n)*dt
						ftcos=cos(ftarg)
						ftsin=sin(ftarg)
						if(TM .gt. 0)then
    						ftrd_z(iff,jj,k) = ftrd_z(iff,jj,k) + Ez(id,j)*ftcos
							ftid_z(iff,jj,k) = ftid_z(iff,jj,k) + Ez(id,j)*ftsin
						endif
						if(TE .gt. 0)then
							ftrd_x(iff,jj,k) = ftrd_x(iff,jj,k) + Ex(id,j)*ftcos
							ftid_x(iff,jj,k) = ftid_x(iff,jj,k) + Ex(id,j)*ftsin
							ftrd_y(iff,jj,k) = ftrd_y(iff,jj,k) + Ey(id,j)*ftcos
							ftid_y(iff,jj,k) = ftid_y(iff,jj,k) + Ey(id,j)*ftsin
						endif
1237				continue
				endif
				endif
1137			continue
		endif

		if((spec_index(k) .eq. 1) .and. (det_dir(k) .eq. 2))then
			ii=0	
			do 1037 i=det_xp(k),det_xp(k)+det_xlen(k)-1,n_skip(k)
				jd=det_yp(k)
				ii=ii+1
				if((jd.ge.j1c).and.(jd.le.j2c).and.(i.ge.i1c).and.(i.le.i2c))then
				if((k .ne. 1) .or. (n .le. ndft_input))then
					do 1337 iff=0,ifmax
						ffreq=wstart + float(iff)*wspan/float(ifmax)
						ftarg=ffreq*float(n)*dt
						ftcos=cos(ftarg)
						ftsin=sin(ftarg)
						if(TM .gt. 0)then
    						ftrd_z(iff,ii,k) = ftrd_z(iff,ii,k) + Ez(i,jd)*ftcos
							ftid_z(iff,ii,k) = ftid_z(iff,ii,k) + Ez(i,jd)*ftsin
						endif
						if(TE .gt. 0)then
							ftrd_x(iff,ii,k) = ftrd_x(iff,ii,k) + Ex(i,jd)*ftcos
							ftid_x(iff,ii,k) = ftid_x(iff,ii,k) + Ex(i,jd)*ftsin
							ftrd_y(iff,ii,k) = ftrd_y(iff,ii,k) + Ey(i,jd)*ftcos
							ftid_y(iff,ii,k) = ftid_y(iff,ii,k) + Ey(i,jd)*ftsin
						endif
1337				continue
				endif
				endif
1037		continue
		endif

1236	continue

endif  !END OF INTEGRATED SPECTRUM DETECTOR


!write(*,*) "finished integrated spectrum detector"

if(fluxmap .eq. 1)then
!----------------------------------------------------------------------
!		 s_x, s_y, s flux map part, gets flux map over grid (2D)  
!----------------------------------------------------------------------
	! writes output only once (before exit)
	! number of time steps detector is on (summing up on steps over half period)
	! num_ - number, ts_ - time steps, _fd, _f_d - flux detector; 
	

! > $2761-1
	if ((n.gt.(n_fluxmap-num_ts_fd-1)).and.(n.le.n_fluxmap)) then 

		do 1240 j = j1c,j2c
		do 1240 i = i1c,i2c
			if (TE.eq.0)then	 ! (TM)
				s_x_f_d(i,j) = s_x_f_d(i,j) - Ez(i,j)*Hy(i,j)
				s_y_f_d(i,j) = s_y_f_d(i,j) + Ez(i,j)*Hx(i,j)
			elseif (TM.eq.0)then ! (TE)
				s_x_f_d(i,j) = s_x_f_d(i,j) + Ey(i,j)*Hz(i,j)
				s_y_f_d(i,j) = s_y_f_d(i,j) - Ex(i,j)*Hz(i,j)
			else							 ! (TE, TM)
				s_x_f_d(i,j) = s_x_f_d(i,j) + (-Ez(i,j)*Hy(i,j)+Ey(i,j)*Hz(i,j))
				s_y_f_d(i,j) = s_y_f_d(i,j) + ( Ez(i,j)*Hx(i,j)-Ex(i,j)*Hz(i,j))
		   	endif

			s_f_d(i,j) = s_f_d(i,j) + sqrt( s_x_f_d(i,j)**2 + s_y_f_d(i,j)**2 )

1240	continue

		! after number of time steps = num_ts_fd, it is time to write flux detector output
		if (n.eq.n_fluxmap) then
			

!DEC$ IF DEFINED (USE_MPI_CH) 
			if(pid .eq. master)then

				s_x_f_d_max = 0.
				s_y_f_d_max = 0.
				s_f_d_max = 0.

				s_x_f_d_fid = 1356 ! file ids ( _fid.. )
				s_y_f_d_fid = 1358
				s_f_d_fid = 1360
		!
				open(s_x_f_d_fid, file= (trim(datadir)//'s_x_f_d.csv'))
				open(s_y_f_d_fid, file= (trim(datadir)//'s_y_f_d.csv'))
				open(s_f_d_fid, file= (trim(datadir)//'s_f_d.csv'))

				allocate(s_x_f_d_r(ie,dje_np),s_y_f_d_r(ie,dje_np),s_f_d_r(ie,dje_np))
	
			endif
		
		
			do 9985 u=0,numprocs-1
				call MPI_BARRIER( comm_a, ierr )
				if(u .ne. 0)then
					nsend = (je_np(u+1) - je_np(u))*ie
						if( pid .eq. u ) then	
							call MPI_ISEND(s_x_f_d(1,je_np(u)+1),nsend, MPI_REAL, int(0), int(u), &
									comm_a, req, ierr )		
						endif						
						if( pid .eq. 0 ) then 
							call MPI_IRECV(s_x_f_d_r(1,1), nsend, MPI_REAL, int(u), &
								int(u), comm_a, req, ierr )	
						endif
						call MPI_WAIT ( req, status, ierr )	
						if( pid .eq. u ) then	
							call MPI_ISEND(s_y_f_d(1,je_np(u)+1),nsend, MPI_REAL, int(0), numprocs+int(u), &
									comm_a, req, ierr )		
						endif						
						if( pid .eq. 0 ) then 
							call MPI_IRECV(s_y_f_d_r(1,1), nsend, MPI_REAL, int(u), &
								numprocs+int(u), comm_a, req, ierr )	
						endif
						call MPI_WAIT ( req, status, ierr )	
						if( pid .eq. u ) then	
							call MPI_ISEND(s_f_d(1,je_np(u)+1),nsend, MPI_REAL, int(0), numprocs+int(u), &
									comm_a, req, ierr )		
						endif						
						if( pid .eq. 0 ) then 
							call MPI_IRECV(s_f_d_r(1,1), nsend, MPI_REAL, int(u), &
								numprocs+int(u), comm_a, req, ierr )	
						endif		
						call MPI_WAIT ( req, status, ierr )	
				else
					if(pid .eq. 0)then
						do 8871 j=1,je_np(1)
						do 8871 i=1,ie	
								s_x_f_d_r(i,j)=s_x_f_d(i,j)
								s_y_f_d_r(i,j)=s_y_f_d(i,j)
								s_f_d_r(i,j)=s_f_d(i,j)
8871					continue
					endif
				endif

				call MPI_BARRIER( comm_a, ierr )
				if(pid.eq.0)then
					write(*,*) 'MPI-send & receive finished for process:', u 
				endif

				if(pid.eq.0)then
					write(*,*) 'Writing flux data from process..... ',u
					if(mod(je_np(u)+1,reso_snapshot_y).eq.0)then
						js_np=reso_snapshot_y
					else
						js_np=mod(je_np(u)+1,reso_snapshot_y)
					endif
					do 9981 j=js_np,je_np(u+1)-je_np(u),reso_snapshot_y
					do 9981 i=1,ie,reso_snapshot_x
							write(s_x_f_d_fid,*) s_x_f_d_r(i,j)
							if(s_x_f_d_r(i,j) .gt. s_x_f_d_max) then
								s_x_f_d_max = s_x_f_d_r(i,j)
							endif
							write(s_y_f_d_fid,*) s_y_f_d_r(i,j)
							if(s_y_f_d_r(i,j) .gt. s_y_f_d_max) then
								s_y_f_d_max = s_y_f_d_r(i,j)
							endif
							write(s_f_d_fid,*) s_f_d_r(i,j)
							if(s_f_d_r(i,j) .gt. s_f_d_max) then
								s_f_d_max = s_f_d_r(i,j)
							endif						
9981				continue
					write(*,*) '          ........Finished.'			
				endif !for pid=0	
9985		continue

			call MPI_BARRIER( comm_a, ierr )

			if(pid.eq.0)then
					close(s_x_f_d_fid)
					close(s_y_f_d_fid)
					close(s_f_d_fid)
					deallocate(s_x_f_d_r,s_y_f_d_r,s_f_d_r)
			endif


!DEC$ ELSE ! (USE_MPI_CH)

				s_x_f_d_max = 0.
				s_y_f_d_max = 0.
				s_f_d_max = 0.

				do 1242 j = 1,je
				do 1242 i = 1,ie

					if (s_x_f_d(i,j).gt.s_x_f_d_max) then
						s_x_f_d_max = s_x_f_d(i,j)
					endif
					if (s_y_f_d(i,j).gt.s_y_f_d_max) then
						s_y_f_d_max = s_y_f_d(i,j)
					endif
					if (s_f_d(i,j).gt.s_f_d_max) then
						s_f_d_max = s_f_d(i,j)
					endif

1242			continue

				s_x_f_d_fid = 1356 ! file ids ( _fid.. )
				s_y_f_d_fid = 1358
				s_f_d_fid = 1360
		!
				open(s_x_f_d_fid, file= (trim(datadir)//'s_x_f_d.csv'))
				open(s_y_f_d_fid, file= (trim(datadir)//'s_y_f_d.csv'))
				open(s_f_d_fid, file= (trim(datadir)//'s_f_d.csv'))

				! normalize output to max value and write to file
				! data in the range [0,1] is written
				write(*,*) " "
				write(*,*) "writing flux data...................."
				do 9994 j=1,je
				do 9994 i=1,ie
					write(s_x_f_d_fid,*)  s_x_f_d(i,j)
					write(s_y_f_d_fid,*)  s_y_f_d(i,j)
					write(s_f_d_fid,*)	  s_f_d(i,j)
9994			continue

				close(s_x_f_d_fid)
				close(s_y_f_d_fid)
				close(s_f_d_fid)
	
!DEC$ ENDIF ! (USE_MPI_CH) 

			s_x_f_d(:,:) = 0.
			s_y_f_d(:,:) = 0.
			s_f_d(:,:) = 0.

		endif ! n .eq. n_fluxmap

	endif !(n.gt.(n_fluxmap-num_ts_fd-1)).and.(n.le.n_fluxmap)
! < $2761-1
endif !fluxmap output

!write(*,*) "passed flux output"

	!----------------------------------------------------------------------
	!		 output snap shot
	!----------------------------------------------------------------------

	if( (mod(n, nsnapshotinterval) .eq. 0).and. (count .lt. nsnapshot )) then
		count=count+1
!write(*,*) 'enter snap shot'
! medium
!DEC$ IF DEFINED (USE_MPI_CH)
	if (gain .gt. 0) then
		do 2244 j=j1_al,j2_al
		do 2244 i=1,ie
			n0_result(i,j)=sum(n0(i,j,:))
			n1_result(i,j)=sum(n1(i,j,:))
			n2_result(i,j)=sum(n2(i,j,:))
			n3_result(i,j)=sum(n3(i,j,:))
2244    continue
	endif
!write(*,*) 'stop 1'
!DEC$ ELSE
	if (gain .gt. 0) then
		do 2244 j=1,je
		do 2244 i=1,ie
			n0_result(i,j)=sum(n0(i,j,:))
			n1_result(i,j)=sum(n1(i,j,:))
			n2_result(i,j)=sum(n2(i,j,:))
			n3_result(i,j)=sum(n3(i,j,:))
2244    continue
	endif
!DEC$ ENDIF ! (USE_MPI_CH)
! end medium

!DEC$ IF DEFINED (USE_MPI_CH) 
!
!		call dbg_dump_2d( 'Ex_0___', Ex, ie, jb, pid_a )			
!
!		call MPI_BARRIER( comm_a, ierr ) ! synchronize entrance;
!
!		if(TM .gt. 0) then		! transfer Ez to master	
!			call array2d_to_master( Ez, ie, je, je_np, int(0), pid_a, numprocs, n, comm_a )
!		endif
!
!		if(TE .gt. 0) then		! transfer Ex, Ey to master	
!			call array2d_to_master( Ex, ie, jb, je_np, int(1), pid_a, numprocs, n, comm_a )
!			call array2d_to_master( Ey, ib, je, je_np, int(0), pid_a, numprocs, n, comm_a )
!		endif
!
!		call dbg_dump_2d( 'Ex_1___', Ex, ie, jb, pid_a )			
!		call exit
!
!DEC$ ENDIF ! (USE_MPI_CH) 

!DEC$ IF DEFINED (USE_MPI_CH)
		nfile1=6100+count
		nfile2=6300+count
		nfile3=6500+count
		nfile4=6700+count

		nfile11=4100+count
		nfile12=4300+count
		nfile13=4500+count
 
write(*,*) "Enter snap output"
		
		if(pid .eq. master)then
			if(TM .gt. 0)then
			allocate(ezr(ie,dje_np))
			open(nfile11, file=trim(afilename(count)))
			endif
	
			if(TE.gt.0)then
			allocate(Hzr(ie,dje_np))
			allocate(eyr(ie,dje_np))
			open(nfile12, file=trim(ffilename(count)))
			open(nfile13, file=trim(gfilename(count)))
			endif
write(*,*) 'stop 2'			
			if (gain .gt. 0) then
				allocate(n0r(ie,dje_np))
				allocate(n1r(ie,dje_np))
				allocate(n2r(ie,dje_np))
				allocate(n3r(ie,dje_np))
				open(nfile1, file=trim(bfilename(count)))
				open(nfile2, file=trim(cfilename(count)))
				open(nfile3, file=trim(dfilename(count)))
				open(nfile4, file=trim(efilename(count)))
			endif	
		endif
		do 9995 u=0,numprocs-1	!processor id
			call MPI_BARRIER( comm_a, ierr )
			if(u .ne. 0)then !not master
				nsend = (je_np(u+1) - je_np(u))*ie
	 			nsendx = (je_np(u+1) - je_np(u)+1)*ie
				nsendy = (je_np(u+1) - je_np(u))*ib
				
				if(TM.gt.0)then
					if( pid .eq. u ) then	
write(*,*) je_np(u)+1
						call MPI_ISEND(ez(1,je_np(u)+1),nsend, MPI_REAL, int(0), int(u), &
								comm_a, req, ierr )		
					endif						
					if( pid .eq. 0 ) then 
						call MPI_IRECV(ezr(1,1), nsend, MPI_REAL, int(u), &
							int(u), comm_a, req, ierr )	
					endif
					call MPI_WAIT ( req, status, ierr )	
				endif
				
				if(TE.gt.0)then
					if( pid .eq. u ) then	
						call MPI_ISEND(Hz(1,je_np(u)+1),nsend, MPI_REAL, int(0), int(u), &
								comm_a, req, ierr )
					endif
					if( pid .eq. 0 ) then 
						call MPI_IRECV(Hzr(1,1), nsend, MPI_REAL, int(u), &
							int(u), comm_a, req, ierr )
					endif
					call MPI_WAIT ( req, status, ierr )	
					if( pid .eq. u ) then	
						call MPI_ISEND(ey(1,je_np(u)+1),nsendy, MPI_REAL, int(0), numprocs+int(u), &
								comm_a, req, ierr )
					endif
					if( pid .eq. 0 ) then 
						call MPI_IRECV(eyr(1,1), nsendy, MPI_REAL, int(u), &
							numprocs+int(u), comm_a, req, ierr )
					endif
					call MPI_WAIT ( req, status, ierr )	
				endif
			
				if (gain .gt. 0) then
					if( pid .eq. u ) then	
						call MPI_ISEND(n0_result(1,je_np(u)+1),nsend, MPI_REAL, int(0), int(u), &
						comm_a, req, ierr )		
					endif						
					if( pid .eq. 0 ) then 
						call MPI_IRECV(n0r(1,1), nsend, MPI_REAL, int(u), &
						int(u), comm_a, req, ierr )	
					endif
					call MPI_WAIT ( req, status, ierr )	
					
					if( pid .eq. u ) then	
						call MPI_ISEND(n1_result(1,je_np(u)+1),nsend, MPI_REAL, int(0), int(u), &
						comm_a, req, ierr )		
					endif						
					if( pid .eq. 0 ) then 
						call MPI_IRECV(n1r(1,1), nsend, MPI_REAL, int(u), &
						int(u), comm_a, req, ierr )	
					endif
					call MPI_WAIT ( req, status, ierr )	
					
					if( pid .eq. u ) then	
						call MPI_ISEND(n2_result(1,je_np(u)+1),nsend, MPI_REAL, int(0), int(u), &
						comm_a, req, ierr )		
					endif						
					if( pid .eq. 0 ) then 
						call MPI_IRECV(n2r(1,1), nsend, MPI_REAL, int(u), &
						int(u), comm_a, req, ierr )	
					endif
					call MPI_WAIT ( req, status, ierr )	
					
					if( pid .eq. u ) then	
						call MPI_ISEND(n3_result(1,je_np(u)+1),nsend, MPI_REAL, int(0), int(u), &
								comm_a, req, ierr )		
					endif						
					if( pid .eq. 0 ) then 
						call MPI_IRECV(n3r(1,1), nsend, MPI_REAL, int(u), &
							int(u), comm_a, req, ierr )	
					endif
					call MPI_WAIT ( req, status, ierr )	
				endif
			else
				if(pid .eq. 0)then
					do 8881 j=1,je_np(1)
						do 8882 i=1,ie
							if(TM.gt.0)then
								ezr(i,j)=ez(i,j)
							endif
							if(TE.gt.0)then
								Hzr(i,j)=Hz(i,j)
								eyr(i,j)=ey(i,j)
							endif
8882					continue
					if(TE.gt.0)then			!gy added line
							eyr(ib,j)=ey(ib,j) 	
					endif                   !gy added line  
8881				continue
					endif
				endif

			call MPI_BARRIER( comm_a, ierr )
			if(pid.eq.0)then
				write(*,*) 'MPI-send & receive finished for process:', u 
			endif

			if(pid.eq.0)then
				write(*,*) 'Writing data from process..... ',u
				if(mod(je_np(u)+1,reso_snapshot_y).eq.0)then
					js_np=reso_snapshot_y
				else
					js_np=mod(je_np(u)+1,reso_snapshot_y)
				endif
				temp1=0.5*cc*epsz*rnf
				do 9991 j=js_np,je_np(u+1)-je_np(u),reso_snapshot_y
				do 9991 i=1,ie,reso_snapshot_x
					!medium
					if(gain .gt. 0) then
						write(nfile1,*) n0r(i,j)
						write(nfile2,*) n1r(i,j)
						write(nfile3,*) n2r(i,j)
						write(nfile4,*) n3r(i,j)
					endif
					!end medium
					if(TM .gt. 0) then				!TM
						if (outputchoice .eq. 1) then !output E
							write(nfile11,*) Ezr(i,j)
							if(Ezr(i,j) .gt. ezmax) then
								ezmax=Ezr(i,j)
							endif
							if(Ezr(i,j) .lt. ezmin) then
								ezmin=Ezr(i,j)
							endif
						endif
						if (outputchoice .eq. 2) then !output I
							iz= temp1*Ezr(i,j)*Ezr(i,j)
							write(nfile11,*) iz 
							if(iz .gt. izmax) then
								izmax=iz
							endif
							if(iz .lt. izmin) then
								izmin=iz
							endif
						endif
						if (outputchoice .eq. 3) then !output log(I), I can not be 0
							if(Ezr(i,j) .eq.0) then
								logiz=-50.0
							else
								logiz=log(temp1*Ezr(i,j)*Ezr(i,j))
							endif
							write(nfile11,*)  logiz
							if(logiz .gt. logizmax) then
								logizmax=logiz
							endif
							if(logiz .lt. logizmin) then
								logizmin=logiz
							endif
						endif
					endif !TM
					if(TE.gt.0)then	!TE
						if (outputchoice .eq. 1) then !output E
							temp2=Hzr(i,j) 
							temp3=(Eyr(i,j)+Eyr(i+1,j))/2.
							write(nfile12,*) temp2 
							write(nfile13,*) temp3 
							if(temp2 .gt. exmax) then
								exmax=temp2
							endif
							if(temp2 .lt. exmin) then
								exmin=temp2
							endif
							if(temp3 .gt. eymax) then
								eymax=temp3
							endif
							if(temp3 .lt. eymin) then
								eymin=temp3
							endif
						endif
						if (outputchoice .eq. 2) then !output I
							temp2=temp1*Hzr(i,j)*Hzr(i,j)
							temp3=temp1*((Eyr(i,j)+Eyr(i+1,j))/2.)**2

							write(nfile12,*) temp2 
							write(nfile13,*) temp3 
							if(temp2 .gt. ixmax) then
								ixmax=temp2
							endif
							if(temp2 .lt. ixmin) then
								ixmin=temp2
							endif
							if(temp3 .gt. iymax) then
								iymax=temp3
							endif
							if(temp3 .lt. iymin) then
								iymin=temp3
							endif
						endif
						if (outputchoice .eq. 3) then !output log(I)
							temp2=log(temp1*Hzr(i,j)*Hzr(i,j))
							temp3=log(temp1*((Eyr(i,j)+Eyr(i+1,j))/2.)**2)

							write(nfile12,*) temp2
							write(nfile13,*) temp3 
							if(temp2 .gt. logixmax) then
								logixmax=temp2
							endif
							if(temp2 .lt. logixmin) then
								logixmin=temp2
							endif
							if(temp3 .gt. logiymax) then
								logiymax=temp3
							endif
							if(temp3 .lt. logiymin) then
								logiymin=temp3
							endif
						endif	!output log(I)
					endif	!TE
9991			continue
				write(*,*) '          .......Finished.'			
			endif !for pid=0	
9995	continue

		call MPI_BARRIER( comm_a, ierr )
		if(pid.eq.0)then
			if(TM .gt. 0)then
				close(nfile11)
				deallocate(ezr)
			endif
			if(TE.gt.0)then
				deallocate(Hzr)
				deallocate(eyr)
				close(nfile12)
				close(nfile13)
			endif
		endif

!DEC$ ELSE ! (USE_MPI_CH)
		! filename counter for snapshot
		nfile11=4100+count
		nfile12=4300+count
		nfile13=4500+count

		nfile1=6100+count
		nfile2=6300+count
		nfile3=6500+count
		nfile4=6700+count

		if(TM .gt. 0)then
		open(nfile11, file=trim(afilename(count)))
		endif
	
		if(TE.gt.0)then
		open(nfile12, file=trim(ffilename(count)))
		open(nfile13, file=trim(gfilename(count)))
		endif

		if (gain .gt. 0) then
			open(nfile1, file=trim(bfilename(count)))
			open(nfile2, file=trim(cfilename(count)))
			open(nfile3, file=trim(dfilename(count)))
			open(nfile4, file=trim(efilename(count)))
		endif
		! output choice ( 1: E, 2: Intensity, 3:Log(Intensity) )
		if (outputchoice .eq. 1) then
			do 9991 j=1,je,reso_snapshot_y
			do 9991 i=1,ie,reso_snapshot_x
				!medium
				if(gain .gt. 0) then
					write(nfile1,*) n0_result(i,j)
					write(nfile2,*) n1_result(i,j)
					write(nfile3,*) n2_result(i,j)
					write(nfile4,*) n3_result(i,j)
				endif
				!end medium

				if(TM .gt. 0) then
					write(nfile11,*) Ez(i,j)
					if(Ez(i,j) .gt. ezmax) then
						ezmax=Ez(i,j)
					endif
					if(Ez(i,j) .lt. ezmin) then
						ezmin=Ez(i,j)
					endif
				endif 
				if(TE.gt.0)then
					write(nfile12,*) Hz(i,j) 
					write(nfile13,*) (Ey(i,j)+Ey(i+1,j))/2. 
					if((Ex(i,j)+Ex(i,j+1))/2. .gt. exmax) then
						exmax=(Ex(i,j)+Ex(i,j+1))/2.
					endif
					if((Ex(i,j)+Ex(i,j+1))/2. .lt. exmin) then
						exmin=(Ex(i,j)+Ex(i,j+1))/2.
					endif
					if((Ey(i,j)+Ey(i+1,j))/2. .gt. eymax) then
						eymax=(Ey(i,j)+Ey(i+1,j))/2.
					endif
					if((Ey(i,j)+Ey(i+1,j))/2. .lt. eymin) then
						eymin=(Ey(i,j)+Ey(i+1,j))/2.
					endif
				endif
				
				
9991		continue
			if(TM .gt. 0)then
				close(nfile11)
			endif
			if(TE.gt.0)then
				close(nfile12)
				close(nfile13)
			endif
			if (gain .gt. 0) then
				close(nfile1)
				close(nfile2)
				close(nfile3)
				close(nfile4)
			endif
		! outputchoice = 2 and 3 need to be updated
		elseif (outputchoice .eq. 2) then
			temp1=0.5*cc*epsz*rnf
			do 9992 j=1,je,reso_snapshot_y
			do 9992 i=1,ie,reso_snapshot_x
				iz= temp1*Ez(i,j)**2
				write(nfile11,*) iz 
				if(iz .gt. izmax) then
					izmax=iz
				endif
				if(iz .lt. izmin) then
					izmin=iz
9992			endif

			close(nfile11)
		elseif (outputchoice .eq. 3) then
			temp1=0.5*cc*epsz*rnf
			do 9993 j=1,je,reso_snapshot_y
			do 9993 i=1,ie,reso_snapshot_x
			    logiz= log(temp1*Ez(i,j)**2)
				write(nfile11,*)  logiz
				if(logiz .gt. logizmax) then
					logizmax=logiz
				endif
				if(logiz .lt. logizmin) then
					logizmin=logiz
9993			endif

			close(nfile11)

		endif
!DEC$ ENDIF ! (USE_MPI_CH) 
write(*,*) "passed output snapshot"
	endif !End of snapshot



500 CONTINUE 

!------------------------------------------------------------------------------
!
!		 END of TIME STEPPING LOOP
!
!------------------------------------------------------------------------------

!DEC$ IF DEFINED (USE_MPI_CH) 
	if( pid .eq. master ) then
!DEC$ ENDIF ! (USE_MPI_CH) 

	close(1345)
	close(122)
	close(123)
	close(124)
	close(125)
	close(126)


	open(102,file= (trim(datadir)//'src-INFO.csv'),ACCESS='APPEND')! access mode have to be 'APPEND'

	!--------------------------------------------------------
	!	writing maxE, minE or maxI, minI or maxLogI, minLogI
	!--------------------------------------------------------
	if(outputchoice .eq. 1) then
		write(102,*) 'maximum & minimum field output--------------------'
		write(102,*) 'maximum Ez',',',ezmax
		write(102,*) 'minimum Ez',',',ezmin
		write(102,*) 'maximum Ex',',',exmax
		write(102,*) 'minimum Ex',',',exmin
		write(102,*) 'maximum Ey',',',eymax
		write(102,*) 'minimum Ey',',',eymin

	elseif(outputchoice .eq.2)then
		write(102,*) 'maximum & minimum field output--------------------'
		write(102,*) 'maximum I',',',izmax
		write(102,*) 'minimum I',',',izmin
	else 
		write(102,*) 'maximum & minimum field output--------------------'
		write(102,*) 'maximum LogI',',',logizmax
		write(102,*) 'minimum LogI',',',logizmin
	endif
	if(fluxmap.eq.1)then
		write(102,*) 'Maximum of flux intensity--------------------'
		write(102,*) 'maximum Sx',',',s_x_f_d_max
		write(102,*) 'minimum Sy',',',s_y_f_d_max
		write(102,*) 'maximum S',',',s_f_d_max
	endif
	CALL FDATE(timeend)
	write(102,*) 'end time :',timeend
	write(102,'(a1)') ' ' ! eof
	close(102)

	write(*,*) 'start time :',timestart
	write(*,*) 'end time :',timeend


!----------------------------------------------------------------
!		 spectrum analysis 
!---------------------------------------------------------------- 

	if ((spectrumanalysis .eq. 1).and.(tempdata .eq. 1)) then 
			! initialize the fft matrix
			do 90 j=1,nobserver
				if(TM .gt. 0)then
				ezinc(j)=0.
				fte2max_z(j)=0.
				endif
				if(TE .gt. 0)then
				exinc(j)=0.
				eyinc(j)=0.
				fte2max_x(j)=0.
				fte2max_y(j)=0.
				endif
			do 90 iff=0,ifmax
				if(TM .gt. 0)then
				ftr_z(j,iff)=0.
				fti_z(j,iff)=0.
				fte2_z(j,iff)=0.
				endif
				if(TE .gt. 0)then
				ftr_x(j,iff)=0.
				fti_x(j,iff)=0.
				fte2_x(j,iff)=0.
				ftr_y(j,iff)=0.
				fti_y(j,iff)=0.
				fte2_y(j,iff)=0.
				endif
90			continue

			if(TE .eq. 0)then
				open(91,file= (trim(datadir)//'observerTM.csv'), access='sequential', &
					action='read', share='denywr')
			elseif(TM .eq. 0)then
				open(91,file= (trim(datadir)//'observerTE.csv'), access='sequential', &
					action='read', share='denywr')
			else
				open(91,file= (trim(datadir)//'observerTM_TE.csv'), access='sequential', &
					action='read', share='denywr')
			endif

			! special time limit for input spectrum
			ndft_input=int(2*(tdelay0+1.5*pwidth(1,1)+tdel(1,1))/dt)

			DO 510 N=1, nmax
  
				if(TE.eq.0)then
					read(91,991) il5,a1,ez01,a1,ezinc(1),a1,ezinc(2),a1,ezinc(3),a1,ezinc(4) &
					,a1,ezinc(5),a1,ezinc(6),a1,ezinc(7),a1,ezinc(8),a1,ezinc(9),a1,ezinc(10)
				elseif(TM.eq.0)then
					read(91,992) il5,a1,ez01,a1,exinc(1),a1,exinc(2),a1,exinc(3),a1,exinc(4) &
					,a1,exinc(5),a1,exinc(6),a1,exinc(7),a1,exinc(8),a1,exinc(9),a1,exinc(10),  &
					a1,eyinc(1),a1,eyinc(2),a1,eyinc(3),a1,eyinc(4) &
					,a1,eyinc(5),a1,eyinc(6),a1,eyinc(7),a1,eyinc(8),a1,eyinc(9),a1,eyinc(10)
				else
					 read(91,993) il5,a1,ez01,a1,ezinc(1),a1,ezinc(2),a1,ezinc(3),a1,ezinc(4), &
					a1,ezinc(5),a1,ezinc(6),a1,ezinc(7),a1,ezinc(8),a1,ezinc(9),a1,ezinc(10),	&
					a1,exinc(1),a1,exinc(2),a1,exinc(3),a1,exinc(4), &
					a1,exinc(5),a1,exinc(6),a1,exinc(7),a1,exinc(8),a1,exinc(9),a1,exinc(10),  &
					a1,eyinc(1),a1,eyinc(2),a1,eyinc(3),a1,eyinc(4), &
					a1,eyinc(5),a1,eyinc(6),a1,eyinc(7),a1,eyinc(8),a1,eyinc(9),a1,eyinc(10)
				endif
				! DFT spectrum analysis	-----------------
				do 2053 iff=0,ifmax
					ffreq=wstart + float(iff)*wspan/float(ifmax)
					ftarg=ffreq*float(n)*dt
					ftcos=cos(ftarg)
					ftsin=sin(ftarg)

				do 2056 j=1,nobserver
					if((j .ne. 1) .or. (n .le. ndft_input))then
						if(TM .gt. 0)then
    						ftr_z(j,iff) = ftr_z(j,iff) + ezinc(j)*ftcos
							fti_z(j,iff) = fti_z(j,iff) + ezinc(j)*ftsin
						endif
						if(TE .gt. 0)then
							ftr_x(j,iff) = ftr_x(j,iff) + exinc(j)*ftcos
							fti_x(j,iff) = fti_x(j,iff) + exinc(j)*ftsin
							ftr_y(j,iff) = ftr_y(j,iff) + eyinc(j)*ftcos
							fti_y(j,iff) = fti_y(j,iff) + eyinc(j)*ftsin
						endif
					endif
2056			continue
2053			continue
510			continue
			close (910)

991 format(i8,a3,e13.6e3	&
			,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3   &
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3)
		   
992 format(i8,a3,e13.6e3  &
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3	&
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3	&
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3	  &
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3)
 

993 format(i8,a3,e13.6e3  &
			,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3   &
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3	&
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3	  &
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3	  &
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3		&
		   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3)

			! spectrum analysis ------------------
			if(TM .gt. 0)then
			open(87, file= (trim(datadir)//'Ez-DFT.csv'))
			endif
			if(TE .gt. 0)then
			open(88, file= (trim(datadir)//'Ex-DFT.csv'))
			open(89, file= (trim(datadir)//'Ey-DFT.csv'))
			endif

			do 3030 iff=0,ifmax
				do 3055 j=1,nobserver
					if(TM .gt. 0)then
						fte2_z(j,iff)=ftr_z(j,iff)**2+fti_z(j,iff)**2
						if (fte2_z(j,iff) .gt. fte2max_z(j)) then
							fte2max_z(j)=fte2_z(j,iff)
						endif
 					endif
					if(TE .gt. 0)then
						fte2_x(j,iff)=ftr_x(j,iff)**2+fti_x(j,iff)**2
						if(fte2_x(j,iff) .gt. fte2max_x(j)) then
							fte2max_x(j)=fte2_x(j,iff)
						endif
						fte2_y(j,iff)=ftr_y(j,iff)**2+fti_y(j,iff)**2
						if (fte2_y(j,iff) .gt. fte2max_y(j)) then
							fte2max_y(j)=fte2_y(j,iff)
						endif
					endif
3055			continue
3030		continue
	
			do 30309 iff=0,ifmax
				ffreq=wstart + float(iff)*wspan/float(ifmax)
       			wavelength0=tpi*cc/ffreq
				do 30319 k=1,nobserver
					if(TM .gt. 0)then
						if (fte2max_z(k) .eq. 0)then
							fte2max_z(k)=1.
						endif
					endif
					if(TE .gt. 0)then
						if (fte2max_x(k) .eq. 0) Then
							fte2max_x(k)=1.
						endif
						if (fte2max_y(k) .eq. 0) Then
							fte2max_y(k)=1.
						endif
					endif
30319			continue
	
				!----output DFT spectrum with normalization and without
				if(TM .gt. 0)then 
				write(87,871)  wavelength0,',',ffreq/tpi	&
					,',',fte2_z(1,iff)/fte2max_z(1),',',fte2_z(2,iff)/fte2max_z(2)     &
					,',',fte2_z(3,iff)/fte2max_z(3),',',fte2_z(4,iff)/fte2max_z(4)  &
					,',',fte2_z(5,iff)/fte2max_z(5),',',fte2_z(6,iff)/fte2max_z(6)  &
					,',',fte2_z(7,iff)/fte2max_z(7),',',fte2_z(8,iff)/fte2max_z(8)  &
					,',',fte2_z(9,iff)/fte2max_z(9),',',fte2_z(10,iff)/fte2max_z(10)&
					,',',fte2_z(1,iff)/fte2_z(1,iff),',',fte2_z(2,iff)/fte2_z(1,iff) &
					,',',fte2_z(3,iff)/fte2_z(1,iff),',',fte2_z(4,iff)/fte2_z(1,iff) &
					,',',fte2_z(5,iff)/fte2_z(1,iff),',',fte2_z(6,iff)/fte2_z(1,iff)  &
					,',',fte2_z(7,iff)/fte2_z(1,iff),',',fte2_z(8,iff)/fte2_z(1,iff) &
					,',',fte2_z(9,iff)/fte2_z(1,iff),',',fte2_z(10,iff)/fte2_z(1,iff) &
					,',',fte2_z(1,iff)/fte2_in(iff),',',fte2_z(2,iff)/fte2_in(iff) &
					,',',fte2_z(3,iff)/fte2_in(iff),',',fte2_z(4,iff)/fte2_in(iff) &
					,',',fte2_z(5,iff)/fte2_in(iff),',',fte2_z(6,iff)/fte2_in(iff)  &
					,',',fte2_z(7,iff)/fte2_in(iff),',',fte2_z(8,iff)/fte2_in(iff)  &
					,',',fte2_z(9,iff)/fte2_in(iff),',',fte2_z(10,iff)/fte2_in(iff)
				endif
				if(TE .gt. 0)then
				write(88,871)  wavelength0,',',ffreq/tpi	&
					,',',fte2_x(1,iff)/fte2max_x(1),',',fte2_x(2,iff)/fte2max_x(2)     &
					,',',fte2_x(3,iff)/fte2max_x(3),',',fte2_x(4,iff)/fte2max_x(4)  &
					,',',fte2_x(5,iff)/fte2max_x(5),',',fte2_x(6,iff)/fte2max_x(6)  &
					,',',fte2_x(7,iff)/fte2max_x(7),',',fte2_x(8,iff)/fte2max_x(8)  &
					,',',fte2_x(9,iff)/fte2max_x(9),',',fte2_x(10,iff)/fte2max_x(10)&
					,',',fte2_x(1,iff)/fte2_x(1,iff),',',fte2_x(2,iff)/fte2_x(1,iff)  &
					,',',fte2_x(3,iff)/fte2_x(1,iff),',',fte2_x(4,iff)/fte2_x(1,iff) &
					,',',fte2_x(5,iff)/fte2_x(1,iff),',',fte2_x(6,iff)/fte2_x(1,iff)  &
					,',',fte2_x(7,iff)/fte2_x(1,iff),',',fte2_x(8,iff)/fte2_x(1,iff) &
					,',',fte2_x(9,iff)/fte2_x(1,iff),',',fte2_x(10,iff)/fte2_x(1,iff) &
					,',',fte2_x(1,iff)/fte2_in(iff),',',fte2_x(2,iff)/fte2_in(iff)  &
					,',',fte2_x(3,iff)/fte2_in(iff),',',fte2_x(4,iff)/fte2_in(iff) &
					,',',fte2_x(5,iff)/fte2_in(iff),',',fte2_x(6,iff)/fte2_in(iff)  &
					,',',fte2_x(7,iff)/fte2_in(iff),',',fte2_x(8,iff)/fte2_in(iff)  &
					,',',fte2_x(9,iff)/fte2_in(iff),',',fte2_x(10,iff)/fte2_in(iff)
				write(89,871)  wavelength0,',',ffreq/tpi	&
					,',',fte2_y(1,iff)/fte2max_y(1),',',fte2_y(2,iff)/fte2max_y(2)     &
					,',',fte2_y(3,iff)/fte2max_y(3),',',fte2_y(4,iff)/fte2max_y(4)  &
					,',',fte2_y(5,iff)/fte2max_y(5),',',fte2_y(6,iff)/fte2max_y(6)  &
					,',',fte2_y(7,iff)/fte2max_y(7),',',fte2_y(8,iff)/fte2max_y(8)  &
					,',',fte2_y(9,iff)/fte2max_y(9),',',fte2_y(10,iff)/fte2max_y(10)  &
					,',',fte2_y(1,iff)/fte2_y(1,iff),',',fte2_y(2,iff)/fte2_y(1,iff)  &
					,',',fte2_y(3,iff)/fte2_y(1,iff),',',fte2_y(4,iff)/fte2_y(1,iff) &
					,',',fte2_y(5,iff)/fte2_y(1,iff),',',fte2_y(6,iff)/fte2_y(1,iff)  &
					,',',fte2_y(7,iff)/fte2_y(1,iff),',',fte2_y(8,iff)/fte2_y(1,iff) &
					,',',fte2_y(9,iff)/fte2_y(1,iff),',',fte2_y(10,iff)/fte2_y(1,iff) &
					,',',fte2_y(1,iff)/fte2_in(iff),',',fte2_y(2,iff)/fte2_in(iff)  &
					,',',fte2_y(3,iff)/fte2_in(iff),',',fte2_y(4,iff)/fte2_in(iff) &
					,',',fte2_y(5,iff)/fte2_in(iff),',',fte2_y(6,iff)/fte2_in(iff)  &
					,',',fte2_y(7,iff)/fte2_in(iff),',',fte2_y(8,iff)/fte2_in(iff)  &
					,',',fte2_y(9,iff)/fte2_in(iff),',',fte2_y(10,iff)/fte2_in(iff)
				endif
30309		continue

			if(TM .gt. 0)then
			close(87)
			endif
			if(TE .gt. 0)then
			close(88)
			close(89)
			endif

	endif !END OF SPECTRUM ANALYSIS (only performed when spectrumanalysis=1 & tempdat=1)

!DEC$ IF DEFINED (USE_MPI_CH) 
	endif 
!DEC$ ENDIF ! (USE_MPI_CH) 

871		format(e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3   &
                   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3  &
				   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3  &
				   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3  &
				   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3  &
				   ,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3,a3,e13.6e3  &
				   ,a3,e13.6e3)

!----------------------------------------------------------------------------------------
!			Integrattion of detector spectrum
!----------------------------------------------------------------------------------------
if(ftdet_num .ge. 1)then
		
	do 1238 k=1,det_num
		if((spec_index(k) .eq. 1) .and. (det_dir(k) .eq. 1))then
			do 1139 j=1,INT(det_ylen(k)/n_skip(k))
			do 1139 iff=0,ifmax
		
				ffreq=wstart + float(iff)*wspan/float(ifmax)
				ftarg=ffreq*float(n)*dt
				ftcos=cos(ftarg)
				ftsin=sin(ftarg)
					
				if(TM .gt. 0)then
    				fte2d_z(iff,j,k) = ftrd_z(iff,j,k)**2 + ftid_z(iff,j,k)**2 
				endif
				if(TE .gt. 0)then
					fte2d_x(iff,j,k) = ftrd_x(iff,j,k)**2 + ftid_x(iff,j,k)**2 
					fte2d_y(iff,j,k) = ftrd_y(iff,j,k)**2 + ftid_y(iff,j,k)**2 
				endif
				
1139		continue
		endif

		if((spec_index(k) .eq. 1) .and. (det_dir(k) .eq. 2))then
			do 1039 i=1,INT(det_xlen(k)/n_skip(k))
			do 1039 iff=0,ifmax
		
				ffreq=wstart + float(iff)*wspan/float(ifmax)
				ftarg=ffreq*float(n)*dt
				ftcos=cos(ftarg)
				ftsin=sin(ftarg)
					
				if(TM .gt. 0)then
    				fte2d_z(iff,i,k) = ftrd_z(iff,i,k)**2 + ftid_z(iff,i,k)**2
				endif
				if(TE .gt. 0)then
					fte2d_x(iff,i,k) = ftrd_x(iff,i,k)**2 + ftid_x(iff,i,k)**2
					fte2d_y(iff,i,k) = ftrd_y(iff,i,k)**2 + ftid_y(iff,i,k)**2
				endif
					
1039		continue
		endif

1238	continue

!write(*,*) "Finished squaring Fourier amplitude"

!DEC$ IF DEFINED (USE_MPI_CH)

	call MPI_BARRIER( comm_a, ierr ) ! synchronize entrance;
	
	do 1040 k=1,det_num
		
		if(TM.ge.1)then
			num_tot=(ifmax+1)*(det_maxlen+1)
			!num_tot=Size(fte2d_z)
			call MPI_Allreduce(fte2d_z(0,0,k),fte2d_zt(0,0,k),num_tot,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr ) 
		endif 
		if(TE.ge.1)then
			!num_tot=Size(fte2d_x)
			num_tot=(ifmax+1)*(det_maxlen+1)
			call MPI_Allreduce(fte2d_x(0,0,k),fte2d_xt(0,0,k),num_tot,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr ) 
			call MPI_Allreduce(fte2d_y(0,0,k),fte2d_yt(0,0,k),num_tot,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr ) 
		endif 
	
1040	continue

!write(*,*) "Finished MPI_allreduce "

	if( pid .eq. master ) then
 
!	integrate along the detector
		do 2003 k=1,det_num
		
			if((spec_index(k) .eq. 1) .and. (det_dir(k) .eq. 1))then
				do 2004 iff= 0, ifmax
					do 2004 j= 1, INT(det_ylen(k)/n_skip(k))
						if(TM.gt.0)then
						fte2d_zt(iff,0,k)=fte2d_zt(iff,0,k)+fte2d_zt(iff,j,k)
						endif
						if(TE.gt.0)then
						fte2d_xt(iff,0,k)=fte2d_xt(iff,0,k)+fte2d_xt(iff,j,k)
						fte2d_yt(iff,0,k)=fte2d_yt(iff,0,k)+fte2d_yt(iff,j,k)
						endif
2004			continue			
			endif
		
			if((spec_index(k) .eq. 1) .and. (det_dir(k) .eq. 2))then
				do 2005 iff= 0, ifmax
				do 2005 j= 1, INT(det_xlen(k)/n_skip(k))
					if(TM.gt.0)then
						fte2d_zt(iff,0,k)=fte2d_zt(iff,0,k)+fte2d_zt(iff,j,k)
					endif
					if(TE.gt.0)then
						fte2d_xt(iff,0,k)=fte2d_xt(iff,0,k)+fte2d_xt(iff,j,k)
						fte2d_yt(iff,0,k)=fte2d_yt(iff,0,k)+fte2d_yt(iff,j,k)
					endif
2005			continue			
			endif

2003	 continue

!write(*,*) "Finished integrating along the detector"

!maximum detection
		do 2006 k=1,det_num
			if(spec_index(k) .eq. 1)then
				do 2007 iff=0,ifmax
					if(TM.gt.0)then
						if(fte2d_zt(iff,0,k).gt.fte2maxd_z(k))then
							fte2maxd_z(k)=fte2d_zt(iff,0,k)
						endif
					endif
					if(TE.gt.0)then
						if(fte2d_xt(iff,0,k).gt.fte2maxd_x(k))then
							fte2maxd_x(k)=fte2d_xt(iff,0,k)
						endif
						if(fte2d_yt(iff,0,k).gt.fte2maxd_y(k))then
							fte2maxd_y(k)=fte2d_yt(iff,0,k)
						endif
					endif
2007			continue
			endif
2006	continue

!write(*,*) "Finished finding maximum"

!data output
		if(TM .gt. 0)then
			open(97, file= (trim(datadir)//'det_spectra_z.csv'))
		endif
		if(TE .gt. 0)then
			open(98, file= (trim(datadir)//'det_spectra_x.csv'))
			open(99, file= (trim(datadir)//'det_spectra_y.csv'))
		endif

		do 2008 iff=0,ifmax
			ffreq=wstart + float(iff)*wspan/float(ifmax)
       		wavelength0=tpi*cc/ffreq
			do 2009 k=1,det_num
				if(TM .gt. 0)then
					if (fte2maxd_z(k) .eq. 0)then
						fte2maxd_z(k)=1.
					endif
				endif
				if(TE .gt. 0)then
					if (fte2maxd_x(k) .eq. 0) Then
						fte2maxd_x(k)=1.
					endif
					if (fte2maxd_y(k) .eq. 0) Then
						fte2maxd_y(k)=1.
					endif
				endif
2009		continue
	
				!----output DFT spectrum with normalization and without
			if(TM .gt. 0)then 
				write(97,906) wavelength0,',',ffreq/tpi ! 
				write(97,907) ( ( ',', fte2d_zt(iff,0,k) ), k = 1,det_num )
				!write(97,907) ( ( ',', fte2d_zt(iff,0,k)/fte2maxd_z(k) ), k = 1,det_num )
				!write(97,907) ( ( ',', fte2d_zt(iff,0,k)/fte2d_zt(iff,0,1) ), k = 1,det_num )
				write(97,'(a1)') '' ! end if line
			endif
			if(TE .gt. 0)then
				write(98,906) wavelength0,',',ffreq/tpi ! 
				write(98,907) ( ( ',', fte2d_xt(iff,0,k) ), k = 1,det_num )
				!write(98,907) ( ( ',', fte2d_xt(iff,0,k)/fte2maxd_x(k) ), k = 1,det_num )
				!write(98,907) ( ( ',', fte2d_xt(iff,0,k)/fte2d_xt(iff,0,1) ), k = 1,det_num )
				write(98,'(a1)') '' ! end if line
				write(99,906) wavelength0,',',ffreq/tpi ! 
				write(99,907) ( ( ',', fte2d_yt(iff,0,k) ), k = 1,det_num )
				!write(99,907) ( ( ',', fte2d_yt(iff,0,k)/fte2maxd_y(k) ), k = 1,det_num )
				!write(99,907) ( ( ',', fte2d_yt(iff,0,k)/fte2d_yt(iff,0,1) ), k = 1,det_num )
				write(99,'(a1)') '' ! end if line
			endif

2008	continue

	endif ! (pid .eq. master)

!DEC$ ELSE ! (USE_MPI_CH)

!	integrate along the detector
		do 3003 k=1,det_num
		
			if((spec_index(k) .eq. 1) .and. (det_dir(k) .eq. 1))then
				do 3004 iff= 0, ifmax
					do 3004 j= 1, INT(det_ylen(k)/n_skip(k))
						if(TM.gt.0)then
						fte2d_z(iff,0,k)=fte2d_z(iff,0,k)+fte2d_z(iff,j,k)
						endif
						if(TE.gt.0)then
						fte2d_x(iff,0,k)=fte2d_x(iff,0,k)+fte2d_x(iff,j,k)
						fte2d_y(iff,0,k)=fte2d_y(iff,0,k)+fte2d_y(iff,j,k)
						endif
3004			continue			
			endif
		
			if((spec_index(k) .eq. 1) .and. (det_dir(k) .eq. 2))then
				do 3005 iff= 0, ifmax
				do 3005 j= 1, INT(det_xlen(k)/n_skip(k))
					if(TM.gt.0)then
						fte2d_z(iff,0,k)=fte2d_z(iff,0,k)+fte2d_z(iff,j,k)
					endif
					if(TE.gt.0)then
						fte2d_x(iff,0,k)=fte2d_x(iff,0,k)+fte2d_x(iff,j,k)
						fte2d_y(iff,0,k)=fte2d_y(iff,0,k)+fte2d_y(iff,j,k)
					endif
3005			continue			
			endif

3003	 continue

!maximum detection
		do 3006 k=1,det_num
			if(spec_index(k) .eq. 1)then
				do 3007 iff=0,ifmax
					if(TM.gt.0)then
						if(fte2d_zt(iff,0,k).gt.fte2maxd_z(k))then
							fte2maxd_z(k)=fte2d_z(iff,0,k)
						endif
					endif
					if(TE.gt.0)then
						if(fte2d_xt(iff,0,k).gt.fte2maxd_x(k))then
							fte2maxd_x(k)=fte2d_x(iff,0,k)
						endif
						if(fte2d_yt(iff,0,k).gt.fte2maxd_y(k))then
							fte2maxd_y(k)=fte2d_y(iff,0,k)
						endif
					endif
3007			continue
			endif
3006	continue

!data output

		if(TM .gt. 0)then
			open(97, file= (trim(datadir)//'det_spectra_z.csv'))
		endif
		if(TE .gt. 0)then
			open(98, file= (trim(datadir)//'det_spectra_x.csv'))
			open(99, file= (trim(datadir)//'det_spectra_y.csv'))
		endif
 
		do 3008 iff=0,ifmax
			ffreq=wstart + float(iff)*wspan/float(ifmax)
       		wavelength0=tpi*cc/ffreq
			do 3009 k=1,det_num
				if(TM .gt. 0)then
					if (fte2maxd_z(k) .eq. 0)then
						fte2maxd_z(k)=1.
					endif
				endif
				if(TE .gt. 0)then
					if (fte2maxd_x(k) .eq. 0) Then
						fte2maxd_x(k)=1.
					endif
					if (fte2maxd_y(k) .eq. 0) Then
						fte2maxd_y(k)=1.
					endif
				endif
3009		continue
	
				!----output DFT spectrum with normalization and without
			if(TM .gt. 0)then 
				write(97,906) wavelength0,',',ffreq/tpi ! 
				write(97,907) ( ( ',', fte2d_z(iff,0,k) ), k = 1,det_num )
				!write(97,907) ( ( ',', fte2d_z(iff,0,k)/fte2maxd_z(k) ), k = 1,det_num )
				!write(97,907) ( ( ',', fte2d_z(iff,0,k)/fte2d_z(iff,0,1) ), k = 1,det_num )
				write(97,'(a1)') '' ! end if line
			endif
			if(TE .gt. 0)then
				write(98,906) wavelength0,',',ffreq/tpi ! 
				write(98,907) ( ( ',', fte2d_x(iff,0,k) ), k = 1,det_num )
				!write(98,907) ( ( ',', fte2d_x(iff,0,k)/fte2maxd_x(k) ), k = 1,det_num )
				!write(98,907) ( ( ',', fte2d_x(iff,0,k)/fte2d_x(iff,0,1) ), k = 1,det_num )
				write(98,'(a1)') '' ! end if line
				write(99,906) wavelength0,',',ffreq/tpi ! 
				write(99,907) ( ( ',', fte2d_y(iff,0,k) ), k = 1,det_num )
				!write(99,907) ( ( ',', fte2d_y(iff,0,k)/fte2maxd_y(k) ), k = 1,det_num )
				!write(99,907) ( ( ',', fte2d_y(iff,0,k)/fte2d_y(iff,0,1) ), k = 1,det_num )
				write(99,'(a1)') '' ! end if line
			endif

3008	continue
!DEC$ ENDIF ! (USE_MPI_CH)

	if(TM .gt. 0)then
		close(97)
	endif
	if(TE .gt. 0)then
		close(98)
		close(99)
	endif

endif  !END OF INTEGRATED SPECTRUM DETECTOR

906 format(e13.6e3,a3,e13.6e3,\)
907 format(a3,e13.6e3,\)

!DEC$ IF DEFINED (USE_MPI_CH) 
	call MPI_FINALIZE(ierr) ! finalize MPICH interface; 
!DEC$ ENDIF ! (USE_MPI_CH) 

end  
