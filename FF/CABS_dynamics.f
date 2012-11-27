
C	***************************************************************
c
c	Ca-Cb-SG   lattice model of protein dynamics and thermodynamics
c	Edited March 2001 by Andrzej Kolinski (Kolinski@chem.uw.edu.pl)
c	Lattice representation of the Ca-trace employs 800 vectors type
c	v=(i,j,k) with |v*v| = 29,....        or 49. The scaling to the 
c	underlying cubic lattice is assumed to be 0.61 Angstroms to the
c	lattice unit(a small correction to this scalling may be needed).
c
C	***************************************************************
c
c
c	NOTES:   1. Derive a NN predictor for contact_order and contact
c	number (needs to be scaled uniformly due to different cut-offs)
c	2. Refine definition of the hydrophobic core elipsoid - may be
c	iportant for the refinement with a limited system's relaxation.
c	3. Improve sampling scheme.  4. Build a program for the initial
c	chains traces & the ditance collection from multiple templates. 
c	5. Implement "smart" large-scale moves via pdb-based rebuilding
c	procedure (local and global moves) - include occasional "image"
c	transformations at high temerature replicas. 6.Consider several
c	(~3) rotamers for the bigest hydrophobic residues (PHE,TYR,TRP)
c	7. Try the phantom chain approach at early stages of folding(?)
c	8. Try multiple local "smart" moves that optimize the Ca-Ca-Cb 
c	contacts before attempting the "soft" energy calculations, also
c	try the short range + H-bonds minimization's large scale moves.
c 

	IMPLICIT INTEGER(I-Z)
	LOGICAL LOOK, GOODC(800,800)
        character*5 struct(5)
	character*3 aa(-1:20), NAME

		PARAMETER(NDIM=400)
		PARAMETER(NREPS=50)
c
c	The recommended
c	number of replicas (T1=T2=1.0) is about 20 for small
c	proteins up to 30 for N>200
c
	
	COMMON /CHAINS/  ICA(NDIM), X(NDIM), Y(NDIM), Z(NDIM)
	COMMON /THREE/   ICONF(800,800)
	COMMON /LENGHTS/ LENF1,LENF,PROD(800,800)
	COMMON /SEQE/  	 SEQ(ndim), SEC(ndim) 	  
	COMMON /VECTORS/ vx(800),vy(800),vz(800)
	COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
	COMMON /BISEC/   CAX(800,800),CAY(800,800),CAZ(800,800)
	COMMON /HB/      HBX(800,800),HBY(800,800),HBZ(800,800)

	common/arandom/  aarand,abrand,acrand,adrand
	COMMON /ENERGY/	 EHBOND, ESC, EREP, BAT, EHBIJ(ndim,ndim)


	COMMON /pair/	apa(ndim,ndim,2,2),app(ndim,ndim,2,2),
     *			apm(ndim,ndim,2,2)
	COMMON /size/	arla(0:19,0:19,2,2),arlm(0:19,0:19,2,2),
     *			arlp(0:19,0:19,2,2)
	COMMON /sizea/	ala(0:19,0:19,2,2),alm(0:19,0:19,2,2),
     *			alp(0:19,0:19,2,2)
	COMMON /KBB/ kbin(800,800)

	COMMON /short/ IBIN(-500:500),asr(ndim,-12:12)
	COMMON /short1/ JBIN(800),bsr(ndim,16) 
	COMMON /short0/ SBIN(200),csr(ndim,8) 

	COMMON /one/ 	acrit, contt, eoinp(0:19,0:100)
	COMMON /shape/	amx,amy,amz
     	COMMON /FR/ FRG(ndim,ndim)
      COMMON /order/acorder,ardscal,SUMCT,SUMCTO,ICNT,ICNTO,dord,cstr
	COMMON /RES/ EREST, MRES(ndim),KRES(ndim,50),awrca(ndim,420)
	COMMON /RCN/ ERESTA, MRESA(ndim),KRESA(ndim,420),arca(ndim,420),
     *  	arcam(ndim,420)
	COMMON /short2/ erestsum, arestsum,arcs(ndim,ndim)
	DIMENSION csre(0:19,0:19,2)	
	COMMON /SG/ GX(800,800,0:19),GY(800,800,0:19),GZ(800,800,0:19)		
	DIMENSION irest(4000),jrest(4000),iresta(8000),jresta(8000)
	DIMENSION vector(-7:7,-7:7,-7:7)
	DIMENSION apabla(0:19,0:19,2,2),apablp(0:19,0:19,2,2)
	DIMENSION apablm(0:19,0:19,2,2)
	DIMENSION asrr(0:19,0:19,-12:12),bsrr(0:19,0:19,16)
	DIMENSION asrh(0:19,0:19,-12:12),bsrh(0:19,0:19,16)
	DIMENSION asre(0:19,0:19,-12:12),bsre(0:19,0:19,16)
	DIMENSION xt(ndim),yt(ndim),zt(ndim)	
	DIMENSION ATREP(NREPS),EREPLICA(NREPS), zreps(ndim,NREPS)
	DIMENSION xreps(ndim,NREPS),yreps(ndim,NREPS),SWAPS(NREPS)
	DIMENSION axalf(0:19),ayalf(0:19),azalf(0:19)
	DIMENSION axbet(0:19),aybet(0:19),azbet(0:19)
	DIMENSION apba(ndim,ndim),apbp(ndim,ndim),apbm(ndim,ndim)
      common/lim/erlim,arlim,erold,ernew,arold,arnew,actre,ectre 
	dimension fax(ndim),fay(ndim),faz(ndim) 
	common /icgg/ icg(ndim)
	dimension icao(ndim),icab(ndim)
	dimension xback(ndim),yback(ndim),zback(ndim)



c	
c	x,y,z - Ca lattice coordinates
c	[vx,vy,vz] 800 possible lattice Ca-Ca vectors 	
c	cax, cay, caz - normalyzed coordinates of (vi-vi+1)
c	cbx, cby, cbz - vectors from Ca to Cb (vi, vi+1 dependent)
c	gx,gy,gz - coordinates of the center of mas of sidechains
c	hbx,hby,hbz   - "hydrogen bond" vectors (vi, vi+1 dependent)
c	SEQ - sequence (LENF-2 numbers, 0=Gly, 1=Ala, ... 19=Trp)
c	SEC - secondary structure 1=coil, 2=helix, 3=turn, 4=beta
c	      (coil and turn treated in the same way)
c	apa,apm,app (pairwise, orientation and position dependent 
c	      potentials, apabla, apablp, apablm- sequence independent)
c	arla, arlm, arlp - cut-off distances for pairwise interactions
c	axalf, ayalf, azalf - coeficients of the side group coordinates
c             in local othonormal coordinate system - for alpha type
c	      main chain conformation (planar angle) 
c	axbet, aybet, azbet - as above for an open angle of the backbone

c	FRG - short range distances deduced from SEC - i.e. helix and
c	      beta target geometry
c	acrit - expected radius of gyration for a molecule


c
c
c	ESTIMATE THE "CONTACT ORDER" FROM THE TABLE   
c
c	  SIZE     H       E       H/E
c  	2   40  0.116   0.324      0.252
c   	3   60  0.119   0.357      0.230
c  	4   80  0.115   0.280      0.212
c   	5  100  0.105   0.259      0.198
c   	6  120  0.132   0.269      0.168
c   	7  140  0.105   0.272      0.176
c   	8  160  0.114   0.186      0.183
c   	9  180  0.116   0.197      0.160
c      10  200  0.107   0.184      0.134
c
c	TARGET NUMBER OF CONTACTS  1.7*N (when an atomic criterion used)
c
c	In this version use rather something like 1.9*N due to the wider
c	interaction well for the side chains.
c
c	INPUT
c	IJKM  10  KKK 20 	  (RANDOM, ICYCLE, PHOT,  N_OF_REPLICAS)
c	1.5  1.0  4.0  0.0 1.5    (T1,T2, REPULSIVE, BURIAL, ALPHAS_ATR)
c	0.75  1.75  1.5 -1.5 0.375(gen, PAIRS, CENTRO, H_bonds,  SHORTR)
c	0.75  0.0                 (PROFLE_BURIAL,SECSTRUC_PAIR_COUPLING)
c	XX.0  5.0  YYY.0  0.125   (Cont_order, scale, N-of_contacts, sc)
c	iiii, 0.5		  (#of_contacts,  scale)
c	jjjj, 0.5		  (#of_distances, scale)
c
c	Strenghts rescaled in this version

	data struct /' coil','helix',' turn',' beta','    ?'/
	data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
     &		     'PRO','MET','ASP','ASN','LEU',
     &		     'LYS','GLU','GLN','ARG',
     &		     'HIS','PHE','TYR','TRP','CYX'/

	
	OPEN(UNIT=7, FILE='SEQ',     STATUS='OLD')
	OPEN(UNIT=6, FILE='OUT',     STATUS='NEW')
	OPEN(UNIT=5, FILE='INP',     STATUS='OLD')
	OPEN(UNIT=10,FILE='FCHAINS',  STATUS='OLD')
	OPEN(UNIT=9, FILE='TRAF',   STATUS='NEW')
	OPEN(UNIT=26,FILE='QUASI3S',  STATUS='OLD')
	OPEN(UNIT=11,FILE='R14',   STATUS='OLD')
      OPEN(UNIT=23,FILE='ENERGY', STATUS='NEW')
	OPEN(UNIT=12,FILE='R15',   STATUS='OLD')
	OPEN(UNIT=13,FILE='R13',   STATUS='OLD')
	OPEN(UNIT=29,FILE='CENTRO',  STATUS='OLD')
	OPEN(UNIT=33,FILE='SIDECENT',STATUS='OLD')
	OPEN(UNIT=40, FILE='TRASG', STATUS='NEW')	
	OPEN(UNIT=66,FILE="TEMP_ENERGY",STATUS="NEW")
c
c	SEQ sequence   (i  ALA  SECONDARY_STRUCT_code=1,2,4)
c	INP input file with scaling parameters and restraints
c	ACHAINS - set of chains for Replicas, built by NatCA.f
c	QUASI3 - statistical pairwise potential, orient. dependent
c	R13, R14, R15 - short ramge statistical potentials
c	R14E, R14H, R15E, R15H are for known (predicted) fragments
c               of secondary structure, R14, R15 are generic
c
c	CENTRO - one-body centrosymmetric
c	PROFILE3 - Statistical, multibody, orientation dependent
c	SIDECENT - Coordinates of the side chains in respect to the
c		   backbone, diferent for Helices and Expanded
c	PAIR3 - "homology" based, protein dependent pairwise
c

	do k=0,19
	read(33,*)
	read(33,*) axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
	enddo
	CLOSE(33)		
c
c
c	Centrosymmetric one body potential, energy of aminoacid
c	i in the shell j of a spherical potential, 3-bins up to
c	Acrit contain about 55% of aminoacids, the rest in two
c	additional bins up to 5*acrit/3
	        do i=0,19
	       	read(29,3009) name, (eoinp(i,j), j=1,4)
 		eoinp(i,0)=eoinp(i,1)
 		eoinp(i,5)=eoinp(i,4) 		
               	do j=6,100
               	eoinp(i,j)=eoinp(i,j-1)+0.25
c	for larger distances energy grows in alinear fashion
               	enddo
                enddo
3009 		FORMAT(A3,4F6.2)
	CLOSE(29)
	
	
	READ(5,*) RANDOM, NCYCLE, PHOT, REPLICAS
	READ(5,*) ATEMP2, ATEMP1, EREP, EKD, BAT
	READ(5,*) ESC_0, ARLO, ENONE, EHBOND, ASRO
	READ(5,*) aprofs
	READ(5,*) ACORDER, ARDSCAL, contt, cstr

	adrand=121.3421721
	abrand=0.3709112
	acrand=0.8021453
	aarand=float(random)/(12.3211+1.5321*float(RANDOM))
	
	write(6,*) arand(ise)

	REWIND(11)

					
	do i=0,19
	do j=0,19
	read(11,*)
c	reading of "average" r14 interactions, it is chiral
	read(11,*) (asrr(i,j,k),k=-12,-5)
	read(11,*) (asrr(i,j,k),k=-4,3)
	read(11,*) (asrr(i,j,k),k=5,12)
	do k=4,1,-1
	asrr(i,j,k)=asrr(i,j,k-1)
	enddo
	enddo
	enddo
 	
	do i=0,19
	do j=0,19
	read(12,*)
c	reading of "average" r15 interactions, k-type of aminoacid
	read(12,*) (bsrr(i,j,k),k=1,16)
c	read(12,*) (bsrr(i,j,k),k=1,8)
c	read(12,*) (bsrr(i,j,k),k=9,16)
	enddo
	enddo
	
	CLOSE(11)
	CLOSE(12)

	OPEN(UNIT=11,FILE='R14H',   STATUS='OLD')
	OPEN(UNIT=12,FILE='R15H',   STATUS='OLD')
	REWIND(11)
	REWIND(12)
	

c
c	When secondary structure known (SEC=2, or SEC=4) use
c	statistical potential thatis derived for specific secondary
c	structure elements H, or E  (read in two following segments)
	do i=0,19
	do j=0,19
	read(11,*)
	read(11,*) (asrh(i,j,k),k=-12,-5)
	read(11,*) (asrh(i,j,k),k=-4,3)
	read(11,*) (asrh(i,j,k),k=5,12)
	do k=4,1,-1
	asrh(i,j,k)=asrh(i,j,k-1)
	enddo
	enddo
	enddo
 	
	do i=0,19
	do j=0,19
	read(12,*)
	read(12,*) (bsrh(i,j,k),k=1,16)	
c	read(12,*) (bsrh(i,j,k),k=1,8)
c	read(12,*) (bsrh(i,j,k),k=9,16)
	enddo
	enddo	

	CLOSE(11)
	CLOSE(12)
	
	OPEN(UNIT=11,FILE='R14E',   STATUS='OLD')
	OPEN(UNIT=12,FILE='R15E',   STATUS='OLD')

	REWIND(11)
	REWIND(12)
	
	
	do i=0,19
	do j=0,19
	read(11,*)
	read(11,*) (asre(i,j,k),k=-12,-5)
	read(11,*) (asre(i,j,k),k=-4,3)
	read(11,*) (asre(i,j,k),k=5,12)
	do k=4,1,-1
	asre(i,j,k)=asre(i,j,k-1)
	enddo
	enddo
	enddo
 	
	do i=0,19
	do j=0,19
	read(12,*)
	read(12,*) (bsre(i,j,k),k=1,16)
c	read(12,*) (bsre(i,j,k),k=1,8)
c	read(12,*) (bsre(i,j,k),k=9,16)
	enddo
	enddo	

	CLOSE(11)
	CLOSE(12)
		
	do i=1,800
	kk=int((sqrt(float(i))*0.61))+1
	if(kk.gt.16) kk=16
	JBIN(I) = kk
	ENDDO	
	
		
	do i=1,500
	kk=int((sqrt(float(i))*0.61))+1
	if(kk.gt.12) kk=12
	IBIN(I) = kk
	IBIN(-I)=-kk
	ENDDO
	IBIN(0)=IBIN(1)


	do ii=0,19

		do j=0,100
		eoinp(ii,j)=enone*eoinp(ii,j)
		enddo
	enddo

c
c	Pairwise interactions apablp ... and cut-off parmeters
c 	arlp, orientation dependent, pairwise specific, sequence
c	independent  (two types of local conformations
c
c	Compact(H)  r13<6A
c	EXPANDED	r13>6A
c	Potentials  HH, EE, HE, EH
c
c
 	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apablp(i,j,1,1),j=0,19)
	enddo
 	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apablp(i,j,2,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apablp(i,j,1,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apablp(i,j,2,1),j=0,19)
	enddo



	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apablm(i,j,1,1),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apablm(i,j,2,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apablm(i,j,1,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apablm(i,j,2,1),j=0,19)
	enddo



	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apabla(i,j,1,1),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apabla(i,j,2,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apabla(i,j,1,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (apabla(i,j,2,1),j=0,19)
	enddo



	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arlp(i,j,1,1),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arlp(i,j,2,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arlp(i,j,1,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arlp(i,j,2,1),j=0,19)
	enddo


	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arlm(i,j,1,1),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arlm(i,j,2,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arlm(i,j,1,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arlm(i,j,2,1),j=0,19)
	enddo


	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arla(i,j,1,1),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arla(i,j,2,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arla(i,j,1,2),j=0,19)
	enddo
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,725) NAME, (arla(i,j,2,1),j=0,19)
	enddo

725	format(a3,1x,20f5.1)



c
C	***************************************************************
C	 PREPARATION OF THE CHAIN UNITS AND THEIR CORRELATION



	NWMAX=0
	DO ix=-7,7
	DO iy=-7,7
	DO iz=-7,7
	vector(ix,iy,iz)=0		
	ir=ix*ix+iy*iy+iz*iz
	if(ir.ge.29) then
	if(ir.le.49) then
c	definition of Ca-Ca lattice vectors
	NWMAX=NWMAX+1
	VX(NWMAX)=ix
	VY(NWMAX)=iy
	VZ(NWMAX)=iz
	VECTOR(ix,iy,iz)=NWMAX
c	VECTOR()= 1,2,3....800
c	ix,iy,iz= -7, -6, -4,-3,  ...3,4,6,7 with  28< ||vx,vy,vz||<50
c	this is an arbitrary choice of the lattice geometry
	endif
	endif
	ENDDO
	ENDDO
	ENDDO

	write(*,*) NWMAX

c
c	READ THE CHAIN -DEFINE VECTOR INDICES
c	
	
	READ(10,*) LENF
c	LENF -chain length = # of AA + two dummy end caps
	LENF1=LENF-1
	LENF2=LENF-2
	LENF3=lenf-3
	AL2=LENF2-0.0001
	AL4=LENF-4.0001
	AL5=LENF-5.0001
	al19=lenf-19.0001

	if(al19.lt.2) al19=2
		
c
c	*******************   CREATING REPLICAS  ********************
c

	REWIND(10)
	do k=1,REPLICAS
	SWAPS(k)=0
	read(10,*)
	DO I=1,LENF
	READ(10,*) X(I),Y(I),Z(I)
	ENDDO

	do i=1,lenf
c	replicas Ca coordinates
	xreps(i,k)=x(i)
	yreps(i,k)=y(i)
	zreps(i,k)=z(i)
	enddo
	enddo

		ICA(LENF)=1
		
c
c	SEQUENCE INPUT
c
c	
	 			do 121 i=2,lenf1
 				read(7,707) k, NAME,SEC(I)
				do j=0,20
 				if(NAME.eq.aa(j)) then
 				SEQ(i)=j
c	sequence SEQ()=0 means GLY,  =19 means TRP
	icg(i)=0
	if(NAME.eq.'ASP'.or.NAME.eq.'GLU') then
		NMINUS=NMINUS+1
		icg(i)=-1
		endif
	if(NAME.eq.'LYS'.or.NAME.eq.'ARG') then
		NPLUS=NPLUS+1
		icg(i)= 1
		endif	
 				go to 121
 				endif
				enddo
 121				continue
 707				format(i5,3x,a3,2i5)
 


 	
c
c	EHBIJ - set-up secondary structure dependent
c	strength of the hyrogen bond network - stroger for helices
c 	and beta-sheets
c
	EHSEC=EHBOND*1.5	
	do i=2,lenf1
	is=sec(i)
	do j=2,lenf1
	if(iabs(i-j).gt.2) then
	EHBIJ(i,j)=0.75*EHBOND
	js=sec(j)

	if(is.ne.2.AND.is*js.eq.4) EHBIJ(i,j)=EHBOND	
	if(iabs(i-j).eq.3.and.is.eq.2.and.js.eq.2) EHBIJ(i,j)=EHSEC
	
	if(is.eq.4.and.js.eq.4) EHBIJ(i,j)=EHSEC
	if(is*js.EQ.8) EHBIJ(i,j)=0
	else
	EHBIJ(i,j)=0.0
	endif
	
	enddo
	enddo

c
c	Combine the generic (asrr) short range potential with
c 	secondary structure dependent potentials (asrh and asre)
c	R14 chiral potential
c
 	do i=1,lenf-3
	do k=-12,12
	asr(i,k)=ASRO*asrr(seq(i+1),seq(i+2),k)
	if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
	if(sec(i).eq.2) then
	asr(i,k)=ASRO*asrh(seq(i+1),seq(i+2),k)
	endif
	endif
	if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
	if(sec(i).eq.4) then
	asr(i,k)=ASRO*asre(seq(i+1),seq(i+2),k)
	endif
	endif
	enddo
	enddo


c
c	R15 potential - strength (ASRO) reduced as above
c	
 	do i=1,lenf-4
	do k=1,16
	bsr(i,k)=ASRO*bsrr(seq(i+1),seq(i+3),k)
	if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
	if(sec(i).eq.2.AND.sec(i+4).eq.2) then
	bsr(i,k)=ASRO*bsrh(seq(i+1),seq(i+3),k)
	endif
	endif
	if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
	if(sec(i).eq.4.AND.sec(i+4).EQ.4) then
	bsr(i,k)=ASRO*bsre(seq(i+1),seq(i+3),k)
	endif
	endif
	enddo
	enddo


c	R13 potential - strength ASRO
c			
	REWIND(13)
	do i=1,200
	kk=int((sqrt(float(i))*0.61))+1
	if(kk.gt.8) kk=8
	SBIN(I) = kk
	ENDDO


	REWIND(13)
	do i=0,19
	do j=0,19
	read(13,*)
c	reading of "average" r13 interactions, k-type of aminoacid
	read(13,*) (bsrr(i,j,k),k=1,8)
	enddo
	enddo
	
	CLOSE(13)
	OPEN(UNIT=13,FILE='R13H',   STATUS='OLD')
	REWIND(13)	
			
	do i=0,19
	do j=0,19
	read(13,*)
c	reading of "average" r13 interactions, k-type of aminoacid
	read(13,*) (bsrh(i,j,k),k=1,8)
	enddo
	enddo
	
	CLOSE(13)
	OPEN(UNIT=13,FILE='R13E',   STATUS='OLD')
	REWIND(13)	

	do i=0,19
	do j=0,19
	read(13,*)
c	reading of "average" r13 interactions, k-type of aminoacid
	read(13,*) (bsre(i,j,k),k=1,8)
	enddo
	enddo
	
	CLOSE(13)
	
 	do i=2,lenf-3
	do k=1,8
	csr(i,k)=ASRO*bsrr(seq(i),seq(i+2),k)
	if(sec(i).eq.2.AND.sec(i+1).eq.2.AND.sec(i+2).eq.2)then
	csr(i,k)=(csr(i,k)+ASRO*bsrh(seq(i),seq(i+2),k))/2.0
	endif
	if(sec(i).eq.4.AND.sec(i+1).eq.4.AND.sec(i+2).eq.4)then
	csr(i,k)=(csr(i,k)+ASRO*bsre(seq(i),seq(i+2),k))/2.0
	endif
	enddo
	enddo	


	do k=1,8
	csr(1,k)=0.0
	csr(lenf-2,k)=0.0
	enddo

	close(26)


	do i=2,lenf1
	ii=SEQ(i)
	do j=2,lenf1
	jj=SEQ(j)

	dd=1.0


c	
	if(iabs(i-j).eq.6) dd=0.50
	if(iabs(i-j).eq.4) dd=0.0

	do k1=1,2
	do k2=1,2			
	apa(i,j,k1,k2)=arlo*(apabla(ii,jj,k1,k2)-dd)
	apm(i,j,k1,k2)=arlo*(apablm(ii,jj,k1,k2)-dd)
	app(i,j,k1,k2)=arlo*(apablp(ii,jj,k1,k2)-dd)
	enddo
	enddo

	enddo
	enddo

	
		d1=1.5 
		d2=2.0  

c	d1+d2 the width of the square well pairwise potential

	do k1=1,2
	do k2=1,2	
	do i=0,19		
	do j=0,19
	d11=d1
	d22=d2


	ala(i,j,k1,k2)=((arla(i,j,k1,k2)+d11)/0.61)**2
	alm(i,j,k1,k2)=((arlm(i,j,k1,k2)+d11)/0.61)**2
	alp(i,j,k1,k2)=((arlp(i,j,k1,k2)+d11)/0.61)**2
      if(apabla(i,j,k1,k2).gt.0)ala(i,j,k1,k2)=(arla(i,j,k1,k2)/0.61)**2
      if(apablm(i,j,k1,k2).gt.0)alm(i,j,k1,k2)=(arlm(i,j,k1,k2)/0.61)**2
      if(apablp(i,j,k1,k2).gt.0)alp(i,j,k1,k2)=(arlp(i,j,k1,k2)/0.61)**2

c
	arla(i,j,k1,k2)=(max(2.,(arla(i,j,k1,k2)-d22))/0.61)**2
	arlm(i,j,k1,k2)=(max(2.,(arlm(i,j,k1,k2)-d22))/0.61)**2
	arlp(i,j,k1,k2)=(max(2.,(arlp(i,j,k1,k2)-d22))/0.61)**2

 	enddo	
	enddo
	enddo
	enddo

								
c
c	READS Side group - side group contacts (from NMR
c	or from therading predictions)
c
	do i=1,lenf
	MRES(i)=0
	MRESA(i)=0
	enddo
	
	READ(5,*) NREST, EREST

	WRITE(6,*)
	WRITE(6,*) 'STRENGTH & # of SG-SG RESTRAINTS:  ',EREST,NREST
	WRITE(6,*)
	if(NREST.gt.0 ) THEN
	do k=1, NREST
	read(5,*) IREST(K),JREST(k),junki
	i=irest(k)+1
	j=jrest(k)+1
	ii=MRES(i)+1
	jj=MRES(j)+1
	MRES(i)=ii
	MRES(j)=jj
c	KRES(i,ii) ii contacts of residue i
	KRES(i,ii)=j
	KRES(j,jj)=i
	irest(k)=i
	jrest(k)=j
	enddo
	endif

c	
c	RESCALE THE STRENGTH OF THE CONTACT RESTRAINTS
c
c	READS THE SHORT RANGE RESTRAINTS - FRAGMENTS
c
	nresta=0
	READ(5,*) NRESTAA, ERESTA
c	
c
	IF(NRESTAA.gt.0 ) THEN
	do k=1, NRESTAA
	read(5,*) IRESTA(K),JRESTA(k), junki, arrr,arrrm, awr
	i=iresta(k)+1
	j=jresta(k)+1
	NRESTA=NRESTA+1
	IRESTA(NRESTA)=i
	JRESTA(NRESTA)=j
	arrr=arrr/0.61
	arrrm=arrrm/0.61
	

c	arca(i,j) distance between i-th and j-th Ca

	ii=MRESA(i)+1
	jj=MRESA(j)+1
	MRESA(i)=ii
	MRESA(j)=jj

	arca(i,ii)=arrr
	arca(j,jj)=arrr
	arcam(j,jj)=arrrm
	arcam(i,ii)=arrrm
	awrca(i,ii)=awr
	awrca(j,jj)=awr


	KRESA(i,ii)=j
	KRESA(j,jj)=i	

	enddo
	ENDIF
c
c	NORMALIZE THE STRENGTH OF THE FLEXIBLE RESTRAINTS ???????
c

	erlim=0.0
	arlim=0.0

c	for a weaker homology use the set below
	erlim=0.5*float(nrest)
	arlim=0.5*float(nresta)
	
 	sec(1)=1
 	sec(lenf)=1


c	Set a generalized shell model- Acrit is radius of gyration
c	the equation is semiempirical regresion against chain length AL2
c	0.61 is the scaling to the lattice

	Acrit= 2.2*exp(0.38*alog(AL2))/0.61
	ACRIT2=Acrit*Acrit


c
c	800 vectors representing possible representations of the
c	alpha carbon virtual bonds.  1 lattice units = 0.61 A
c

c	===============================================================	


c
c	Prepare the set of good pairs of vectors - exclude too narrow
c	pairs and the colinear ones - to enable Cb definition
c
	DO I=1,800
	ix=vx(i)
	iy=vy(i)
	iz=vz(i)

	DO J=1,800
	goodc(i,j)=.TRUE.
	jx=vx(j)
	jy=vy(j)
	jz=vz(j)
	prod(i,j)= ix*jx+iy*jy+iz*jz
	ICONF(i,j)=(ix+jx)**2+(iy+jy)**2+(iz+jz)**2	
	if(iconf(i,j).lt.45) GOODC(i,j)=.FALSE.	
	if(iconf(i,j).gt.145) GOODC(i,j)=.FALSE.	
	
	kx=iy*jz-iz*jy
	ky=jx*iz-ix*jz
	kz=ix*jy-iy*jx
	product=(kx*kx+ky*ky+kz*kz)
	IF(product.EQ.0) goodc(I,J)=.FALSE.			
c
c	Prepare Cb positions - this is temporary (USE A FILE INSTEAD)
c	ideal tetrahedral conformation of Ca assumed -  should be OK

	IF(GOODC(i,j)) THEN

	if(iconf(i,j).lt.95) then
	kbin(i,j)=1
	else
	kbin(i,j)=2
	endif

	cx=ix+jx
	cy=iy+jy
	cz=iz+jz
	a=sqrt(cx*cx+cy*cy+cz*cz)
	cx=cx/a
	cy=cy/a
	cz=cz/a
	
	ax=ix-jx
	ay=iy-jy
	az=iz-jz
	
	a=sqrt(ax*ax+ay*ay+az*az)
	ax=ax/a
	ay=ay/a
	az=az/a
	
	CAX(i,j)=ax
	CAY(i,j)=ay
	CAZ(i,j)=az
c
c	Here define also the H-bond (Ca-based vectors)- 4.5A
c	
	b=sqrt(float(product))
	bx=float(kx)/b
	by=float(ky)/b
	bz=float(kz)/b
	
	HBX(i,j)=8.25*bx
	HBY(i,j)=8.25*by	
	HBZ(i,j)=8.25*bz
	
	DO k=0,19
	IF(iconf(i,j).LT.75) THEN
c	write down the helical Side Group positions
c	for a helical and turn-like conformations 	
	GX(i,j,k)=(axalf(k)*cx+ayalf(k)*bx+azalf(k)*ax)/0.61
	GY(i,j,k)=(axalf(k)*cy+ayalf(k)*by+azalf(k)*ay)/0.61
	GZ(i,j,k)=(axalf(k)*cz+ayalf(k)*bz+azalf(k)*az)/0.61
	ELSE
c	for expanded conformations
	IF(iconf(i,j).gt.110) THEN
	GX(i,j,k)=(axbet(k)*cx+aybet(k)*bx+azbet(k)*ax)/0.61
	GY(i,j,k)=(axbet(k)*cy+aybet(k)*by+azbet(k)*ay)/0.61
	GZ(i,j,k)=(axbet(k)*cz+aybet(k)*bz+azbet(k)*az)/0.61
	else
c	an approximation of the orientational averaging-consider fixing it 
c	(75-110)  is the intermediate range of iconf
	fb=float(iconf(i,j)-75)/35.0
	fa=1-fb
	axb=axbet(k)*fb+fa*axalf(k)
	ayb=aybet(k)*fb+fa*ayalf(k)
	azb=azbet(k)*fb+fa*azalf(k)
	GX(i,j,k)=(axb*cx+ayb*bx+azb*ax)/0.61
	GY(i,j,k)=(axb*cy+ayb*by+azb*ay)/0.61
	GZ(i,j,k)=(axb*cz+ayb*bz+azb*az)/0.61
	endif
	ENDIF

	ENDDO
	CBX(i,j)=GX(i,j,1)
	CBY(i,j)=GY(i,j,1)
	CBZ(i,j)=GZ(i,j,1)
	ENDIF	
	ENDDO

	ENDDO

	write(6,*) contt, cstr, kkkk
	
c
c	COMPUTE THE SECONDARY FRAGMENT BIASES
c
	do i=1,lenf
	do j=1,lenf
	FRG(j,i)=0.0
	enddo
	enddo

	do i=2,lenf-7
	  q=0
	  do j=i+1,i+5
	  if(sec(j).eq.4) q=q+4
	  enddo
c	length of six bond beta-type fragment
	  if(q.eq.20) FRG(i,i+6)=19.1/0.61
	enddo

	do i=2,lenf-8
	  q=0
	  do j=i,i+7
	  if(sec(j).eq.2) q=q+2
	  enddo
c	length of eight bond Helix-type fragment
	  if(q.eq.16) FRG(i,i+7)=10.75/0.61
	enddo	
	

c
c	THE MOST EXTERNAL LOOP FOR THE MC DYNAMICS
c	twenty temperatures for the annealing - setting 
c	atemp2=atemp1 reduces program to a more standard REMC
c	dtem1 -annealing increment
	dtem1=(atemp2-atemp1)/20.
	
	dtem=5.8/float(lenf)
	dtem=0.05
c	dtem=0.1
c	dtem distance between replicas	
	aatemp2=atemp2

	DO IDDDUM=1,20
	inreps=0
		
	aatemp1=aatemp2-dtem1	

	atemp=aatemp1

	do i=1, REPLICAS
	
	ATREP(i)=atemp
C	atemp=atemp+dtem*((float(i-1+REPLICAS))/float(REPLICAS))**3
	atemp=atemp+dtem
c	atemp temperature of a replica, lowest replica temp=aatemp1
	EREPLICA(i)=0.0
	enddo
	
	aatemp2=aatemp1
	
	I2flip=0
	I3flip=0
	I6flip=0
	I8flip=0


	
	DO 7000 ICYCLE=1,NCYCLE

	IF(ICYCLE.GT.1) THEN
	
	If(REPLICAS.GT.1) THEN
    
	DO ITEMP=REPLICAS-1,1,-1
        	WRITE(66,666) ATREP(ITEMP),EREPLICA(ITEMP)
	ENDDO
666	FORMAT(2F13.3)
	DO ITEMP=REPLICAS-1,1,-1
		DO ITEMP1=REPLICAS, ITEMP+1, -1
c
c	exchanging replicas
c
	betam=1.0/ATREP(itemp)
	betan=1.0/ATREP(itemp1)
	delta=(betan-betam)*(EREPLICA(itemp)-EREPLICA(itemp1))

			if(exp(-delta).gt.arand(ise)) then
			DO i=1,lenf
			xt(i)=xreps(i,itemp)
			yt(i)=yreps(i,itemp)
			zt(i)=zreps(i,itemp)
			xreps(i,itemp)=xreps(i,itemp1)
			yreps(i,itemp)=yreps(i,itemp1)
			zreps(i,itemp)=zreps(i,itemp1)
			xreps(i,itemp1)=xt(i)
			yreps(i,itemp1)=yt(i)
			zreps(i,itemp1)=zt(i)
			enddo
			SWAPS(itemp)=SWAPS(itemp)+1
			SWAPS(itemp1)=SWAPS(itemp1)+1			
			attt=EREPLICA(itemp)
			EREPLICA(itemp)=EREPLICA(itemp1)
			EREPLICA(itemp1)=attt
			inreps=inreps+1	
			endif
		ENDDO
	ENDDO
	ENDIF
	ENDIF

c
c	******************* ITERATE THE REPLICAS ********************
c
	

	
	DO ITEMP=REPLICAS, 1, -1
	

			DO i=1,lenf
			x(i)=xreps(i,itemp)
			y(i)=yreps(i,itemp)
			z(i)=zreps(i,itemp)
			enddo					
	
		DO I=1,LENF1
		J=I+1
		WX=X(J)-X(I)
		WY=Y(J)-Y(I)
		WZ=Z(J)-Z(I)
		ICA(I)=VECTOR(WX,WY,WZ)
c
c	Define chain of vectors ICA() =1,2,...800
c
		ENDDO

c
c	?????? RESCALING OF THE GENERIC TERMS  ????????????????????????
c	
	ATEMP=ATREP(itemp)
	
c
c
c	SCALING UP AT HIGH TEMP THE SHORT RANGE INTERACTIONS (GENERIC)
c
c	ESC=ESC_0 *sqrt(atemp/aatemp1)*(aatemp1/ATEMP1)
	ESC=ESC_0 *sqrt(aatemp1/ATEMP1)






	ESC=ESC_0



c	a new cycle starts here (internal "time" loop)


	DO 7777 IPHOT=1, PHOT

c
c	CENTER THE CHAIN at [0,0,0] point of the lattice
c
	sx=0
	sy=0
	sz=0
	do i=1,lenf
	sx=sx+x(i)
	sy=sy+y(i)
	sz=sz+z(i)
	enddo
	sx=nint(sx/float(lenf))
	sy=nint(sy/float(lenf))
	sz=nint(sz/float(lenf))
	do i=1,lenf
	x(i)=x(i)-sx
	y(i)=y(i)-sy
	z(i)=z(i)-sz
	enddo

	asx=0
	asy=0
	asz=0
	do i=2, lenf1
           	k=seq(i)
        		ii=ICA(i-1)
       		jj=ICA(i)
        	asx=asx+(x(i)+GX(ii,jj,k))**2
        	asy=asy+(y(i)+GY(ii,jj,k))**2
        	asz=asz+(z(i)+GZ(ii,jj,k))**2
	enddo
c
c	amx, amy, amz dimenions of an elipsoid containing
c	(presumably) the protein hyrophobic core
c	
	amx=(asx)/al2
	amy=(asy)/al2
	amz=(asz)/al2

        ICNTO=0
        SUMCTO=0 		
        ICNT=0
        SUMCT=0 		

	erestsum=0.0
	arestsum=0.0
	erold=0.0
	arold=0.0		

	energ=EHB(2,lenf1,1)+ESHORT(2,lenf1,10)

	erestsum=ernew
	arestsum=arnew
c
c	Compute energy from "contact order" (acorder)
c	and the expected total number of contacts (contt)
c		
	if(sumct.gt.0) energ=energ+ardscal*abs(ICNT/float(SUMCT)-acorder) 
     *	+cstr*abs(float(SUMCT)-contt)
     
        ICNTO=ICNT
        SUMCTO=SUMCT 	

 	DO 7775 IDUM=1,LENF2
c
c	TWO BOND KINKS- THE SMALLEST MOVES ATTEMPTED 10-times 
c	more frequently
c
	DO 7774 IONEBEAD=1,10
	
	j=INT(arand(ise)*AL4)+3
	i=j-1
	ii=ica(i)
	jj=ica(j)
	k=J+1


 3307	ix=int(arand(ise)*6.9999)-3 +x(j)
 	iy=int(arand(ise)*6.9999)-3 +y(j)
 	iz=int(arand(ise)*6.9999)-3 +z(j)
 	ir=(x(i)-ix)**2+(y(i)-iy)**2+(z(i)-iz)**2
 	if(ir.gt.49) go to 3307
 	if(ir.lt.29) go to 3307
 	ir=(x(k)-ix)**2+(y(k)-iy)**2+(z(k)-iz)**2
 	if(ir.gt.49) go to 3307	
  	if(ir.lt.29) go to 3307
  	
	nv1=vector((-x(i)+ix),(-y(i)+iy),(-z(i)+iz))
	if(nv1.eq.ii) go to 7774
	nv2=vector((x(k)-ix),(y(k)-iy),(z(k)-iz))
	if(.NOT.GOODC(nv1,nv2)) GO to 3307
	
c
c	nv1, nv2 - possible new values of ICA(j-1),ICA(j)
c	j- is the moved Ca unit
	
	ica1=ica(i-1)
c	check planar angles with the rest of the chain
	if(.NOT.GOODC(ica1,nv1)) go to 3307
	ica3=ica(k)
	if(.NOT.GOODC(nv2,ica3)) go to 3307
	
			px=X(i)+vx(nv1)
			py=Y(i)+vy(nv1)
			pz=Z(i)+vz(nv1)
	
	kx=x(j)
	ky=y(j)
	kz=z(j)			
	ICA(I)=nv1
	ICA(j)=nv2	
	x(j)=px
	y(j)=py
	z(j)=pz
		
	IF(LOOK(i,k)) THEN
c	No steric overlaps of Ca and Cb detected for the new conformer
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE
c	

	
		ENEW=EHB(i,k,1)+ESHORT(i,k,1)

		ICA(I)=ii
		ICA(j)=jj	
		x(j)=kx
		y(j)=ky
		z(j)=kz
	
		EOLD=EHB(i,k,-1)+ESHORT(i,k,-1)
						
		DE=ENEW-EOLD+dord
		IF(DE.GT.0.0) THEN
		if(arand(ise).gt.EXP(-DE/ATEMP)) then
    		ICNT=ICNTO
        	SUMCT=SUMCTO 		
		go to 7774
		endif
		ENDIF
c
c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION
c
        ICNTO=ICNT
        SUMCTO=SUMCT 	

		ENERG=ENERG+DE
		I2flip=I2flip+1
		ICA(I)=nv1
		ICA(j)=nv2	
		x(j)=px
		y(j)=py
		z(j)=pz	
	erestsum=ectre
	arestsum=actre	
	
		go to 7774
	ELSE
c
c	Move rejected - restore the initial state of the chain		
	ICA(I)=ii
	ICA(j)=jj	
	x(j)=kx
	y(j)=ky
	z(j)=kz
	ENDIF

 7774	CONTINUE	
c
c	THREE-BOND MOVE -permutation
c	

 	I=INT(arand(ise)*AL5)+2
	a=arand(ise)
	
	ii=ica(i)
	j=i+1
	jj=ica(j)
	k=i+2
	kk=ica(k)
	l=i+3
	
	if(arand(ise).gt.0.5) then
	nv1=kk	
	nv2=jj	 
	nv3=ii
	if(nv1.eq.ii) go to 8885
	else
		if(arand(ise).gt.0.5) then
		nv1=kk	
		nv2=ii	 
		nv3=jj
		if(nv1.eq.ii) go to 8885
		else
		nv1=jj	
		nv2=kk	 
		nv3=ii
		if(nv1.eq.ii) go to 8885
		endif	
	endif
			
	ica1=ica(i-1)
	if(.NOT.GOODC(ica1,nv1)) go to 8885
	if(.NOT.GOODC(nv1,nv2)) go to 8885
	if(.NOT.GOODC(nv2,nv3)) go to 8885	
	if(.NOT.GOODC(nv3,ica(l))) go to 8885		

c	Store the old conformation
c
	jx=x(j)
	jy=y(j)
	jz=z(j)
	
	kx=x(k)
	ky=y(k)
	kz=z(k)
	
	NX=x(i)+vx(nv1)
	NY=y(i)+vy(nv1)
	NZ=z(i)+vz(nv1)
		MX=NX+vx(nv2)
		MY=NY+vy(nv2)
		MZ=NZ+vz(nv2)
	X(j)=NX
	Y(j)=NY
	Z(J)=NZ
	X(k)=MX
	Y(k)=MY
	Z(k)=MZ
		ICA(i)=NV1
		ICA(j)=NV2
		ICA(k)=NV3
		
	IF(LOOK(i,l)) THEN
c	No overalps of Cas and Cbs
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE
c
		
		ENEW=EHB(i,l,1)+ESHORT(i,l,1)

		ICA(i)=ii
		ICA(j)=jj
		ICA(k)=kk
		x(j)=jx
		y(j)=jy
		z(j)=jz
		x(k)=kx
		y(k)=ky
		z(k)=kz

		EOLD=EHB(i,l,-1)+ESHORT(i,l,-1)
		
		DE=ENEW-EOLD+dord	

		IF(DE.GT.0.0) THEN
		if(arand(ise).gt.EXP(-DE/ATEMP)) then
 	   	ICNT=ICNTO
        	SUMCT=SUMCTO 		
		go to 8885
		endif
		ENDIF
		
		
c
c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION
c
        ICNTO=ICNT
        SUMCTO=SUMCT 	


		ENERG=ENERG+DE	
		X(j)=NX
		Y(j)=NY
		Z(J)=NZ
		X(k)=MX
		Y(k)=MY
		Z(k)=MZ
		ICA(i)=NV1
		ICA(j)=NV2
		ICA(k)=NV3
		I3flip=I3flip+1
	erestsum=ectre
	arestsum=actre	
	
		go to 8885

	ELSE
c
c	Move rejected - restore the initial state of the chain			
	ICA(i)=ii
	ICA(j)=jj
	ICA(k)=kk
	x(j)=jx
	y(j)=jy
	z(j)=jz
	x(k)=kx
	y(k)=ky
	z(k)=kz

	ENDIF	


 8885	CONTINUE

	go to 7775


c	REMOVED LARGE SCALE MOVES











c	"REPTATION" -type moves, two bond progress (4 to 22 bonds)


 	j=INT(arand(ise)*AL19)+3
 	ppp=INT(arand(ise)*18.9999)+2
	i=j-1
	k=j+ppp
	L=k+1
	if(L.gt.lenf1) go to 4001
		
 		do ii=i,k
		icao(ii)=ica(ii)
		enddo
  
		IF(arand(ise).gt.0.5) THEN
c		a move forward
		if(.not.GOODC(ica(l-3),ica(l))) go to 4001
		if(.not.GOODC(ica(i-1),ica(k-1))) go to 4001
		if(.not.GOODC(ica(k),ica(i))) go to 4001
		ii=ica(k-1)
		kk=ica(k)
			do pp=k,j+1,-1
			ica(pp)=ica(pp-2)
			enddo
			ica(i)=ii
			ica(j)=kk
		ELSE

c		a move backward
		if(.not.GOODC(ica(i-1),ica(j+1))) go to 4001
		if(.not.GOODC(ica(k),ica(i))) go to 4001
		if(.not.GOODC(ica(j),ica(l))) go to 4001
		ii=ica(i)
		kk=ica(j)
			do pp=i,k-2
			ica(pp)=ica(pp+2)
			enddo
			ica(k-1)=ii
			ica(k)=kk	

		ENDIF

		do kkkk=i,k
		icab(kkkk)=ica(kkkk)
		enddo
		do kkkk=j,k
		pp=kkkk-1
		ii=ica(pp)
		x(kkkk)=x(pp)+vx(ii)
		y(kkkk)=y(pp)+vy(ii)
		z(kkkk)=z(pp)+vz(ii)
		enddo
		
	IF(LOOK(i,l)) THEN
c	No overalps of Cas and Cbs
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE
c	

		ENEW=EHB(i,l,1)+ESHORT(i,l,1)


		do kkkk=i,k
		ica(kkkk)=icao(kkkk)
		enddo
		do kkkk=j,k
		pp=kkkk-1
		ii=ica(pp)
		x(kkkk)=x(pp)+vx(ii)
		y(kkkk)=y(pp)+vy(ii)
		z(kkkk)=z(pp)+vz(ii)
		enddo

		EOLD=EHB(i,l,-1)+ESHORT(i,l,-1)
		DE=ENEW-EOLD+dord	

		IF(DE.GT.0.0) THEN
		if(arand(ise).gt.EXP(-DE/ATEMP)) then
 	   	ICNT=ICNTO
        	SUMCT=SUMCTO 		
		go to 4001
		endif
		ENDIF
		
		
c
c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION
c
        ICNTO=ICNT
        SUMCTO=SUMCT 	

		ENERG=ENERG+DE

		do kkkk=i,k
		ica(kkkk)=icab(kkkk)
		enddo
		do kkkk=j,k
		pp=kkkk-1
		ii=ica(pp)
		x(kkkk)=x(pp)+vx(ii)
		y(kkkk)=y(pp)+vy(ii)
		z(kkkk)=z(pp)+vz(ii)
		enddo			

		I8flip=I8flip+1
	erestsum=ectre
	arestsum=actre	
		go to 4001

	ELSE
c
c	Move rejected - restore the initial state of the chain
c


		do kkkk=i,k
		ica(kkkk)=icao(kkkk)
		enddo
		do kkkk=j,k
		pp=kkkk-1
		ii=ica(pp)
		x(kkkk)=x(pp)+vx(ii)
		y(kkkk)=y(pp)+vy(ii)
		z(kkkk)=z(pp)+vz(ii)
		enddo			

	ENDIF	

c
c
c	FOUR-to-22 BOOND MOVES 
c

 4001	j=INT(arand(ise)*AL19)+3
 	ppp=INT(arand(ise)*18.9999)+2
	i=j-1
	k=j+ppp
	L=k+1
	if(L.gt.lenf1) go to 7775
		
 8809		IF(arand(ise).gt.0.25) THEN
 		mx=INT(arand(ise)*2.999999)-1
		my=INT(arand(ise)*2.999999)-1
		mz=INT(arand(ise)*2.999999)-1
		if((mx*mx+my*my+mz*mz).eq.0) go to 8809	
		ELSE	
 		mx=INT(arand(ise)*4.999999)-2
		my=INT(arand(ise)*4.999999)-2
		mz=INT(arand(ise)*4.999999)-2
		if((mx*mx+my*my+mz*mz).eq.0) go to 8809	
		ENDIF
	
	ii=ica(i)
	wx=vx(ii)+mx
	wy=vy(ii)+my
	wz=vz(ii)+mz
	ir=wx*wx+wy*wy+wz*wz
	if(ir.lt.29) go to 8809
	if(ir.gt.49) go to 8809
	iii=vector(wx,wy,wz)
	if(.NOT.GOODC(ica(i-1),iii)) go to 7775
	if(.NOT.GOODC(iii,ica(j))) go to 7775
	
	kk=ica(k)
	wx=vx(kk)-mx
	wy=vy(kk)-my
	wz=vz(kk)-mz
	ir=wx*wx+wy*wy+wz*wz
	if(ir.lt.29) go to 7775
	if(ir.gt.49) go to 7775	
	kkk=vector(wx,wy,wz)
	if(.NOT.GOODC(ica(k-1),kkk)) go to 7775
	if(.NOT.GOODC(kkk,ica(l))) go to 7775	


		ICA(i)=iii
		ICA(k)=kkk
		do kkkk=j,k
		x(kkkk)=x(kkkk)+mx
		y(kkkk)=y(kkkk)+my
		z(kkkk)=z(kkkk)+mz
		enddo
		
	IF(LOOK(i,l)) THEN
c	No overalps of Cas and Cbs
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE
c	

		ENEW=EHB(i,l,1)+ESHORT(i,l,1)

		ICA(i)=ii
		ICA(k)=kk
		do kkkk=j,k
		x(kkkk)=x(kkkk)-mx
		y(kkkk)=y(kkkk)-my
		z(kkkk)=z(kkkk)-mz
		enddo

		EOLD=EHB(i,l,-1)+ESHORT(i,l,-1)
		DE=ENEW-EOLD+dord	

		IF(DE.GT.0.0) THEN
		if(arand(ise).gt.EXP(-DE/ATEMP)) then
 	   	ICNT=ICNTO
        	SUMCT=SUMCTO 		
		go to 7775
		endif
		ENDIF
		
		
c
c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION
c
        ICNTO=ICNT
        SUMCTO=SUMCT 	

		ENERG=ENERG+DE
			
		do kkkk=j,k
		x(kkkk)=x(kkkk)+mx
		y(kkkk)=y(kkkk)+my
		z(kkkk)=z(kkkk)+mz
		enddo
		ICA(i)=iii
		ICA(k)=kkk
		I6flip=I6flip+1
	erestsum=ectre
	arestsum=actre	
		go to 7775

	ELSE
c
c	Move rejected - restore the initial state of the chain
c			
		ICA(i)=ii
		ICA(k)=kk
		do kkkk=j,k
		x(kkkk)=x(kkkk)-mx
		y(kkkk)=y(kkkk)-my
		z(kkkk)=z(kkkk)-mz
		enddo

	ENDIF	

c

c
c	TWO BOND CHAIN END MOVES- random selection of two bonds
c
c	N TERMINUS
c
 7775	CONTINUE


c
c	Accumulate the contact map averages
c   

 	JV3=ICA(3)
 60	NV2=INT(arand(ise)*799.99 )+1
	NV1=INT(arand(ise)*799.99 )+1

		if(.NOT.GOODC(nv1,nv2)) go to 60
		if(.NOT.GOODC(nv2,jv3)) go to 60
		
		
		kx=x(1)
		ky=y(1)
		kz=z(1)
		px=x(2)
		py=y(2)
		pz=z(2)
		ica1=ica(1)
		ica2=ica(2)

	X(2)=X(3)-vx(nv2)
	Y(2)=Y(3)-vy(nv2)
	Z(2)=Z(3)-vz(nv2)

		X(1)=X(2)-vx(nv1)
		Y(1)=Y(2)-vy(nv1)
		Z(1)=Z(2)-vz(nv1)
	ICA(1)=nv1
	ICA(2)=nv2
		
	if(LOOK(2,3)) THEN
	
c
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE
c	
		
		ENEW=EHB(2,3,1)+ESHORT(2,3,1)

		x(1)=kx
		y(1)=ky
		z(1)=kz
		x(2)=px
		y(2)=py
		z(2)=pz
		ICA(1)=ica1
		ICA(2)=ica2	

		EOLD=EHB(2,3,-1)+ESHORT(2,3,-1)

		DE=ENEW-EOLD+dord	

		IF(DE.GT.0.0) THEN
		if(arand(ise).gt.EXP(-DE/ATEMP)) then	
	    ICNT=ICNTO
        	SUMCT=SUMCTO 	
		go to 81
		endif
		ENDIF
		
c
c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION
c
        ICNTO=ICNT
        SUMCTO=SUMCT 	

		ENERG=ENERG+DE	
		X(2)=X(3)-vx(nv2)
		Y(2)=Y(3)-vy(nv2)
		Z(2)=Z(3)-vz(nv2)

		X(1)=X(2)-vx(nv1)
		Y(1)=Y(2)-vy(nv1)
		Z(1)=Z(2)-vz(nv1)
		ICA(1)=nv1
		ICA(2)=nv2
	erestsum=ectre
	arestsum=actre	
	
		go to 81
	ELSE
c
c	Move rejected - restore the initial state of the chain	 	
	x(1)=kx
	y(1)=ky
	z(1)=kz
	x(2)=px
	y(2)=py
	z(2)=pz
	ICA(1)=ica1
	ICA(2)=ica2
	ENDIF

c
c	C TERMINUS
c
 81	JV3=ICA(LENF3)
 80	NV2=INT(arand(ise)*799.99 )+1
	NV1=INT(arand(ise)*799.99 )+1
		if(.NOT.GOODC(nv2,nv1)) go to 80
		if(.NOT.GOODC(jv3,nv2)) go to 80

		
		kx=x(lenf)
		ky=y(lenf)
		kz=z(lenf)
		px=x(lenf1)
		py=y(lenf1)
		pz=z(lenf1)
		ica1=ica(lenf1)
		ica2=ica(lenf2)

	X(lenf1)=X(lenf2)+vx(nv2)
	Y(lenf1)=Y(lenf2)+vy(nv2)
	Z(lenf1)=Z(lenf2)+vz(nv2)

		X(lenf)=X(lenf1)+vx(nv1)
		Y(lenf)=Y(lenf1)+vy(nv1)
		Z(lenf)=Z(lenf1)+vz(nv1)
		
	ICA(lenf1)=nv1
	ICA(lenf2)=nv2
		
	if(LOOK(lenf2,lenf1)) THEN
c
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE
c

		ENEW=EHB(lenf2,lenf1,1)+ESHORT(lenf2,lenf1,1)

		x(lenf)=kx
		y(lenf)=ky
		z(lenf)=kz
		x(lenf1)=px
		y(lenf1)=py
		z(lenf1)=pz
		ICA(lenf1)=ica1
		ICA(lenf2)=ica2
		
		EOLD=EHB(lenf2,lenf1,-1)+ESHORT(lenf2,lenf1,-1)

		DE=ENEW-EOLD+dord	

		IF(DE.GT.0.0) THEN
		if(arand(ise).gt.EXP(-DE/ATEMP)) then
 	   	ICNT=ICNTO
        	SUMCT=SUMCTO 	
		go to 7777
		endif
		ENDIF
				
c
c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION
c
        ICNTO=ICNT
        SUMCTO=SUMCT 	
			
		ENERG=ENERG+DE
		X(lenf1)=X(lenf2)+vx(nv2)
		Y(lenf1)=Y(lenf2)+vy(nv2)
		Z(lenf1)=Z(lenf2)+vz(nv2)

		X(lenf)=X(lenf1)+vx(nv1)
		Y(lenf)=Y(lenf1)+vy(nv1)
		Z(lenf)=Z(lenf1)+vz(nv1)
		
		ICA(lenf1)=nv1
		ICA(lenf2)=nv2
	erestsum=ectre
	arestsum=actre

		go to 7777
		

	ELSE
c
c	Move rejected - restore the initial state of the chain		
	x(lenf)=kx
	y(lenf)=ky
	z(lenf)=kz
	x(lenf1)=px
	y(lenf1)=py
	z(lenf1)=pz
	ICA(lenf1)=ica1
	ICA(lenf2)=ica2

	ENDIF

 7777	CONTINUE
 
 
 	EREPLICA(itemp)=energ
 		

c
c	end of ITEMP's replica	
c
	do i=1,lenf
	xreps(i,itemp)=x(i)
	yreps(i,itemp)=y(i)
	zreps(i,itemp)=z(i)
	enddo
 
	ENDDO
c

c	Search for the lowest energy replica
	energ=EREPLICA(1)
	ITEMP=1
	do i=1, REPLICAS
	if(ereplica(i).lt.energ) then
	energ=ereplica(i)
	ITEMP=i
	endif
	enddo

	do i=1,lenf
	x(i)=xreps(i,itemp)
	y(i)=yreps(i,itemp)
	z(i)=zreps(i,itemp)
	enddo

		DO I=1,LENF1
		J=I+1
		WX=X(J)-X(I)
		WY=Y(J)-Y(I)
		WZ=Z(J)-Z(I)
		ICA(I)=VECTOR(WX,WY,WZ)
c
c	Define chain of vectors ICA() =1,2,...800
c
		ENDDO


	iiii=(IDDDUM-1)*10+ICYCLE
	
	write(23,*) energ
	WRITE(9,716) IIII, LENF, energ, aarand,abrand
	WRITE(40,716) IIII, LENF2,ENERG, aarand, abrand

 716	format(2I6, F12.2, 2F7.4)
	write(9,713)    (x(i),y(i),z(i),i=1,lenf)
 713	format(12I5)
			
	c=float(ICYCLE)*PHOT*LENF2*REPLICAS
	a=I2flip/c/10.0
	b=I3flip/c
	d=I6flip/c
	dd=I8flip/c
	
	if(icycle.gt.1) then
	axr=float(inreps)/((icycle-1)*float((replicas-1)*replicas)/2.0)
	else
	axr=0.0
	endif
	
	kkk=0
	do i=2,lenf1
	
        k=seq(i)
        ii=ICA(i-1)
        jj=ICA(i)
        ax=(x(i)+GX(ii,jj,k))
        ay=(y(i)+GY(ii,jj,k))
        az=(z(i)+GZ(ii,jj,k))
        
	WRITE(40,7161) (i-1), ax, ay, az, aa(k)


	fff=ax*ax/amx+ay*ay/amy+az*az/amz
	if(fff.lt.3.0) kkk=KKK+1
	enddo

 7161	FORMAT(i5, 3F7.1, 5X,A3)
	
	amxx=sqrt(amx)
	amyy=sqrt(amy)
	amzz=sqrt(amz)
      	Write(6,714)ICYCLE,a,b,d,dd,ENERG,acrit,amxx,amyy,amzz,axr,
     * 	kkk, SUMCT, itemp
 714	FORMAT(i5, 4F7.4,f10.2,4f5.1,f6.3, i4,i4, i6)	
	
	
c	 	
 	
 7000 	CONTINUE
 
 	ENDDO


	
c	===============================================================	
 
c
c	CHECK THE GEOMETRY and WRITE THE FINAL CHAIN
c 
	
	
 	DO J=1,LENF
	I=J-1
	if(j.ne.1.and.j.ne.lenf) THEN
	II=ICA(I)
	JJ=ICA(J)
	IF(.NOT.GOODC(II,JJ)) THEN
	WRITE(6,8112)I,J,vx(ii),vy(ii),vz(ii),vx(jj),vy(jj),vz(jj)
 8112		FORMAT(5X,'WARNING2 -WRONG INPUT CHAIN - VECTORS ',8I4)
 		STOP
		ENDIF
		ENDIF
	ENDDO 
	
c
c	PUT ALSO TEST FOR OVERLAPS OF Ca with Cb - initial conformations
c	from program NatCa.f may sometimes have overlaps - they relax
c	rapidly
	kcc=0
	kcb=0
	kbb=0
	do i=2,lenf3
	ix=x(i)
	iy=y(i)
	iz=z(i)
	bx=ix+CBX(ica(i-1),ica(i))
	by=iy+CBY(ica(i-1),ica(i))
	bz=iz+CBZ(ica(i-1),ica(i))
	do j=i+2,lenf1
	jx=x(j)
	jy=y(j)
	jz=z(j)
	cx=jx+CBX(ica(j-1),ica(j))
	cy=jy+CBY(ica(j-1),ica(j))
	cz=jz+CBZ(ica(j-1),ica(j))
	ir=(ix-jx)**2+(iy-jy)**2+(iz-jz)**2
	if(ir.lt.43) kcc=kcc+1
	if(SEQ(i).NE.0  .AND. SEQ(j).NE.0) then
	ar=(bx-cx)**2+(by-cy)**2+(bz-cz)**2
	if(ar.lt.36.0) kbb=kbb+1
	endif
	if(SEQ(j).NE.0) THEN
	ar=(ix-cx)**2+(iy-cy)**2+(iz-cz)**2
	if(ar.lt.36.0) kcb=kcb+1
	ENDIF
	if(SEQ(i).NE.0) THEN		
	ar=(bx-jx)**2+(by-jy)**2+(bz-jz)**2			
	if(ar.lt.36.0) kcb=kcb+1
	endif
	enddo
	enddo
	
	
	write(6,819) energ
c
c	Compute from scratch the final energy for the consistency
c	comparison with the cumulated value of energy - TEST

        ICNTO=0
        SUMCTO=0 		
        ICNT=0
        SUMCT=0 
	
	erestsum=0.0
	arestsum=0.0		
	erold=0.0
	arold=0.0		
		

	energ=EHB(2,lenf1,1)+ESHORT(2,lenf1,10)

	erestsum=ernew
	arestsum=arnew	
		
	if(sumct.gt.0) energ=energ+ardscal*abs(ICNT/float(SUMCT)-acorder) 
     *	+cstr*abs(float(SUMCT)-contt)
     	
        ICNTO=ICNT
        SUMCTO=SUMCT 	  

	write(6,*)
	write(6,819) energ
 819	FORMAT('FINAL ENERGY', f10.2)
	
	write(6,*)
	write(6,*) 'Ca-Ca overlaps  ', kcc	
	write(6,*) 'Ca-Cb overlaps  ', kcb		
	write(6,*) 'Cb-Cb overlaps  ', kbb	


c
	OPEN(UNIT=19,FILE='ACHAINS_NEW',STATUS='UNKNOWN')
C
	REWIND(UNIT=19)
	do k=1,REPLICAS
	write(19,*) LENF
	DO I=1,LENF
	WRITE(19,8000) XREPS(I,k),YREPS(I,k),ZREPS(I,k)
	ENDDO
	enddo
	
 8000	format(1x,3i5)

	
	WRITE(6,*)
	WRITE(6,*) '   	FINAL REPLICA DISTRIBUTION'
	do i=1, REPLICAS
	write(6,5005) i,ATREP(i),EREPLICA(i), SWAPS(i)
	enddo
 5005	format('   replica #',I4,'   T =', f6.2,'  E =',f8.1, I10)	
 	
	write(6,*)
 	write(6,*)
 	WRITE(6,*) ' Number of restraints = ',NREST
 	IF(NREST.GT.0) THEN
 	nsat=0
 	do k=1, NREST
 	i=irest(k)
 	j=jrest(k)
 	
        kkk=seq(i)
        ii=ICA(i-1)
        jj=ICA(i)
               
        ax=(xreps(i,1)+GX(ii,jj,kkk))
        ay=(yreps(i,1)+GY(ii,jj,kkk))
        az=(zreps(i,1)+GZ(ii,jj,kkk))
         	
       	kk=ica(j-1)
        ll=ica(j)
        jjj=seq(j)
        bx=xreps(j,1)+GX(kk,ll,jjj)
        by=yreps(j,1)+GY(kk,ll,jjj)
        bz=zreps(j,1)+GZ(kk,ll,jjj) 
 
  	ar=(ax-bx)**2+(ay-by)**2+(az-bz)**2
  	ic=0
  	if(ar.lt.120.0) ic=1
  	nsat=nsat+ic
  	write (6,8765) k, i,j, ar, ic
 8765	format(3i4,F8.0,i4)
 	enddo
 	write(6,*) 'satisfied ',nsat,'   per ',NREST,'  of contacts'
	write(6,*)
	write(6,*)
 	ENDIF	 	

 	STOP
 	END	
	
c
c	END OF THE MAIN PROGRAM
c
c	***************************************************************		
	
	
	
	FUNCTION LOOK(jj,kk)
	IMPLICIT INTEGER(I-Z)
	LOGICAL LOOK

		PARAMETER(NDIM=400)

	COMMON /SEQE/  	 SEQ(ndim), SEC(ndim) 	 	
	COMMON /CHAINS/  ICA(NDIM), X(NDIM), Y(NDIM), Z(NDIM)
	COMMON /LENGHTS/ LENF1,LENF,PROD(800,800) 
	COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
c
c	CHECKS Ca-Ca, Ca-Cb and Cb-Cb distances and overlaps
c

	DO k=jj,kk

	px=x(k)
	py=y(k)
	pz=z(k)
	if(SEQ(k).NE.0) THEN
	nv1=ica(k-1)
	nv2=ica(k)
	bx=CBX(nv1,nv2)+px
	by=CBY(nv1,nv2)+py
	bz=CBZ(nv1,nv2)+pz
	ENDIF
	
	istart=k-2
	if(istart.gt.1) THEN
	do i=2,istart
	irx=(px-x(i))**2
	if(irx.lt.120) then
		iry=(py-y(i))**2
		if(iry.lt.120) then
		ir=irx+iry+(pz-z(i))**2
		if(ir.lt.120) then
		
c	detect detailed overlaps Ca-Ca, Cb-Cb, Ca-Cb, Cb-Ca

		if(ir.lt.43) then
		LOOK=.FALSE.
		RETURN
		endif
		
		if(SEQ(k).NE.0) THEN
		ar=(bx-x(i))**2+(by-y(i))**2+(bz-z(i))**2
		if(ar.lt.36.0) then
		LOOK=.FALSE.
		RETURN
		endif
		ENDIF
		
		IF(SEQ(i).NE.0) THEN
		icai=ICA(i-1)
		icaj=ICA(i)
		ax= CBX(icai,icaj)+x(i)
		ay= CBY(icai,icaj)+y(i)
		az= CBZ(icai,icaj)+z(i)
		
		IF(SEQ(k).NE.0) THEN
		ar=(ax-bx)**2+(ay-by)**2+(az-bz)**2
		if(ar.lt.36.0) then
		LOOK=.FALSE.
		RETURN
		endif
		ENDIF
		
		ar=(ax-px)**2+(ay-py)**2+(az-pz)**2		
		if(ar.lt.36.0) then
		LOOK=.FALSE.
		RETURN
		endif
		ENDIF		
		
		endif
		endif
	endif
	enddo
	ENDIF
	
	
	iend=max(k+2,kk+1)
	if(iend.lt.lenf) THEN
	do i=iend,lenf1
	irx=(px-x(i))**2
	if(irx.lt.120) then
		iry=(py-y(i))**2
		if(iry.lt.120) then
		ir=irx+iry+(pz-z(i))**2
		if(ir.lt.120) then
		
c	detect detailed overlaps Ca-Ca, Cb-Cb, Ca-Cb, Cb-Ca

		if(ir.lt.43) then
		LOOK=.FALSE.
		RETURN
		endif

		
		if(SEQ(k).NE.0) THEN
		ar=(bx-x(i))**2+(by-y(i))**2+(bz-z(i))**2
		if(ar.lt.36.0) then
		LOOK=.FALSE.
		RETURN
		endif
		ENDIF
		
		if(SEQ(i).NE.0) THEN
		icai=ICA(i-1)
		icaj=ICA(i)
		ax= CBX(icai,icaj)+x(i)
		ay= CBY(icai,icaj)+y(i)
		az= CBZ(icai,icaj)+z(i)
		
		IF(SEQ(k).NE.0) THEN		
		ar=(ax-bx)**2+(ay-by)**2+(az-bz)**2
		if(ar.lt.36.0) then
		LOOK=.FALSE.
		RETURN
		endif
		ENDIF
				
		ar=(ax-px)**2+(ay-py)**2+(az-pz)**2		
		if(ar.lt.36.0) then
		LOOK=.FALSE.
		RETURN
		endif	
		
		endif
		endif
		ENDIF	

	endif
	enddo
	ENDIF


	ENDDO	
		
c
c	NO HARD CORE OVERLAPPS DETECTED
c	
	LOOK=.TRUE.
	
	RETURN
	END
	


c	***************************************************************	
c	A random number generator by AK - tested
	FUNCTION arand(ise)
	common/arandom/a,b,c,d
	g=(a+b+c)*d
	g=g-aint(g)
	a=b
	b=c
	c=g
	arand=g
	return
	end
	
C	***************************************************************	



	FUNCTION EHB(jjjj,kkkk,ISTAT)
	
	IMPLICIT INTEGER(I-Z)

		PARAMETER(NDIM=400)	

	COMMON /SEQE/  	 SEQ(ndim), SEC(ndim) 		
	COMMON /CHAINS/  ICA(NDIM), X(NDIM), Y(NDIM), Z(NDIM)
	COMMON /LENGHTS/ LENF1,LENF,PROD(800,800) 
	COMMON /HB/      HBX(800,800),HBY(800,800),HBZ(800,800)
	COMMON /BISEC/   CAX(800,800),CAY(800,800),CAZ(800,800)
	COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
	COMMON /ENERGY/	 EHBOND, ESC, EREP, BAT, EHBIJ(ndim,ndim)

	COMMON /pair/	apa(ndim,ndim,2,2),app(ndim,ndim,2,2),
     *			apm(ndim,ndim,2,2)
	COMMON /size/	arla(0:19,0:19,2,2),arlm(0:19,0:19,2,2),
     *			arlp(0:19,0:19,2,2)
	COMMON /sizea/	ala(0:19,0:19,2,2),alm(0:19,0:19,2,2),
     *			alp(0:19,0:19,2,2)
	COMMON /KBB/ kbin(800,800)

	COMMON /one/ 	acrit, contt, eoinp(0:19,0:100)
      COMMON /order/acorder,ardscal,SUMCT,SUMCTO,ICNT,ICNTO,dord,cstr 
	COMMON /SG/ GX(800,800,0:19),GY(800,800,0:19),GZ(800,800,0:19)
	common /icgg/ icg(ndim)	
	LOGICAL contact
c
c	OTHER INTERACTIONS CAN BE ADDED HERE
c
c

	EHB=0.0
	EREP4=EREP/4.0

	if(ISTAT.gt.0) THEN
	ICNT=ICNTO
	SUMCT=SUMCTO
	ENDIF

						 
	DO k=jjjj,kkkk
	im=seq(k)
	px=x(k)
	py=y(k)
	pz=z(k)

		nv1=ica(k-1)
		nv2=ica(k)

		bx=HBX(nv1,nv2)
		by=HBY(nv1,nv2)
		bz=HBZ(nv1,nv2)
		hx=CAX(nv1,nv2)
		hy=CAY(nv1,nv2)
		hz=CAZ(nv1,nv2)	
		agx=px+GX(nv1,nv2,im)	
		agy=py+GY(nv1,nv2,im)	
		agz=pz+GZ(nv1,nv2,im)

		kb= kbin(nv1,nv2)

	istart=k-2
	if(istart.gt.1) THEN
	do 1001 i=2,istart

	irx=(px-x(i))**2
	if(irx.lt.400) then
	iry=(py-y(i))**2
	if(iry.lt.400) then
	ir=irx+iry+(pz-z(i))**2
	if(ir.lt.400) then
	
		xi=x(i)
		yi=y(i)
		zi=z(i)

	ica1=ica(i-1)
	ica2=ica(i)
	ib=kbin(ica1,ica2)
	in=seq(i)

		bgx=xi+GX(ica1,ica2,in)	
		bgy=yi+GY(ica1,ica2,in)	
		bgz=zi+GZ(ica1,ica2,in)

c	Repulsion between the peptide bond and Ca


	IF(iabs(i-k).gt.3) THEN
	k1=k+1
	i1=i+1
	dx=(px+x(k1))/2.0-xi
	dy=(py+y(k1))/2.0-yi
	dz=(pz+z(k1))/2.0-zi
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.54.0) EHB=EHB+EREP

	dx=px-(xi+x(i1))/2.0
	dy=py-(yi+y(i1))/2.0
	dz=pz-(zi+z(i1))/2.0
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.54.0) EHB=EHB+EREP
	
	k1=k-1
	i1=i-1
	dx=(px+x(k1))/2.0-xi
	dy=(py+y(k1))/2.0-yi
	dz=(pz+z(k1))/2.0-zi
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.54.0) EHB=EHB+EREP

	dx=px-(xi+x(i1))/2.0
	dy=py-(yi+y(i1))/2.0
	dz=pz-(zi+z(i1))/2.0
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.54.0) EHB=EHB+EREP
	
		

c	side group-Ca distance
		if(im.gt.1) then
		aks=(agx-xi)**2+(agy-yi)**2+(agz-zi)**2
		if(aks.lt.43.0) EHB=EHB+EREP
		endif

c	side group-Ca distance
		if(in.gt.1) then
		aks=(bgx-px)**2+(bgy-py)**2+(bgz-pz)**2
		if(aks.lt.43.0) EHB=EHB+EREP
		endif
	ENDIF

c	side group- side group distance		
		ar=(agx-bgx)**2+(agy-bgy)**2+(agz-bgz)**2
						
		fx=CAX(ica1,ica2)
		fy=CAY(ica1,ica2)
		fz=CAZ(ica1,ica2)

c	contact orientation (parallel, othogonal, or antiparallel)		
		aor=hx*fx+hy*fy+hz*fz
		
			IF(aor.gt.0.5) THEN
			if(ar.gt.arlp(in,im,ib,kb)) THEN
			if(ar.lt.alp(in,im,ib,kb)) then 
					EHB=EHB+app(i,k,ib,kb) 
c	contact order update
        				ICNT=ICNT+istat*iabs(i-k)
c	total number of contact update
        				SUMCT=SUMCT+ISTAT		     
					endif
			else
			EHB=EHB+EREP
c	repulsons from side group overlapps weaker than for Ca or Cb
			endif			
			ELSE
			IF(aor.lt.-0.5) THEN
			if(ar.gt.arla(in,im,ib,kb)) THEN
			if(ar.lt.ala(in,im,ib,kb)) then
					EHB=EHB+apa(i,k,ib,kb) 
        				ICNT=ICNT+istat*iabs(i-k)
        				SUMCT=SUMCT+ISTAT
					endif
			else
			EHB=EHB+EREP
			endif			
			ELSE
			if(ar.gt.arlm(in,im,ib,kb)) THEN
			if(ar.lt.alm(in,im,ib,kb)) then
					EHB=EHB+apm(i,k,ib,kb) 
        				ICNT=ICNT+istat*iabs(i-k)
        				SUMCT=SUMCT+ISTAT
					endif
			else
			EHB=EHB+EREP
			endif						
			ENDIF
			ENDIF



	IF(aor.gt.0.0) THEN	
		idist=iabs(i-k)
		if(idist.gt.2) THEN					
		if(ir.lt.100) THEN

		IF(idist.ne.4) THEN
			IF(idist.gt.3) then	
			if(sec(k).eq.2.OR.sec(i).eq.2) go to 1001
			endif

		cx=HBX(ica1,ica2)
		cy=HBY(ica1,ica2)
		cz=HBZ(ica1,ica2)
			
c	COOPERATIVITY OF THE H-bonds, plates of contacting fragments
c		must be almost parallel

	
	if(prod(nv1,ica1).GT.0.AND.prod(nv2,ica2).GT.0) THEN

c	chains must be parallel or (****)
c
			IF(idist.eq.3) then	
			if(sec(k-1).eq.4.OR.sec(i+1).eq.4) go to 1001
			endif
	if(idist.gt.3.AND.idist.lt.20) go to 1001
c	distance dependent attraction of Ca in the H-bond type
c	geometry
	id=1
	ELSE
	if(prod(nv1,ica2).gt.0.OR.prod(nv2,ica1).gt.0) go to 1001
	id=-1
	ENDIF
	
		ax=(px-xi)
		ay=(py-yi)
		az=(pz-zi)


	delta=0.0
	k1=k+1
	i1=i+id
	dx=(px+x(k1)-xi-x(i1))/2.0
	dy=(py+y(k1)-yi-y(i1))/2.0
	dz=(pz+z(k1)-zi-z(i1))/2.0
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.80.0) delta=delta+0.5

	k1=k-1
	i1=i-id
	dx=(px+x(k1)-xi-x(i1))/2.0
	dy=(py+y(k1)-yi-y(i1))/2.0
	dz=(pz+z(k1)-zi-z(i1))/2.0
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.80.0) delta=delta+0.5



		ar=(bx-ax)**2+(by-ay)**2+(bz-az)**2
		if(ar.lt.10.0) then
			EHB=EHB+EHBIJ(i,k)*(0.5+delta)
			else
		ar=(bx+ax)**2+(by+ay)**2+(bz+az)**2
		if(ar.lt.10.0) then
			EHB=EHB+EHBIJ(i,k)*(0.5+delta)
			endif
			endif	

		ar=(cx-ax)**2+(cy-ay)**2+(cz-az)**2
		if(ar.lt.10.0) then
			EHB=EHB+EHBIJ(i,k)*(0.5+delta)
			else
		ar=(cx+ax)**2+(cy+ay)**2+(cz+az)**2
		if(ar.lt.10.0) then
			EHB=EHB+EHBIJ(i,k)*(0.5+delta) 
			endif
			endif

		ENDIF
		ENDIF
		ENDIF
		ENDIF
	endif
	endif
	endif
 1001	CONTINUE
	ENDIF
	
	
	iend=max(k+2,kkkk+1)
	if(iend.lt.lenf) THEN
	do 1002 i=iend,lenf1

	irx=(px-x(i))**2
	if(irx.lt.400) then
	iry=(py-y(i))**2
	if(iry.lt.400) then
	ir=irx+iry+(pz-z(i))**2
	if(ir.lt.400) then
	
		xi=x(i)
		yi=y(i)
		zi=z(i)
	ica1=ica(i-1)
	ica2=ica(i)
	ib=kbin(ica1,ica2)
	in=seq(i)
		bgx=xi+GX(ica1,ica2,in)	
		bgy=yi+GY(ica1,ica2,in)	
		bgz=zi+GZ(ica1,ica2,in)

c	Repulsion between the peptide bond and Ca

	IF(iabs(i-k).gt.3) THEN
	k1=k+1
	i1=i+1
	dx=(px+x(k1))/2.0-xi
	dy=(py+y(k1))/2.0-yi
	dz=(pz+z(k1))/2.0-zi
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.54.0) EHB=EHB+EREP
	
	dx=px-(xi+x(i1))/2.0
	dy=py-(yi+y(i1))/2.0
	dz=pz-(zi+z(i1))/2.0
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.54.0) EHB=EHB+EREP

	k1=k-1
	i1=i-1
	dx=(px+x(k1))/2.0-xi
	dy=(py+y(k1))/2.0-yi
	dz=(pz+z(k1))/2.0-zi
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.54.0) EHB=EHB+EREP
	
	dx=px-(xi+x(i1))/2.0
	dy=py-(yi+y(i1))/2.0
	dz=pz-(zi+z(i1))/2.0
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.54.0) EHB=EHB+EREP


		if(im.gt.1) then
		aks=(agx-xi)**2+(agy-yi)**2+(agz-zi)**2
		if(aks.lt.43.0) EHB=EHB+EREP
		endif	
		
		if(in.gt.1) then
		aks=(bgx-px)**2+(bgy-py)**2+(bgz-pz)**2
		if(aks.lt.43.0) EHB=EHB+EREP
		endif
	ENDIF	
	
		ar=(agx-bgx)**2+(agy-bgy)**2+(agz-bgz)**2

		fx=CAX(ica1,ica2)
		fy=CAY(ica1,ica2)
		fz=CAZ(ica1,ica2)	

		aor=hx*fx+hy*fy+hz*fz
	

	
			IF(aor.gt.0.5) THEN
			if(ar.gt.arlp(im,in,kb,ib)) THEN
			if(ar.lt.alp(im,in,kb,ib)) then
					EHB=EHB+app(k,i,kb,ib) 
        				ICNT=ICNT+istat*iabs(i-k)
        				SUMCT=SUMCT+ISTAT
					endif
				
			else
			EHB=EHB+EREP
			endif			
			ELSE
			IF(aor.lt.-0.5) THEN
			if(ar.gt.arla(im,in,kb,ib)) THEN
			if(ar.lt.ala(im,in,kb,ib)) then
					EHB=EHB+apa(k,i,kb,ib) 
        				ICNT=ICNT+istat*iabs(i-k)
        				SUMCT=SUMCT+ISTAT
					endif
			else
			EHB=EHB+EREP
			endif			
			ELSE
			if(ar.gt.arlm(im,in,kb,ib)) THEN
			if(ar.lt.alm(im,in,kb,ib)) then
					EHB=EHB+apm(k,i,kb,ib) 
        				ICNT=ICNT+istat*iabs(i-k)
        				SUMCT=SUMCT+ISTAT
					endif
			else
			EHB=EHB+EREP
			endif						
			ENDIF
			ENDIF
			
	IF(aor.gt.0.0) THEN	
		idist=iabs(i-k)
		if(idist.gt.2) THEN					
		if(ir.lt.100) THEN

		IF(idist.ne.4) THEN
			IF(idist.gt.3) then	
			if(sec(k).eq.2.OR.sec(i).eq.2) go to 1002
			endif			

		cx=HBX(ica1,ica2)
		cy=HBY(ica1,ica2)
		cz=HBZ(ica1,ica2)
			
c		COOPERATIVITY


	if(prod(nv1,ica1).GT.0.AND.prod(nv2,ica2).GT.0) THEN
			IF(idist.eq.3) then	
			if(sec(k+1).eq.4.OR.sec(i-1).eq.4) go to 1002
			endif
	if(idist.gt.3.AND.idist.lt.20) go to 1002
	id=1
	ELSE
	if(prod(nv1,ica2).gt.0.OR.prod(nv2,ica1).gt.0) go to 1002
	id=-1
	ENDIF		
	
		ax=(px-xi)
		ay=(py-yi)
		az=(pz-zi)


	delta=0.0
	k1=k+1
	i1=i+id
	dx=(px+x(k1)-xi-x(i1))/2.0
	dy=(py+y(k1)-yi-y(i1))/2.0
	dz=(pz+z(k1)-zi-z(i1))/2.0
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.80.0) delta=delta+0.5

	k1=k-1
	i1=i-id
	dx=(px+x(k1)-xi-x(i1))/2.0
	dy=(py+y(k1)-yi-y(i1))/2.0
	dz=(pz+z(k1)-zi-z(i1))/2.0
	arrr=dx*dx+dy*dy+dz*dz
	if(arrr.lt.80.0) delta=delta+0.5



		ar=(bx-ax)**2+(by-ay)**2+(bz-az)**2
		if(ar.lt.10.0) then
			EHB=EHB+EHBIJ(i,k)*(0.5+delta)
			else
		ar=(bx+ax)**2+(by+ay)**2+(bz+az)**2
		if(ar.lt.10.0) then
			EHB=EHB+EHBIJ(i,k)*(0.5+delta)
			endif
			endif	

		ar=(cx-ax)**2+(cy-ay)**2+(cz-az)**2
		if(ar.lt.10.0) then
			EHB=EHB+EHBIJ(i,k)*(0.5+delta)
			else
		ar=(cx+ax)**2+(cy+ay)**2+(cz+az)**2
		if(ar.lt.10.0) then
			EHB=EHB+EHBIJ(i,k)*(0.5+delta) 
			endif
			endif

		ENDIF
		ENDIF
		ENDIF
		ENDIF
		
	endif		
	endif
	endif
 1002	CONTINUE
	ENDIF

	ENDDO
		
	if(istat.lt.0) then	
        a=0.0
	b=0.0
	if(sumct.gt.0) b= abs(ICNT/float(SUMCT)- acorder) 
        if(sumcto.gt.0) a= abs(ICNTO/float(SUMCTO)- acorder) 
        d= abs(float(SUMCT)-contt)
        c= abs(float(SUMCTO)-contt) 
        dord=ardscal*(b-a)+cstr*(d-c)
        endif	

	RETURN
	END
	
C	***************************************************************	



	FUNCTION ESHORT(iiii,jjjj, ISTAT)
	
	IMPLICIT INTEGER(I-Z)
		PARAMETER(NDIM=400)	
	
	COMMON /CHAINS/  ICA(NDIM), X(NDIM), Y(NDIM), Z(NDIM)
	COMMON /LENGHTS/ LENF1,LENF,PROD(800,800) 
	COMMON /SEQE/  	 SEQ(ndim), SEC(ndim) 	  
	COMMON /VECTORS/ vx(800),vy(800),vz(800)
	COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
	COMMON /HB/      HBX(800,800),HBY(800,800),HBZ(800,800)
	COMMON /BISEC/   CAX(800,800),CAY(800,800),CAZ(800,800)
	COMMON /ENERGY/	 EHBOND, ESC, EREP, BAT, EHBIJ(ndim,ndim)
	COMMON /short/ IBIN(-500:500),asr(ndim,-12:12)
	COMMON /short1/ JBIN(800),bsr(ndim,16) 
	COMMON /short0/ SBIN(200),csr(ndim,8) 
	COMMON /one/ 	acrit, contt, eoinp(0:19,0:100)
	COMMON /shape/	amx,amy,amz	
     	COMMON /FR/ FRG(ndim,ndim)
	COMMON /short2/ erestsum, arestsum,arcs(ndim,ndim)
	COMMON /THREE/   ICONF(800,800)
	COMMON /SG/ GX(800,800,0:19),GY(800,800,0:19),GZ(800,800,0:19)	
	COMMON /RES/ EREST, MRES(ndim),KRES(ndim,50),awrca(ndim,420)
	COMMON /RCN/ ERESTA, MRESA(ndim),KRESA(ndim,420),arca(ndim,420), 
     *		arcam(ndim,420) 
	COMMON /size/	arla(0:19,0:19,2,2),arlm(0:19,0:19,2,2),
     *			arlp(0:19,0:19,2,2)
      common/lim/erlim,arlim,erold,ernew,arold,arnew,actre,ectre 


c	
	     				
	ESHORT=0.0
	ESC2= ESC/2.0
	ESC1= ESC*4.0
	ESC4= ESC/8.0

	if(istat.lt.0) THEN	
	arold=0.0
	erold=0.0	
	else
	ernew=0.0
	arnew=0.0
	endif
	
	iiiii=iiii+1
	jjjjj=jjjj-1
	if(iiii.eq.2) iiiii=2
	if(jjjj.eq.lenf1) jjjjj=lenf1
	
	DO i=iiiii,jjjjj  
	
c	RESTRAINTS - SHORT RANGE & LONG RANGE
c              	
          
        MM=MRESA(i)
        if(MM.GT.0) THEN
        ix=x(i)
        iy=y(i)
        iz=z(i)
        do kkk=1,MM
        j=kresa(i,kkk)
        if(j.lt.i.OR.j.gt.jjjjj) then
        ir=(ix-x(j))**2+(iy-y(j))**2+(iz-z(j))**2
        ar=sqrt(float(ir))
        arc1=arca(i,kkk)
        arc2=arcam(i,kkk)  
			if(ar.LT.arc1) then
			eshort=eshort+eresta*(arc1-ar)*awrca(i,kkk)
			if(istat.lt.0) then
        		arold=arold+(arc1-ar)
        		else
        		arnew=arnew+(arc1-ar)
        		endif
			endif
                        
			if(ar.GT.arc2) then
                        eshort=eshort+eresta*(ar-arc2)*awrca(i,kkk)
                        if(istat.lt.0) then
                        arold=arold+(ar-arc2)
                        else
                        arnew=arnew+(ar-arc2)
                        endif
                        endif

	endif
        enddo
        endif   
        
c	RESTRAINTS - long RANGE (CONTACTS OF Ca-Ca)
c              	
          
        MM=MRES(i)
        if(MM.GT.0) THEN
        ix=x(i)
        iy=y(i)
        iz=z(i)

		ax=ix+gx(ica(i-1),ica(i),seq(i))			
		ay=iy+gy(ica(i-1),ica(i),seq(i))
		az=iz+gz(ica(i-1),ica(i),seq(i))
				
        do kkk=1,MM
        j=kres(i,kkk)
        if(j.lt.i.OR.j.gt.jjjjj) then
	
		bx=x(j)+gx(ica(j-1),ica(j),seq(j))
		by=y(j)+gy(ica(j-1),ica(j),seq(j))
		bz=z(j)+gz(ica(j-1),ica(j),seq(j))

		ar=(ax-bx)**2+(ay-by)**2+(az-bz)**2

		ar=sqrt(ar)-arla(seq(i),seq(j),1,1)
			if(ar.gt.4.0) then
			eshort=eshort+erest*sqrt(ar)
			ar=ar-4.0
			if(istat.lt.0) then
        		erold=erold+ar
        		else
        		ernew=ernew+ar
        		endif
			endif
        endif
        enddo
        endif   

	ENDDO

        do i=iiii,jjjj
        k=seq(i)
        ii=ICA(i-1)
        jj=ICA(i)
        
        aract2=float(x(i)**2+y(i)**2+z(i)**2)+0.01
        aract=sqrt(aract2)

        ff=aract/acrit
        ika=int(ff/0.333333)
c	one-body centrosymmetric
        ESHORT=ESHORT+eoinp(k,ika)


c	and the r13 potentials
	kkk=iconf(ii,jj)
	ESHORT=ESHORT+csr(i-1,SBIN(kkk))
	enddo
               
	ii1=iiii-4
	if(ii1.lt.2) ii1=2
	ii2=jjjj+1
	if(ii2.gt.lenf1-4) ii2=lenf1-4
	do i=ii1,ii2
	i1=i+1
	i2=i+2
	i3=i+3
	j=i+4
	icam4=ica(i-1)
	icam3=ica(i)
	icam2=ica(i1)
	icam1=ica(i2)
	icam=ica(i3)
	icaj=ica(j)

c	ff=1.0
c	if(sec(i+1).eq.1) ff=ff-0.33
c	if(sec(i+2).eq.1) ff=ff-0.33 
c	if(sec(i+3).eq.1) ff=ff-0.33

	ax=float(x(i1)+x(i2)+x(i3)+x(j))/4.0
	ay=float(y(i1)+y(i2)+y(i3)+y(j))/4.0
	az=float(z(i1)+z(i2)+z(i3)+z(j))/4.0	
	ff=min(1.0,((acrit*acrit)/(ax*ax+ay*ay+az*az+0.01)))
		
	
c
c	Modiffied the generic stiffness algorithm (5/09/03)

c
	a=CAX(icam3,icam2)+CAX(icam2,icam1)+CAX(icam1,icam)+CAX(icam,icaj) 
	b=CAY(icam3,icam2)+CAY(icam2,icam1)+CAY(icam1,icam)+CAY(icam,icaj)
	c=CAZ(icam3,icam2)+CAZ(icam2,icam1)+CAZ(icam1,icam)+CAZ(icam,icaj)
	aor4=max(sqrt(a*a+b*b+c*c), 0.5) -0.5
	if(aor4.gt.0.5) aor4=0.5
	ESHORT=ESHORT-ff*ESC*(0.5-aor4)

c
c	Additional generic (13/07/05)
c
	a=CAX(icam3,icam2)*CAX(icam2,icam1) 
	b=CAY(icam3,icam2)*CAY(icam2,icam1)
	c=CAZ(icam3,icam2)*CAZ(icam2,icam1)
	aor1=max(0.0, (a+b+c))

	a=CAX(icam2,icam1)*CAX(icam1,icam) 
	b=CAY(icam2,icam1)*CAY(icam1,icam)
	c=CAZ(icam2,icam1)*CAZ(icam1,icam)
	aor2=max(0.0, (a+b+c))

	ESHORT=ESHORT +ESC*(aor1+aor2)


     		r14=(x(i3)-x(i))**2+(y(i3)-y(i))**2+(z(i3)-z(i))**2
     		WX1=VX(icam3)
		WY1=VY(icam3)
		WZ1=VZ(icam3)
		WX2=VX(icam2)
		WY2=VY(icam2)
		WZ2=VZ(icam2)
		WX3=VX(icam1)
		WY3=VY(icam1)
		WZ3=VZ(icam1)
		px=wy1*wz2-wy2*wz1
		py=wz1*wx2-wz2*wx1
		pz=wx1*wy2-wx2*wy1
		ihand=px*wx3+py*wy3+pz*wz3	
		if(ihand.lt.0) 	r14=-r14
		IBIN4=IBIN(r14)
c	contribution from chiral r14 potential
 
	 ESHORT=ESHORT+asr(i,IBIN4)
     
	ix=x(j)-x(i)
	iy=y(j)-y(i)
	iz=z(j)-z(i)	
		 
	iii=ix**2+iy**2+iz**2
c	contribution from r15 potential
 	
	ESHORT=ESHORT+bsr(i,JBIN(iii))

c
c	GENERIC BIAS TOWARDS (***) or (###)
c
	
	
	if(iii.gt.80.AND.iii.lt.160) THEN

	if(sec(i+1).ne.4.AND.sec(i+2).ne.4.AND.sec(i+3).ne.4) then
	if(sec(i).ne.4.AND.sec(j).ne.4) then

	a=CAX(icam3,icam2)*CAX(icam1,icam) 
	b=CAY(icam3,icam2)*CAY(icam1,icam)
	c=CAZ(icam3,icam2)*CAZ(icam1,icam)
	aor2=a+b+c		
	if(aor2.lt.0.0) then

c	a helix (***)		

	if(IBIN4.gt.4.AND.IBIN4.lt.8) ESHORT=ESHORT-ESC-ff*ESC2
	endif
	endif
	endif
	
	ELSE
c	a beta (###)
	if(iii.gt.360.AND.iii.lt.530) 	then

	if(sec(i+1).ne.2.AND.sec(i+2).ne.2.AND.sec(i+3).ne.2) then
	if(sec(i).ne.2.AND.sec(j).ne.2) then
	
	a=CAX(icam3,icam2)*CAX(icam1,icam) 
	b=CAY(icam3,icam2)*CAY(icam1,icam)
	c=CAZ(icam3,icam2)*CAZ(icam1,icam)
	aor2=a+b+c		
	if(aor2.gt.0.5) then


	a=CAX(icam3,icam2)*CAX(icam2,icam1) 
	b=CAY(icam3,icam2)*CAY(icam2,icam1)
	c=CAZ(icam3,icam2)*CAZ(icam2,icam1)
	aor3=a+b+c		
	if(aor3.lt.0.0) then
	a=CAX(icam2,icam1)*CAX(icam1,icam) 
	b=CAY(icam2,icam1)*CAY(icam1,icam)
	c=CAZ(icam2,icam1)*CAZ(icam1,icam)
	aor4=a+b+c		
	if(aor4.lt.0.0) then
	if(iabs(IBIN4).gt.8) ESHORT=ESHORT-ESC-ff*ESC2
	endif
	endif
	endif
	endif
	endif
	endif
				
	ENDIF
	
			
	enddo
	
C	***************************************************************	
c	THE OLD STUFF, MODIFIED FOR THE NEW LATTICE	
c	

	I1=IIII-14
	if(i1.lt.2) i1=2
	i2=JJJJ
	if(i2.gt.lenf1-15) i2=lenf1-15
	
	do i=i1,i2
		j=i+5
		k=i+10
		wx1=x(j)-x(i)
		wy1=y(j)-y(i)
		wz1=z(j)-z(i)
		wx2=x(k)-x(j)
		wy2=y(k)-y(j)
		wz2=z(k)-z(j)
	if(wx1*wx2+wy1*wy2+wz1*wz2.lt.0) THEN
		kk=k+5
		wx3=x(kk)-x(k)
		wy3=y(kk)-y(k)
		wz3=z(kk)-z(k)	
		if(wx1*wx3+wy1*wy3+wz1*wz3.gt.0) ESHORT=ESHORT+ESC1
c	penalty for "crumpling"
	ENDIF
	enddo
	


	
C	***************************************************************	
c	bias to predicted secondary structure geometry
		
	I1=IIII-7
	if(i1.lt.2) i1=2
	i2=JJJJ
	if(i2.gt.lenf1-7) i2=lenf1-7
	
	do i=i1,i2
	k=i+7
	m=k-6
	if(sec(m+3).ne.2) then	
    	ag=FRG(m,k)
 	if(ag.gt.0.01) then
 	ff=float((x(m)-x(k))**2+(y(m)-y(k))**2+(z(m)-z(k))**2)
	gg=sqrt(ff)
	df=abs(ag-gg)
	if(df.gt.2.0) ESHORT=ESHORT+ESC4*(df-1.0)+ESC2
        endif
        endif
        
        m=k-7
	if(sec(m+3).ne.4) then	
        ag=FRG(m,k)
        if(ag.gt.0.01) then
	ff=float((x(m)-x(k))**2+(y(m)-y(k))**2+(z(m)-z(k))**2)
	gg=sqrt(ff)
	df=abs(ag-gg)
	if(df.gt.1.0) ESHORT=ESHORT+ESC4*df+ESC2
        endif 
        endif 

        enddo

	IF(istat.eq.10) then
c	The INITIAL (and final) energy of the restraints
	if(ernew.gt.erlim) ESHORT=ESHORT+erest*(ernew-erlim)
	if(arnew.gt.arlim) ESHORT=ESHORT+eresta*(arnew-arlim)
	endif
	
	if(istat.lt.0) then
	ectre=erestsum+ernew-erold
	actre=arestsum+arnew-arold
c	the old energy of restraints (opposite signs due to -EOLD)        
	if(erestsum.gt.erlim) ESHORT=ESHORT+erest*(erestsum-erlim)
	if(arestsum.gt.arlim) ESHORT=ESHORT+eresta*(arestsum-arlim)

c	the new energy of restraints	
	if(ectre.gt.erlim) ESHORT=ESHORT-erest*(ectre-erlim)
	if(actre.gt.arlim) ESHORT=ESHORT-eresta*(actre-arlim)			
	endif
		 	                       

	RETURN
	END
		
	
C	***************************************************************	
	
