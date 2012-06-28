**==aa0001.spg  processed by SPAG 5.11R  at 19:34 on 18 Apr 2002
c
c Chain building for the CABS model using the aligned residues of
c      a template given in a PDB file (CAtempal.pdb)
c	Modified by A. Szilagyi from the build_chain.f program for
c      the SICHO model by A. Kolinski
c  Some errors in the original program were corrected
c      by A. Szilagyi on 04/18/2002
 
c Original program says:
c		Caution:  The conformations could be trapped
c		CHECK THE PROPER EXIT FROM THE PROGRAM- OUT
c
 
c   g77-3.3 -O3 build_cabs61.f 
      implicit none
C*** Start of declarations inserted by SPAG
      real aax , aay , aaz , abox2 , al2 , ar , ar1 , ar2 , armin , 
     &     asx , asy , asz , ax , axi , ay , ayi , az , azi , bx , by
      real bz , centx , centy , centz , dd , dx , dy , dz , f
      integer i , Ica , Iconf , id , idel , ii , ip , ir , ix , iy , 
     &        iz , j , jj , jp , jstart , jx , jy , jz , k , kk
      integer kkk , kx , ky , kz , Lenf , Lenf1 , Lenf2 , Length2 , 
     &        nal , NBOX , NDIM , nmin , nwmax , NWMAX0 , res1 , res2 , 
     &  text , vector , Vx , Vy, lstart,lend,l,ll,NLOC,
     &  resx1, resx2
      integer Vz , wx , wy , wz , X , Y , Z,zflag
C*** End of declarations inserted by SPAG
 
      parameter (NBOX=600)
      parameter (NDIM=1000)
      parameter (NLOC=10000)
      parameter (NWMAX0=800)

      logical LOOK , goodc(NWMAX0,NWMAX0) , dummy(NDIM)

c instead of a huge box that occupies a lot of memory,
c here let's use a list of (x,y,z,m) quadruplets
      
      common /WORK  / Xs(NLOC),Ys(NLOC),Zs(NLOC),Ms(NLOC),Nset
      integer*2 Xs,Ys,Zs,Ms,Nset
      common /CHAIN / Ica(-2:NDIM) , X(-2:NDIM) , Y(-2:NDIM) , 
     &                Z(-2:NDIM)
 
      common /THREE / Iconf(NWMAX0,NWMAX0)
      common /LENGHTS/ Lenf2 , Lenf1 , Lenf , Length2(NWMAX0)
 
      common Vx(NWMAX0) , Vy(NWMAX0) , Vz(NWMAX0)
      dimension vector(-6:6,-6:6,-6:6) , ax(NDIM) , ay(NDIM) , az(NDIM)
        
c     MJamroz - read command line arguments
c     lngh: 1 - target len, 2 - template len, 3 - model number 
      CHARACTER*4 tmptxt
      CHARACTER*13 chainout,alignin

      IF (iargc().NE.3) THEN
        write(*,*) "Usage: <target atm> <template atm> <mdl idx>"
      END IF
      CALL getarg(3,tmptxt)
      chainout='CHAIN'//tmptxt
      alignin='ALIGN'//tmptxt
      CALL getarg(1,tmptxt)
      read(tmptxt,'(i4)' ) Lenf2
      CALL getarg(2,tmptxt)
      read(tmptxt,'(i4)' ) nal
      resx1=1
      resx2=Lenf2 
      Nset=0
C	     ***********************************************************
C	    PREPARATION OF THE CHAIN UNITS AND THEIR CORRELATION
      nwmax = 0
 
      do 100 ix = -7 , 7
        do 50 iy = -7 , 7
          do 20 iz = -7 , 7
               vector(ix,iy,iz) = 0
               ir = ix*ix + iy*iy + iz*iz
               if ( ir.ge.29 ) then
                 if ( ir.le.49) then
                     nwmax = nwmax + 1
                     Vx(nwmax) = ix
                     Vy(nwmax) = iy
                     Vz(nwmax) = iz
                     vector(ix,iy,iz) = nwmax
                  endif
               endif
 20         continue
 50      continue
 100  continue
	
c
       

 
      open (unit=5,file=alignin,status='OLD')
      open (unit=10,file=chainout,status='UNKNOWN')
      open (unit=11,file='PDB_'//chainout,status='UNKNOWN')
      rewind (5)
      rewind(11)
      rewind (10)
	
c

      do 200 i = 1 , NWMAX0
         ix = Vx(i)
         iy = Vy(i)
         iz = Vz(i)
         do 150 j = 1 , NWMAX0
            jx = Vx(j)
            jy = Vy(j)
            jz = Vz(j)
            Iconf(i,j) = (ix+jx)**2 + (iy+jy)**2 + (iz+jz)**2
 150     continue
 200  continue
 
C	..............................................................
C
      do 300 i = 1 , NWMAX0
         ix = Vx(i)
         iy = Vy(i)
         iz = Vz(i)
         do 250 j = 1 , NWMAX0
           goodc(i,j)=.true.
           if ( Iconf(i,j).lt.45) goodc(i,j) = .false.
           if ( Iconf(i,j).gt.145) goodc(i,j) = .false.
               jx = Vx(j)
               jy = Vy(j)
               jz = Vz(j)
c
               kx = iy*jz - iz*jy
               ky = jx*iz - ix*jz
               kz = ix*jy - iy*jx
c
               kkk = kx*kx + ky*ky + kz*kz
c
c	correction for 440 and 330 pairs of vectors (colinear)
c
               if ( kkk.eq.0 ) goodc(i,j) = .false.
 
 250     continue
 300  continue
 
 
c
C       INPUT  INPUT   INPUT   INPUT   INPUT  INPUT  INPUT
C	--------------------------------------------------
C	 SET UP OF THE VECTOR REPRESENTATION OF THE CHAIN
C
 
c      do 400 i = 1 , NBOX
c         do 350 k = 1 , NBOX
c            do 320 j = 1 , NBOX
c               Xyz(j,k,i) = 0
c 320        continue
c 350     continue
c 400  continue
	
      abox2 = NBOX/2
	
      bx = 0.
      by = 0.
      bz = 0.
 
c     read (5,*) Lenf2 , nal , resx1 , resx2


      
      res1=1
      res2=Lenf2
      Lenf = Lenf2 + 2
      Lenf1 = Lenf - 1
      al2 = Lenf2
 
      do 500 i = 2 , Lenf1
         dummy(i) = .false.
 500  continue
      dummy(1) = .true.
      dummy(Lenf) = .true.
 
99001 format (a22,i4,4x,3F8.3)
 
 
      read (5,99001) text , ii , aax , aay , aaz
	
      jp = ii
      jstart = jp + 1
	
      rewind (5)
c      read (5,*)
	
      do 600 i = 1 , nal
         read (5,99001) text , ii , aax , aay , aaz
         j = ii + 1
	
         ax(j) = aax
         ay(j) = aay
         az(j) = aaz
c
c	THERE IS THE INSERTION PART, NEEDS TO BE TREATED CAREFULLY
c	( build "a bubble" outside the main core of the alignment)
c
         if ( j.ne.jp+1 ) then
            kk = j - jp
            dx = (ax(j)/0.61-ax(jp))/kk
            dy = (ay(j)/0.61-ay(jp))/kk
            dz = (az(j)/0.61-az(jp))/kk
	
            do 520 k = jp + 1 , ii
               ax(k) = ax(k-1) + dx
               ay(k) = ay(k-1) + dy
               az(k) = az(k-1) + dz
               dummy(k) = .true.
 520        continue
         endif
	
         ax(j) = ax(j)/0.61
         ay(j) = ay(j)/0.61
         az(j) = az(j)/0.61
         bx = bx + ax(j)
         by = by + ay(j)
         bz = bz + az(j)
         jp = j
 600  continue
	
      centx = bx/nal
      centy = by/nal
      centz = bz/nal
	
      asx = abox2 - bx/nal
      asy = abox2 - by/nal
      asz = abox2 - bz/nal
	
      do 700 i = 2 , Lenf1
         ax(i) = ax(i) + asx
         ay(i) = ay(i) + asy
         az(i) = az(i) + asz
 700  continue

      write(*,*) 'centxyz=',centx,centy,centz
c
c	ADD THE ENDS IF NECESSARY
c
      if ( res1.eq.1 ) then
        if ( jstart.ne.2 ) then
          idel = jstart - res1
c 3.5 seems to be a quite good value           
          dd = 3.5*MIN(150,idel)/FLOAT(idel)
c            dd=2.0
          asx = ax(jstart) - abox2
          asy = ay(jstart) - abox2
          asz = az(jstart) - abox2
          ar = SQRT(asx**2+asy**2+asz**2+0.0001)
          asx = asx/ar
          asy = asy/ar
          asz = asz/ar
          write(*,*) 'asx,asy,asz=',asx,asy,asz
          
          do 720 k = jstart , 2 , -1
            ax(k-1) = ax(k) + asx*dd
            ay(k-1) = ay(k) + asy*dd
            az(k-1) = az(k) + asz*dd
            dummy(k-1) = .true.
 720      continue
        endif
      endif
 
c
c	SECOND END (C-terminus)  WHEN NEEDED
c
 
      if ( res2.eq.Lenf2 ) then
         if ( j.lt.Lenf1 ) then
            idel = res2 + 1 - j
            dd = 3.5*MIN(150,idel)/FLOAT(idel)
c            dd=2.0
            asx = ax(j) - abox2
            asy = ay(j) - abox2
            asz = az(j) - abox2
            ar = SQRT(asx**2+asy**2+asz**2+0.0001)
            asx = asx/ar
            asy = asy/ar
            asz = asz/ar
            write(*,*) 'asx,asy,asz=',asx,asy,asz
            do 740 k = j , Lenf2
               ax(k+1) = ax(k) + asx*dd
               ay(k+1) = ay(k) + asy*dd
               az(k+1) = az(k) + asz*dd
               dummy(k+1) = .true.
 740        continue
         endif
      endif

	
c
c	MAKE THE LOOPS CONVEX in RESPECT to the TEMPLATE
c             old algorithm, not very good
c      do 800 k = 1 , 5
c         do 750 j = 4 , Lenf - 3
c            if ( dummy(j) ) then
c               ar1 = (ax(j)-ax(j-1))**2 + (ay(j)-ay(j-1))
c     &           **2 + (az(j)-az(j-1))**2
c               ar2 = (ax(j)-ax(j+1))**2 + (ay(j)-ay(j+1))
c     &           **2 + (az(j)-az(j+1))**2
c               if ( (ar1+ar2).lt.16.0 ) then
c                  f = (ax(j-1)-ax(j+1))**2 + (ay(j-1)-ay(j+1))
c     &                **2 + (az(j-1)-az(j+1))**2
c                  if ( f.lt.16.0 ) then
c                    asx = ax(j) - abox2
c                    asy = ay(j) - abox2
c                    asz = az(j) - abox2
c                    asx = asx/2.
c                    asy = asy/2.
c                    asz = asz/2.
c                    ax(j) = ax(j) + asx/2.0
c                    ay(j) = ay(j) + asy/2.0
c                    az(j) = az(j) + asz/2.0
c                    write(*,*) 'convexization ',j,asx/2.0,asy/2.0,asz/2.0
c                  endif
c               endif
c            endif
c 750     continue
c 800  continue

c     MAKE THE LOOPS CONVEX WITH RESPECT TO THE TEMPLATE
c        new algorithm by A. Szilagyi, hopefully better
c rewritten again on 5/9/2002 and 5/13/2002

      l=jstart
      do 750 while (l.lt.j)
        do while (.not.dummy(l) .and. l.lt.j)
          l=l+1
        enddo
        if (l.ge.j) goto 800
        lstart=l
        do while (dummy(l) .and. l.lt.j)
          l=l+1
        enddo
        if (l.ge.j) goto 800
        lend=l-1
        ll=lend-lstart+1
        write(*,*)'loop:',lstart,lend

c this goto skips the next section
        
        goto 749

c adjust the loop stems if seems necessary
c this section is optional, doesn't seem to help much
c but maybe in some cases it does        
        if (lstart.gt.2) then
          aax=ax(lstart-1)-abox2
          aay=ay(lstart-1)-abox2
          aaz=az(lstart-1)-abox2
          asx=ax(lstart-2)-abox2
          asy=ay(lstart-2)-abox2
          asz=az(lstart-2)-abox2
          ar=SQRT(aax*aax+aay*aay+aaz*aaz)
          ar1=SQRT(asx*asx+asy*asy+asz*asz)
          if (ar .lt. ar1) then
            ax(lstart-1)=abox2+aax*ar1*ar1/ar/ar
            ay(lstart-1)=abox2+aay*ar1*ar1/ar/ar
            az(lstart-1)=abox2+aaz*ar1*ar1/ar/ar
            write(*,*)'modified loop start:',lstart-1
          endif
        endif
        if (lend.lt.j-2) then
          aax=ax(lend+1)-abox2
          aay=ay(lend+1)-abox2
          aaz=az(lend+1)-abox2
          asx=ax(lend+2)-abox2
          asy=ay(lend+2)-abox2
          asz=az(lend+2)-abox2
          ar=SQRT(aax*aax+aay*aay+aaz*aaz)
          ar1=SQRT(asx*asx+asy*asy+asz*asz)
          if (ar .lt. ar1) then
            ax(lend+1)=abox2+aax*ar1*ar1/ar/ar
            ay(lend+1)=abox2+aay*ar1*ar1/ar/ar
            az(lend+1)=abox2+aaz*ar1*ar1/ar/ar
            write(*,*)'modified loop end:',lend+1
          endif
        endif
        asx=(ax(lend+1)-ax(lstart-1))/(ll+1)
        asy=(ay(lend+1)-ay(lstart-1))/(ll+1)
        asz=(az(lend+1)-az(lstart-1))/(ll+1)
        do l=lstart,lend
          ax(l)=ax(lstart-1)+(l-lstart+1)*asx
          ay(l)=ay(lstart-1)+(l-lstart+1)*asy
          az(l)=az(lstart-1)+(l-lstart+1)*asz
        enddo

c now make loop 'convex'
        
 749    do l=lstart,lend
          asx=ax(l)-abox2
          asy=ay(l)-abox2
          asz=az(l)-abox2
          ar=SQRT(asx**2+asy**2+asz**2)
c 3.5 appears to be a really good value          
          asx=3.5*asx/ar
          asy=3.5*asy/ar
          asz=3.5*asz/ar
          dd=min(l-lstart+1,lend+1-l)
          ax(l)=ax(l)+asx*dd
          ay(l)=ay(l)+asy*dd
          az(l)=az(l)+asz*dd
        enddo
        l=lend+1
 750  continue
        
        

c
c write out pdb file with real coordinates
c (dummy coordinates as well)
c all residue names set to GLY for simplicity      

 800  do 850 i=resx1+1,resx2+1
c        if (.not. dummy(i)) then
        write (11,99004) i-1,ax(i)*0.61,ay(i)*0.61,az(i)*0.61
c        endif
 850  continue
99004 format('ATOM     1   CA  GLY A',I4,'    ',3F8.3)

c
c	SIMULATE THE EXCLUDED VOLUME OF THE ALIGNED PART
c
      do 900 i = 3 , Lenf1
         if ( .not.dummy(i) ) then
            ix = NINT(ax(i))
            iy = NINT(ay(i))
            iz = NINT(az(i))
            if ( .not.LOOK(ix,iy,iz,i) ) call SET(ix,iy,iz,i)
         endif
 900  continue

c
c	Generate the second entry of CHAIN
c
      ix = NINT(ax(2))
      iy = NINT(ay(2))
      iz = NINT(az(2))
      X(2) = ix
      Y(2) = iy
      Z(2) = iz
      if ( .not.LOOK(ix,iy,iz,2) ) call SET(ix,iy,iz,2)

c
c	Generate the third entry of CHAIN
c
      nmin = 0
      armin = 1000000.
	
      if ( .not.dummy(3) ) then
         jx = NINT(ax(3))
         jy = NINT(ay(3))
         jz = NINT(az(3))
         call REM(jx,jy,jz)
      endif
	
      do 1000 k = 1 , nwmax
         jx = Vx(k)
         jy = Vy(k)
         jz = Vz(k)
         kx = ix + jx
         ky = iy + jy
         kz = iz + jz
         if ( .not.(LOOK(kx,ky,kz,3)) ) then
            ar = (kx-ax(3))**2 + (ky-ay(3))**2 + (kz-az(3))**2
            if ( ar.lt.armin ) then
               nmin = k
               armin = ar
            endif
         endif
 1000 continue

      if (nmin.ne.0) then
        kx = ix + Vx(nmin)
        ky = iy + Vy(nmin)
        kz = iz + Vz(nmin)
        call SET(kx,ky,kz,3)
        X(3) = kx
        Y(3) = ky
        Z(3) = kz
      endif
      if ( nmin.eq.0 ) then
         write (10,*) '01 error on the bead # 3'
         stop
      endif
c
c	Generate the middle part of CHAIN
c
      do 1200 i = 4 , Lenf1
	
         if ( .not.dummy(i) ) then
            ix = NINT(ax(i))
            iy = NINT(ay(i))
            iz = NINT(az(i))
            call REM(ix,iy,iz)
         endif
	
 1050    continue
         ix = X(i-1)
         iy = Y(i-1)
         iz = Z(i-1)
         ip = nmin
         nmin = 0
         armin = 1000000.
         axi = ax(i)
         ayi = ay(i)
         azi = az(i)
         zflag=1
         do 1100 k = 1 , nwmax
            jx = Vx(k)
            jy = Vy(k)
            jz = Vz(k)
            if ( goodc(ip,k) ) then
              kx = ix + jx
              ky = iy + jy
              kz = iz + jz
              if ( .not.(LOOK(kx,ky,kz,i)) ) then
                ar = (kx-axi)**2 + (ky-ayi)**2 + (kz-azi)**2
                if ( ar.lt.armin ) then
                  nmin = k
                  armin = ar
                endif
              endif
            endif
 1100     continue
 
         if ( nmin.eq.0 ) then
            do 1120 id = i + 1 , Lenf1
               if ( .not.dummy(id) ) then
                  ix = NINT(ax(id))
                  iy = NINT(ay(id))
                  iz = NINT(az(id))
                  call REM(ix,iy,iz)
                  dummy(id) = .true.
                  goto 1050
               endif
 1120       continue
         endif
 
         if ( nmin.eq.0 ) then
           write (10,*) '02 error on the bead # ' , i
           stop
         endif
c         write(6,*)'i,Vx,Vy,Vz,nmin=',i,Vx(nmin),Vy(nmin),Vz(nmin),nmin
         kx = ix + Vx(nmin)
         ky = iy + Vy(nmin)
         kz = iz + Vz(nmin)
         
         call SET(kx,ky,kz,i)
         X(i) = kx
         Y(i) = ky
         Z(i) = kz
 1200 continue
c
c	Generate the last entry of CHAIN
c
      ix = kx
      iy = ky
      iz = kz
      ip = nmin
      nmin = 0
 
      do 1300 k = 1 , nwmax
         jx = Vx(k)
         jy = Vy(k)
         jz = Vz(k)
         if ( goodc(ip,k) ) then
            kx = ix + jx
            ky = iy + jy
            kz = iz + jz
            if ( .not.(LOOK(kx,ky,kz,Lenf)) ) then
               nmin = k
               goto 1400
            endif
         endif
 
 1300 continue
 
 1400 continue

      if ( nmin.eq.0 ) then
         write (10,*) '03 error on the bead # LENF'
         stop
       endif
       
       kx = ix + Vx(nmin)
       ky = iy + Vy(nmin)
       kz = iz + Vz(nmin)
       call SET(kx,ky,kz,Lenf)
       X(Lenf) = kx
       Y(Lenf) = ky
       Z(Lenf) = kz

c
c	Generate the first entry of CHAIN
c
      ix = X(2)
      iy = Y(2)
      iz = Z(2)
      jx = X(2) - X(3)
      jy = Y(2) - Y(3)
      jz = Z(2) - Z(3)
      ip = vector(jx,jy,jz)
      write (*,*) "ip=" , ip
      nmin = 0
 
      do 1500 k = 1 , nwmax
         jx = Vx(k)
         jy = Vy(k)
         jz = Vz(k)
         if ( goodc(ip,k) ) then
            kx = ix - jx
            ky = iy - jy
            kz = iz - jz
            if ( .not.(LOOK(kx,ky,kz,1)) ) then
               nmin = k
               goto 1600
            endif
         endif
 
 1500 continue
 
 1600 continue
      if ( nmin.eq.0 ) then
        write (10,*) '04 error on the bead # 1'
        stop
      endif
      
      kx = ix - Vx(nmin)
      ky = iy - Vy(nmin)
      kz = iz - Vz(nmin)
      call SET(kx,ky,kz,1)
      X(1) = kx
      Y(1) = ky
      Z(1) = kz

*** write out the CHAIN file ***	
1650  write (10,*) resx2-resx1+3
      do 1700 i = resx1 , resx2+2
        write (10,99002) X(i) , Y(i) , Z(i) , i
 1700 continue
99002 format (4I4)

      close(10)
c*** write out PDB format file but lattice model!
c      write(11,"('TER')")
c 1750 do 1755 i=resx1,resx2+2
c        write (11,99005) i-1,X(i)*0.61,Y(i)*0.61,Z(i)*0.61
c 1755 continue
c      write(11,"('TER')")
c99005 format('ATOM     1   CA  GLY B',I4,'    ',3F8.3)
      
*      stop
c      write(6,*)'i X Y Z Ica Iconf ip'
c      do 1800 i = 1 , Lenf1
c        j = i + 1
c        write(6,*)'i=',i,' j=',j
c        write(6,*)'X(i)=',X(i)
c        write(6,*)'Y(i)=',Y(i)
c        write(6,*)'Z(i)=',Z(i)
c        write(6,*)'X(j)=',X(j)
c        write(6,*)'Y(j)=',Y(j)
c        write(6,*)'Z(j)=',Z(j)
c        wx = X(j) - X(i)
c        wy = Y(j) - Y(i)
c        wz = Z(j) - Z(i)
c        Ica(i) = vector(wx,wy,wz)
c        write(6,*)'Ica(i)=',Ica(i)
c         if ( i.gt.1 ) then
c           jj = Ica(i-1)
c           write(6,*)'Ica(i-1)=jj=',jj
c           write(6,*)'Iconf(jj,Ica(i))',Iconf(jj,Ica(i))
c           ip = wx*wx + wy*wy + wz*wz
c           write(6,*)'ip=',ip,' sqrt(ip)=',SQRT(float(ip))
cc            write (6,99003) i , X(i) , Y(i) , Z(i) , Ica(i) , 
cc     &                      Iconf(jj,Ica(i)) , ip
c         endif
	
c 1800  continue
       
c 99003    format (7I4)

c      write(6,*)'Vectors:'
c      do 1900 j = 1 , Lenf1
c         jj = Ica(j)
c         write (6,*) j , Vx(jj) , Vy(jj) , Vz(jj)
c 1900 continue

c      write(6,*)'Ca-Ca overlaps:'
c      do i=1,Lenf1-2
c        do j=i+2,Lenf1
c          wx=X(j)-X(i)
c          wy=Y(j)-Y(i)
c          wz=Z(j)-Z(i)
c          ip=wx*wx+wy*wy+wz*wz
c          if (ip .lt. 16) then
c            write(6,*) i,j,ip
c          endif
c        enddo
c      enddo
      stop
      end
**==look.spg  processed by SPAG 5.11R  at 19:34 on 18 Apr 2002
c
C	**************************************************************
	
	
      function LOOK(I,J,K,rn)
      implicit none
C*** Start of declarations inserted by SPAG
      integer I , Ica , J , K , NBOX , NDIM , X , Y , Z,NLOC
      integer mx,my,mz,rn
      integer*2 pfind,s
C*** End of declarations inserted by SPAG
	
      parameter (NDIM=1000)
      parameter (NLOC=10000)
	
      common /CHAIN / Ica(-2:NDIM) , X(-2:NDIM) , Y(-2:NDIM) , 
     &                Z(-2:NDIM)
      logical LOOK
      common /WORK  / Xs(NLOC),Ys(NLOC),Zs(NLOC),Ms(NLOC),Nset
      integer*2 Xs,Ys,Zs,Ms,Nset
 
 
      LOOK = .true.
 
c      if ( Xyz(I,J,K).gt.0 ) return
c 
c      if ( Xyz(I+1,J,K).gt.0 ) return
c      if ( Xyz(I-1,J,K).gt.0 ) return
c      if ( Xyz(I,J+1,K).gt.0 ) return
c      if ( Xyz(I,J-1,K).gt.0 ) return
c      if ( Xyz(I,J,K+1).gt.0 ) return
c      if ( Xyz(I,J,K-1).gt.0 ) return
c 
c      if ( Xyz(I+1,J+1,K).gt.0 ) return
c      if ( Xyz(I+1,J-1,K).gt.0 ) return
c      if ( Xyz(I-1,J+1,K).gt.0 ) return
c      if ( Xyz(I-1,J-1,K).gt.0 ) return
c 
c      if ( Xyz(I,J+1,K+1).gt.0 ) return
c      if ( Xyz(I,J+1,K-1).gt.0 ) return
c      if ( Xyz(I,J-1,K+1).gt.0 ) return
c      if ( Xyz(I,J-1,K-1).gt.0 ) return
c 
c      if ( Xyz(I+1,J,K+1).gt.0 ) return
c      if ( Xyz(I+1,J,K-1).gt.0 ) return
c      if ( Xyz(I-1,J,K+1).gt.0 ) return
c      if ( Xyz(I-1,J,K-1).gt.0 ) return
c
      DO mx=-4,4
        DO my=-4,4
          DO mz=-4,4
            s=pfind(I+mx,J+my,K+mz)
            if (mx*mx + my*my + mz*mz .lt. 32) then
              if (s.gt.0 .and. s.ne.rn-1 .and. s.ne.rn+1) return
            endif
          enddo
        enddo
      enddo
        
      LOOK = .false.
 
      return
      end
**==set.spg  processed by SPAG 5.11R  at 19:34 on 18 Apr 2002
	
	
C	**************************************************************
	
      subroutine SET(I,J,K,M)
      implicit none
C*** Start of declarations inserted by SPAG
      integer I , Ica , J , K , M , NBOX , NDIM , NWMAX0 , X , Y , Z
      integer NLOC
C*** End of declarations inserted by SPAG
	
      parameter (NLOC=10000)
      parameter (NDIM=1000)
      parameter (NWMAX0=800)
      common /CHAIN / Ica(-2:NDIM) , X(-2:NDIM) , Y(-2:NDIM) , 
     &                Z(-2:NDIM)
      common /WORK  / Xs(NLOC),Ys(NLOC),Zs(NLOC),Ms(NLOC),Nset
      integer*2 Xs,Ys,Zs,Ms,Nset
 
      call padd(I,J,K,M)
 
      call padd(I+1,J,K,M)
      call padd(I-1,J,K,M)
      call padd(I,J+1,K,M)
      call padd(I,J-1,K,M)
      call padd(I,J,K+1,M)
      call padd(I,J,K-1,M)
 
      call padd(I+1,J+1,K,M)
      call padd(I+1,J-1,K,M)
      call padd(I-1,J+1,K,M)
      call padd(I-1,J-1,K,M)
 
      call padd(I,J+1,K+1,M)
      call padd(I,J+1,K-1,M)
      call padd(I,J-1,K+1,M)
      call padd(I,J-1,K-1,M)
 
      call padd(I+1,J,K+1,M)
      call padd(I+1,J,K-1,M)
      call padd(I-1,J,K+1,M)
      call padd(I-1,J,K-1,M)
 
 
      return
      end
**==rem.spg  processed by SPAG 5.11R  at 19:34 on 18 Apr 2002
	
C	**************************************************************
	
      subroutine REM(I,J,K)
      implicit none
C*** Start of declarations inserted by SPAG
      integer I , Ica , J , K , NBOX , NDIM , X , Y , Z,NLOC
C*** End of declarations inserted by SPAG
	
      parameter (NLOC=10000)
      parameter (NDIM=1000)
      common /CHAIN / Ica(-2:NDIM) , X(-2:NDIM) , Y(-2:NDIM) , 
     &                Z(-2:NDIM)
      common /WORK  / Xs(NLOC),Ys(NLOC),Zs(NLOC),Ms(NLOC),Nset
      integer*2 Xs,Ys,Zs,Ms,Nset
 
      call prem(I,J,K) 
 
      call prem(I+1,J,K) 
      call prem(I-1,J,K) 
      call prem(I,J+1,K) 
      call prem(I,J-1,K) 
      call prem(I,J,K+1) 
      call prem(I,J,K-1) 
 
      call prem(I+1,J+1,K) 
      call prem(I+1,J-1,K) 
      call prem(I-1,J+1,K) 
      call prem(I-1,J-1,K) 
 
      call prem(I,J+1,K+1) 
      call prem(I,J+1,K-1) 
      call prem(I,J-1,K+1) 
      call prem(I,J-1,K-1) 
 
      call prem(I+1,J,K+1) 
      call prem(I+1,J,K-1) 
      call prem(I-1,J,K+1) 
      call prem(I-1,J,K-1) 
 
      return
      end
	
c  function to find a set position in the list
c if found, return its value; if not found, return 0
      
      integer*2 function pfind(i,j,k)

      integer i,j,k,pl,ph,pm,dx,dy,dz
      parameter (NLOC=10000)
      common /WORK  / Xs(NLOC),Ys(NLOC),Zs(NLOC),Ms(NLOC),Nset
      integer*2 Xs,Ys,Zs,Ms,Nset

      pl=0
      ph=Nset+1
      do while (ph-pl.gt.1)
        pm=(ph+pl)/2
        dx=i-Xs(pm)
        dy=j-Ys(pm)
        dz=k-Zs(pm)
        if (dx.eq.0 .and. dy.eq.0 .and. dz.eq.0) then
          pfind=Ms(pm)
          return
        endif
        if (dx.gt.0 .or. (dx.eq.0 .and. dy.gt.0) .or.
     &      (dx.eq.0 .and. dy.eq.0 .and. dz.gt.0)) then
          pl=pm
        else
          ph=pm
        endif
      enddo
      
      pfind=0
      return
      end

c subroutine to add a position to the set position list
c if not already present
      
      subroutine padd(i,j,k,m)

      integer i,j,k,m,pl,ph,pm,dx,dy,dz,p
      parameter (NLOC=10000)
      common /WORK  / Xs(NLOC),Ys(NLOC),Zs(NLOC),Ms(NLOC),Nset
      integer*2 Xs,Ys,Zs,Ms,Nset,pfind

      pl=0
      ph=Nset+1
      do while (ph-pl.gt.1)
        pm=(ph+pl)/2
        dx=i-Xs(pm)
        dy=j-Ys(pm)
        dz=k-Zs(pm)
        if (dx.eq.0 .and. dy.eq.0 .and. dz.eq.0) return
        if (dx.gt.0 .or. (dx.eq.0 .and. dy.gt.0) .or.
     &      (dx.eq.0 .and. dy.eq.0 .and. dz.gt.0)) then
          pl=pm
        else
          ph=pm
        endif
      enddo

      if (Nset.eq.0) then
        pm=1
        goto 2000
      endif
      
      if (i.gt.Xs(pm) .or. (i.eq.Xs(pm) .and. j.gt.Ys(pm))
     &  .or. (i.eq.Xs(pm) .and. j.eq.Ys(pm) .and. k.gt.Zs(pm)))
     &    then
        pm=pm+1
      endif

      if (pm.le.Nset) then
        do p=Nset,pm,-1
          Xs(p+1)=Xs(p)
          Ys(p+1)=Ys(p)
          Zs(p+1)=Zs(p)
          Ms(p+1)=Ms(p)
        enddo
      endif
      
2000  Xs(pm)=i
      Ys(pm)=j
      Zs(pm)=k
      Ms(pm)=m
      Nset=Nset+1

      return
      end

c subroutine to remove a position from the set position list
c if present

      subroutine prem(i,j,k)

      integer i,j,k,p,q,dx,dy,dz,pl,ph,pm
      parameter (NLOC=10000)
      common /WORK  / Xs(NLOC),Ys(NLOC),Zs(NLOC),Ms(NLOC),Nset
      integer*2 Xs,Ys,Zs,Ms,Nset

      pl=0
      ph=Nset+1
      do while (ph-pl.gt.1)
        pm=(ph+pl)/2
        dx=i-Xs(pm)
        dy=j-Ys(pm)
        dz=k-Zs(pm)
        if (dx.eq.0 .and. dy.eq.0 .and. dz.eq.0) then
          if (pm.lt.Nset) then
            do p=pm+1,Nset
              Xs(p-1)=Xs(p)
              Ys(p-1)=Ys(p)
              Zs(p-1)=Zs(p)
              Ms(p-1)=Ms(p)
            enddo
          endif
          Nset=Nset-1
          return
        endif
          
        if (dx.gt.0 .or. (dx.eq.0 .and. dy.gt.0) .or.
     &      (dx.eq.0 .and. dy.eq.0 .and. dz.gt.0)) then
          pl=pm
        else
          ph=pm
        endif
      enddo

      return
      end
      
