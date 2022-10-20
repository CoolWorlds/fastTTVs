C ==============================================================
C ==============================================================
C ==============================================================
C
C    MANDEL: Lightcurve generator file.
C
C ==============================================================
C === MODULE MANDEL ===

      MODULE mandelmod
      implicit none
      CONTAINS
C ==============================================================
C === SUBROUTINE: OCCULTQUAD ===

      SUBROUTINE occultquad(z0,u1,u2,p0,muo1,mu0,nz)
C  This routine computes the lightcurve for occultation
C  of a quadratically limb-darkened source without microlensing.
C  Please cite Mandel & Agol (2002) if you make use of this routine
C  in your research.  Please report errors or bugs to agol@tapir.caltech.edu
      implicit none
      integer i,nz
      double precision z0(nz),u1,u2,p,muo1(nz),mu0(nz),
     &       lambdad(nz),etad(nz),lambdae(nz),
     &       pi,x1,x2,x3,z,p0,omega,kap0,kap1,q,Kk,Ek,n,
     &       ellpic_bulirsch,tol,
     &       mpp1,pm1,pp1,abs1mp,abspm5p5


C
C Input:
C
C rs   radius of the source (set to unity)
C z0   impact parameter in units of rs
C p    occulting star size in units of rs
C u1   linear    limb-darkening coefficient (gamma_1 in paper)
C u2   quadratic limb-darkening coefficient (gamma_2 in paper)
C
C Output:
C
C muo1 fraction of flux at each z0 for a limb-darkened source
C mu0  fraction of flux at each z0 for a uniform source
C
C Limb darkening has the form:
C  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
C 
C To use this routine
C
C Now, compute pure occultation curve:
      omega=1.d0-u1/3.d0-u2/6.d0
      pi=dacos(-1.d0)
      tol = 1d-14

C to mesh with fitting routines
      p = dabs(p0)

C For equalities (mpp1 != 1.d0-p)
      mpp1 = 1.d0-p
      pm1 = p-1.d0
      pp1 = 1.d0+p
      abs1mp = dabs(1.d0-p)
      abspm5p5 = 0.5d0+dabs(p-0.5d0)

C Loop over each impact parameter:
      do i=1,nz

C substitute z=z0(i) to shorten expressions
        z=z0(i)
        x1=(p-z)*(p-z)
        x2=(p+z)*(p+z)
        x3=p*p-z*z

C tolerance for double precision equalities
C special case integrations
        if(dabs(p-z).lt.tol) z = p
        if(dabs((p-1)-z).lt.tol) z = p-1.d0
        if(dabs((1-p)-z).lt.tol) z = 1.d0-p
        if(z.lt.tol) z = 0.d0        

C Case 1 - the star is unocculted
        if((z.ge.pp1).or.(p.le.0.d0)) then
          lambdad(i)=0.d0
          etad(i)=0.d0
          lambdae(i)=0.d0
          goto 10
        endif

C Case 11 - the  source is completely occulted:
        if(p.ge.1.d0.and.z.le.pm1) then
           lambdad(i)=0.d0
           etad(i)=0.5d0
           lambdae(i)=1.d0
           goto 10
        endif

C Case 2, 7, 8 - ingress/egress (uniform disk only)
        if(z.ge.abs1mp.and.z.lt.pp1) then
           kap1=dacos(max(min((1.d0-p*p+z*z)/2.d0/z,1.d0),-1.d0))
           kap0=dacos(max(min((p*p+z*z-1.d0)/2.d0/p/z,1.d0),-1.d0))
           lambdae(i)=(p*p*kap0+kap1-0.5d0*dsqrt(max(4.d0*z*z-
     &          (1.d0+z*z-p*p)**2,0.d0)))/pi
           etad(i) = 1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*
     &          kap0-(1.d0+5.d0*p*p+z*z)/4.d0*
     &          dsqrt((1.d0-x1)*(x2-1.d0)))
        endif

C Case 5, 6, 7 - the edge of planet lies at origin of star
        if(z.eq.p) then
           if(p.lt.0.5d0) then 
C             Case 5
              q=2.d0*p
              Call ellke(q,Ek,Kk)
              lambdad(i)=1.d0/3.d0+2.d0/9.d0/pi*(4.d0*(2.d0*p*p-1.d0)*
     &             Ek+(1.d0-4.d0*p*p)*Kk)
              etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
              lambdae(i) = p*p
           else if(p.gt.0.5d0) then
C             Case 7
              q=0.5d0/p
              Call ellke(q,Ek,Kk)
              lambdad(i)=1.d0/3.d0+16.d0*p/9.d0/pi*(2.d0*p*p-1.d0)*
     &             Ek-(32.d0*p**4-20.d0*p*p+3.d0)/9.d0/pi/p*Kk
           else
C             Case 6
              lambdad(i)=1.d0/3.d0-4.d0/pi/9.d0
              etad(i)=3.d0/32.d0
           endif
           goto 10
        endif

C Case 2, Case 8 - ingress/egress (with limb darkening)
        if(((z.gt.abspm5p5).and.(z.lt.pp1)).or.((p.gt.0.5d0).
     &       and.(z.gt.abs1mp).and.(z.lt.p))) then

           q=dsqrt((1.d0-x1)/(x2-x1))
           Call ellke(q,Ek,Kk)
           n=1.d0/x1-1.d0

           lambdad(i)=2.d0/9.d0/pi/dsqrt(x2-x1)*
     &         (((1.d0-x2)*(2.d0*x2+x1-3.d0)-3.d0*
     &            x3*(x2-2.d0))*Kk+(x2-x1)*
     &           (z*z+7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*
     &            ellpic_bulirsch(n,q))
           
           goto 10
        endif
        
C     Case 3, 4, 9, 10 - planet completely inside star
        if(p.le.1.d0.and.z.le.mpp1) then
           etad(i) = p*p/2.d0*(p*p+2.d0*z*z)
           lambdae(i) = p*p

C     Case 4 - edge of planet hits edge of star
           if(z.eq.mpp1) then
              lambdad(i) = 2.d0/3.d0/pi*dacos(1.d0-2.d0*p)-
     &             4.d0/9.d0/pi*dsqrt(p*(1.d0-p))*
     &             (3.d0+2.d0*p-8.d0*p*p)

              if(p.gt.0.5d0) lambdad(i) = lambdad(i)-2.d0/3.d0
              goto 10
           endif

C     Case 10 - origin of planet hits origin of star
           if(z.eq.0) then
              lambdad(i) = -2.d0/3.d0*(1.d0-p*p)**1.5d0
              goto 10
           endif

C     Case 3, Case 9 - anywhere in between
           q=dsqrt((x2-x1)/(1.d0-x1))
           Call ellke(q,Ek,Kk)
           n=x2/x1-1.d0
           lambdad(i)=2.d0/9.d0/pi/dsqrt(1.d0-x1)*((1.d0-5.d0*z*z+p*p+
     &          x3*x3)*Kk+(1.d0-x1)*(z*z+7.d0*p*p-4.d0)*
     &          Ek-3.d0*x3/x1*ellpic_bulirsch(n,q))
        endif

 10     continue

C limb darkened flux
        if(p.gt.z) then
           muo1(i)=1.d0-((1.d0-u1-2.d0*u2)*lambdae(i)+(u1+2.d0*u2)*
     &          (lambdad(i)+2.d0/3.d0)+u2*etad(i))/omega
        else 
           muo1(i)=1.d0-((1.d0-u1-2.d0*u2)*lambdae(i)+(u1+2.d0*u2)*
     &          lambdad(i)+u2*etad(i))/omega
        endif

C uniform disk
        mu0(i)=1.d0-lambdae(i)
      enddo

      return
      END SUBROUTINE occultquad
      END MODULE mandelmod

C ==============================================================
C === FUNCTION: ELLPIC_BULIRSCH ===
C
      function ellpic_bulirsch(n,k)
      implicit none
      double precision n,k,kc,p,m0,c,d,e,f,g,pi,ellpic_bulirsch
c Computes the complete elliptical integral of the third kind using
c the algorithm of Bulirsch (1965):
      pi=dacos(-1.d0); kc=dsqrt(1d0-k*k); p=n+1d0
      if (p.lt.0d0) print *, "Negative p"
      m0=1d0; c=1d0; p=dsqrt(p); d=1d0/p; e=kc
 10   continue
      f = c; c = d/p+c; g = e/p; d= 2d0*(f*g+d)
      p = g + p; g = m0; m0 = kc + m0; 
      if (dabs(1d0-kc/g).gt.1.d-8) then 
         kc = 2d0*dsqrt(e); e=kc*m0
         goto 10
      endif
      ellpic_bulirsch = 0.5d0*pi*(c*m0+d)/(m0*(m0+p))
      return
      end

C ==============================================================
C === SUBROUTINE: ELLKE ===
C
      subroutine ellke(k,ek,kk)
      implicit none
      double precision k,m1,a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,
     &       ee1,ee2,ek1,ek2,logm1,ek,kk
C Computes Hasting's polynomial approximation for the complete 
C elliptic integral of the first (ek) and second (kk) kind 

      m1=1.d0-k*k
      logm1=dlog(m1)

      a1=0.44325141463d0
      a2=0.06260601220d0
      a3=0.04757383546d0
      a4=0.01736506451d0
      b1=0.24998368310d0
      b2=0.09200180037d0
      b3=0.04069697526d0
      b4=0.00526449639d0
      ee1=1.d0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*(-logm1)
      ek=ee1+ee2

      a0=1.38629436112d0
      a1=0.09666344259d0
      a2=0.03590092383d0
      a3=0.03742563713d0
      a4=0.01451196212d0
      b0=0.5d0
      b1=0.12498593597d0
      b2=0.06880248576d0
      b3=0.03328355346d0
      b4=0.00441787012d0
      ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*logm1
      kk=ek1-ek2 

      return
      end

C ==============================================================
C === SUBROUTINE: OCCULTUNIFORM ===
C
      SUBROUTINE occultuniform(b0,w,muo1,nb)
      implicit none
      INTEGER :: i,nb
      REAL(8) :: muo1(nb),w,b0(nb),z,lambdae,kap0,kap1
      REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
      IF(ABS(w-0.5D0).LT.1.0D-3) w=0.5D0
C  This routine computes the lightcurve for occultation
C  of a uniform source without microlensing  (Mandel & Agol 2002).
C Input:
C 
C  rs   radius of the source (set to unity)
C  b0   impact parameter in units of rs
C  w    occulting star size in units of rs
C 
C Output:
C  muo1 fraction of flux at each b0 for a uniform source
C 
C  Now, compute pure occultation curve:
      DO i=1,nb
C  substitute z=b0(i) to shorten expressions
        z=b0(i)
C  the source is unocculted:
C  Table 3, I.
        IF(z.GE.1.d0+w) THEN
          muo1(i)=1.0D0
          GOTO 1
        ENDIF
C  the  source is completely occulted:
C  Table 3, II.
        IF(w.GE.1.d0.AND.z.LE.w-1.0D0) THEN
          muo1(i)=0.0D0
          GOTO 1
        ENDIF
C  the source is partly occulted and the occulting object crosses the limb:
C  Equation (26):
        IF(z.GE.ABS(1.D0-w).AND.z.LE.1.d0+w) THEN
          kap1=ACOS(MIN((1.D0-w*w+z*z)/2.D0/z,1.D0))
          kap0=ACOS(MIN((w*w+z*z-1.D0)/2.D0/w/z,1.D0))
          lambdae=w*w*kap0+kap1
          lambdae=(lambdae-0.5D0*SQRT(MAX(4.D0*z*z-(1.D0+z*z-w*w)**2,
     &            0.D0)))/pi
          muo1(i)=1.0D0-lambdae
        ENDIF
C  the occulting object transits the source star (but doesn't
C  completely cover it):
        IF(z.LE.1.d0-w) muo1(i)=1.0D0-w*w
 1      CONTINUE
      ENDDO
C muo1=1.d0-lambdae
      return
      END

C ==============================================================
