

MODULE module_mp_wsm5

   USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm
   USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep
   USE module_mp_radar

   REAL, PARAMETER, PRIVATE :: dtcldcr     = 120. 
   REAL, PARAMETER, PRIVATE :: n0r = 8.e6         
   REAL, PARAMETER, PRIVATE :: avtr = 841.9       
   REAL, PARAMETER, PRIVATE :: bvtr = 0.8         
   REAL, PARAMETER, PRIVATE :: r0 = .8e-5         
   REAL, PARAMETER, PRIVATE :: peaut = .55        
   REAL, PARAMETER, PRIVATE :: xncr = 3.e8        
   REAL, PARAMETER, PRIVATE :: xmyu = 1.718e-5    
   REAL, PARAMETER, PRIVATE :: avts = 11.72       
   REAL, PARAMETER, PRIVATE :: bvts = .41         
   REAL, PARAMETER, PRIVATE :: n0smax =  1.e11    
   REAL, PARAMETER, PRIVATE :: lamdarmax = 8.e4   
   REAL, PARAMETER, PRIVATE :: lamdasmax = 1.e5   
   REAL, PARAMETER, PRIVATE :: lamdagmax = 6.e4   
   REAL, PARAMETER, PRIVATE :: dicon = 11.9       
   REAL, PARAMETER, PRIVATE :: dimax = 500.e-6    
   REAL, PARAMETER, PRIVATE :: n0s = 2.e6         
   REAL, PARAMETER, PRIVATE :: alpha = .12        
   REAL, PARAMETER, PRIVATE :: pfrz1 = 100.       
   REAL, PARAMETER, PRIVATE :: pfrz2 = 0.66       
   REAL, PARAMETER, PRIVATE :: qcrmin = 1.e-9     
   REAL, PARAMETER, PRIVATE :: eacrc = 1.0        
   REAL, SAVE ::                                      &
             qc0, qck1, pidnc,                        &
             bvtr1,bvtr2,bvtr3,bvtr4,g1pbr,           &
             g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,    &
             precr1,precr2,xmmax,roqimax,bvts1,       &
             bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,     &
             g5pbso2,pvts,pacrs,precs1,precs2,pidn0r, &
             pidn0s,xlv1,pacrc,pi,                    &
             rslopermax,rslopesmax,rslopegmax,        &
             rsloperbmax,rslopesbmax,rslopegbmax,     &
             rsloper2max,rslopes2max,rslopeg2max,     &
             rsloper3max,rslopes3max,rslopeg3max



CONTAINS


  SUBROUTINE wsm5(th, q, qc, qr, qi, qs                            &
                 ,den, pii, p, delz                                &
                 ,delt,g, cpd, cpv, rd, rv, t0c                    &
                 ,ep1, ep2, qmin                                   &
                 ,XLS, XLV0, XLF0, den0, denr                      &
                 ,cliq,cice,psat                                   &
                 ,rain, rainncv                                    &
                 ,snow, snowncv                                    &
                 ,sr                                               &
                 ,refl_10cm, diagflag, do_radar_ref                &
                 ,has_reqc, has_reqi, has_reqs                     &  
                 ,re_cloud, re_ice,   re_snow                      &  
                 ,ids,ide, jds,jde, kds,kde                        &
                 ,ims,ime, jms,jme, kms,kme                        &
                 ,its,ite, jts,jte, kts,kte                        &
                                                                   )

  IMPLICIT NONE


  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                          &
                                                             th,  &
                                                              q,  &
                                                              qc, &
                                                              qi, &
                                                              qr, &
                                                              qs
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                             pii, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                              rd, &
                                                              rv, &
                                                             t0c, &
                                                            den0, &
                                                             cpd, &
                                                             cpv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime , jms:jme ),                           &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr

  INTEGER, INTENT(IN)::                                           &
                                                        has_reqc, &
                                                        has_reqi, &
                                                        has_reqs
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),                     &
        INTENT(INOUT)::                                           &
                                                        re_cloud, &
                                                          re_ice, &
                                                         re_snow

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::     &  
                                                       refl_10cm


  REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL,                 &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv

  INTEGER                  ::   ids2,ide2, jds2,jde2, kds2,kde2
  REAL, DIMENSION( 64 , kte ) ::   t_,q_,p_,delz_,den_
  REAL, DIMENSION( 64 , kte, 2 ) ::   qci_, qrs_
  REAL, DIMENSION( 64 ) :: rain_,rainncv_,sr_,snow_,snowncv_
  INTEGER ::               i,j,k
  INTEGER ::               ii,ic,ip
  INTEGER :: iii
  INTEGER :: ntds
  INTEGER, EXTERNAL :: omp_get_max_threads


      REAL, DIMENSION(kts:kte):: qv1d, t1d, p1d, qr1d, qs1d, dBZ
      LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
      INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref



  REAL, DIMENSION( kts:kte ) :: den1d
  REAL, DIMENSION( kts:kte ) :: qc1d
  REAL, DIMENSION( kts:kte ) :: qi1d
  REAL, DIMENSION( kts:kte ) :: re_qc, re_qi, re_qs



  
  
  
  
!$OMP PARALLEL DO   &
!$OMP PRIVATE ( ii,k,ic,ip,i,j,t_,q_,p_,delz_,qci_,qrs_,den_,rain_,rainncv_,sr_,snow_,snowncv_ ) &
!$OMP SCHEDULE(dynamic)
      DO ip = 1,((1+(ite-its+1)/64)*64)*(jte-jts+1),64
       
       j  = jts+(ip-1)/((1+(ite-its+1)/64)*64)
       IF ( j .ge. jts .and. j .le. jte ) THEN
        ii = its+mod((ip-1),((1+(ite-its+1)/64)*64))
         DO k=kts,kte
          DO ic=1,min(64,ite-ii+1)
            i = ii+ic -1
            t_(ic,k)=th(i,k,j)*pii(i,k,j)
            q_(ic,k)=q(i,k,j)
            p_(ic,k)=p(i,k,j)
            delz_(ic,k)=delz(i,k,j)
            qci_(ic,k,1) = qc(i,k,j)
            qci_(ic,k,2) = qi(i,k,j)
            qrs_(ic,k,1) = qr(i,k,j)
            qrs_(ic,k,2) = qs(i,k,j)
            den_(ic,k) = den(i,k,j)
          ENDDO
         ENDDO
         DO ic=1,min(64,ite-ii+1)
            i = ii+ic -1
            rain_(ic) = rain(i,j)
            rainncv_(ic) = rainncv(i,j)
            sr_(ic) = sr(i,j)
            snow_(ic) = snow(i,j)
            snowncv_(ic) = snowncv(i,j)
         ENDDO
         IF ( min(64,ite-ii+1) .gt. 0 ) THEN
           CALL wsm52D(T=t_, Q=q_, QCI=qci_, QRS=qrs_                &
                      ,DEN=den_                                      &
                      ,P=p_, DELZ=delz_                              &
                      ,DELT=delt,G=g, CPD=cpd, CPV=cpv, RD=rd        &
                      ,RV=rv, T0C=t0c                                &
                      ,EP1=ep1, EP2=ep2, QMIN=qmin                   &
                      ,XLS=XLS, XLV0=XLV0, XLF0=XLF0                 &
                      ,DEN0=den0, DENR=denr                          &
                      ,CLIQ=cliq,CICE=cice,PSAT=psat                 &
                      ,LON=ii,LAT=j                                  &
                      ,RAIN=rain_  ,RAINNCV=rainncv_                 &
                      ,SR=sr_                                        &
                      ,SNOW=snow_      ,SNOWNCV=snowncv_             &
                      ,NX0=64, NK0=kte-kts+1                      &
                      ,IRESTRICT=min(64,ite-ii+1)                 &
                      ,DOIT=.TRUE., KTS=kts,KTE=kte                  )
         ENDIF
         DO K=kts,kte
          DO ic=1,min(64,ite-ii+1)
            i = ii+ic -1
            th(i,k,j)=t_(ic,k)/pii(i,k,j)
            q(i,k,j) = q_(ic,k)
            qc(i,k,j) = qci_(ic,k,1)
            qi(i,k,j) = qci_(ic,k,2)
            qr(i,k,j) = qrs_(ic,k,1)
            qs(i,k,j) = qrs_(ic,k,2)
          ENDDO
         ENDDO
         DO ic=1,min(64,ite-ii+1)
            i = ii+ic -1
            rain(i,j) = rain_(ic)
            rainncv(i,j) = rainncv_(ic)
            sr(i,j) = sr_(ic)
            snow(i,j) = snow_(ic)
            snowncv(i,j) = snowncv_(ic)
         ENDDO
       ENDIF
      ENDDO
      !$OMP END PARALLEL DO

   IF ( PRESENT (diagflag) ) THEN
     IF (diagflag .and. do_radar_ref == 1) then
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( i,j,k,t1d,p1d,qv1d,qr1d,qs1d,dBZ )
      DO j=jts,jte
        DO I=its,ite
          DO K=kts,kte
            t1d(k)=th(i,k,j)*pii(i,k,j)
            p1d(k)=p(i,k,j)
            qv1d(k)=q(i,k,j)
            qr1d(k)=qr(i,k,j)
            qs1d(k)=qs(i,k,j)
          ENDDO
          call refl10cm_wsm5 (qv1d, qr1d, qs1d,                    &
                              t1d, p1d, dBZ, kts, kte, i, j)
          DO k = kts, kte
                   refl_10cm(i,k,j) = MAX(-35., dBZ(k))
          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
     ENDIF
   ENDIF

  END SUBROUTINE wsm5


  SUBROUTINE wsm52D(t, q                                          & 
                   ,qci, qrs, den, p, delz                        &
                   ,delt,g, cpd, cpv, rd, rv, t0c                 &
                   ,ep1, ep2, qmin                                &
                   ,XLS, XLV0, XLF0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lon,lat                                       &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,snow,snowncv                                  &
                   ,nx0,nk0,irestrict                             &
                   ,doit,kts,kte                                  )


  IMPLICIT NONE
























  INTEGER,      INTENT(IN   )    ::   nx0,nk0,irestrict            &
                                     ,lon,lat,kts,kte


  REAL, DIMENSION( 1:64 , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                              t
  REAL, DIMENSION( 1:64 , kts:kte, 2 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qci, &
                                                             qrs
  REAL, DIMENSION( 1:64 , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                               q
  REAL, DIMENSION( 1:64 , kts:kte ),                           &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                             cpd, &
                                                             cpv, &
                                                             t0c, &
                                                            den0, &
                                                              rd, &
                                                              rv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( 1:64 ),                                     &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
  REAL, DIMENSION( 1:64          ),     OPTIONAL,              &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv

  LOGICAL, INTENT(IN) :: doit    

  REAL, DIMENSION( 1:64 , kts:kte , 2) ::                      &
                                                              rh, &
                                                              qs, &
                                                            falk, &
                                                            fall

  REAL, DIMENSION( 1:64 ,       1 , 2) ::                      &
                                                           work1
  REAL, DIMENSION( 1:64 ,       1    ) ::                      &
                                                           work2

  REAL, DIMENSION( 1:64 ,  2) ::                               &
                                                          rslope_v, &
                                                         rslope2_v, &
                                                         rslope3_v, &
                                                         rslopeb_v



  REAL, DIMENSION( 1:64 , kts:kte ) ::                         &
                                                           falkc, &
                                                           fallc, &
                                                              xl, &
                                                             cpm, &
                                                          denfac, &
                                                             xni, &
                                                          denqrs, &
                                                          denqci, &
                                                          n0sfac, &
                                                          workrs, &
                                                          work1c, &
                                                          work2c
  REAL, DIMENSION( 1:64 ) ::                                   &
                                                          delqrs
  REAL, DIMENSION( 1:64 , 1 ) ::                         &
                                                           pigen, &
                                                           pidep, &
                                                           psdep, &
                                                           praut, &
                                                           psaut, &
                                                           prevp, &
                                                           psevp, &
                                                           pracw, &
                                                           psacw, &
                                                           psaci, &
                                                           pcond, &
                                                           psmlt

  INTEGER, DIMENSION( 1:64 ) ::                                &
                                                           mstep, &
                                                           numdt
  REAL, DIMENSION(1:64) ::                             tstepsnow
  REAL, DIMENSION(1:64) ::                             rmstep
  REAL dtcldden, rdelz, rdtcld
  LOGICAL, DIMENSION( 1:64 ) ::                        flgcld
  REAL, DIMENSION(1:64) :: xal, xbl
  REAL  ::                                                        &
            cpmcal, xlcal, diffus,                                &
            viscos, xka, venfac, conden, diffac,                  &
            x, y, z, a, b, c, d, e,                               &
            qdt, holdrr, holdrs, supcol, supcolt, pvt,            &
            coeres, supsat, dtcld, xmi, eacrs, satdt,             &
            vt2i,vt2s,acrfac,                                     &
            qimax, diameter, xni0, roqi0,                         &
            xlwork2, factor, source,        &
            value, xlf, pfrzdtc, pfrzdtr, supice,  holdc, holdci

  REAL, DIMENSION(64) :: fallsum, fallsum_qsi  
  REAL, DIMENSION(64) :: supsat_v,satdt_v,coeres_v


  REAL, DIMENSION( 1:64 )           ::                    tvec1
  REAL ::                                                    temp 
  INTEGER :: i, j, k, mstepmax,                                   &
            iprt, latd, lond, loop, loops, ifsat, n, idim, kdim

  REAL  :: dldti, xb, xai, xbi, xa, hvap, cvap, hsub, dldt, ttp
  REAL :: tr_v(1:64),logtr_v(1:64),supcol_v(1:64),supcolt_v(1:64),xlf_v(1:64),temp_v(1:64)
  REAL :: diameter_v(64),supice_v(64)
  INTEGER :: ifsat_v(64)

  LOGICAL*4 :: lmask(64)








      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)












!DIR$ ASSUME_ALIGNED t:64
!DIR$ ASSUME_ALIGNED qci:64
!DIR$ ASSUME_ALIGNED qrs:64
!DIR$ ASSUME_ALIGNED q:64
!DIR$ ASSUME_ALIGNED den:64
!DIR$ ASSUME_ALIGNED p:64
!DIR$ ASSUME_ALIGNED delz:64
!DIR$ ASSUME_ALIGNED rain:64
!DIR$ ASSUME_ALIGNED rainncv:64
!DIR$ ASSUME_ALIGNED sr:64
!DIR$ ASSUME_ALIGNED snow:64
!DIR$ ASSUME_ALIGNED snowncv:64

      if ( irestrict .le. 0 .OR. .NOT. doit ) return
      idim = 64-1+1
      kdim = kte-kts+1
      lmask = .FALSE.
      do i = 1, min(irestrict,64)
        lmask(i) = .TRUE.
      enddo




      do k = kts, kte
        WHERE(lmask)
          qci(:,k,1) = max(qci(:,k,1),0.0)
          qrs(:,k,1) = max(qrs(:,k,1),0.0)
          qci(:,k,2) = max(qci(:,k,2),0.0)
          qrs(:,k,2) = max(qrs(:,k,2),0.0)
        ENDWHERE
      enddo






      do k = kts, kte
        WHERE(lmask)
          cpm(:,k) = (cpd*(1.-max(q(:,k),qmin))+max(q(:,k),qmin)*cpv)
          xl(:,k) = (xlv0-xlv1*((t(:,k))-t0c))
        ENDWHERE
      enddo




      WHERE(lmask)
        rainncv(:) = 0.
        sr(:) = 0.
        tstepsnow(:) = 0.
      ENDWHERE
      IF(PRESENT (snowncv) .AND. PRESENT (snow)) THEN
        WHERE( lmask ) snowncv(:) = 0.
      ENDIF




      loops = max(nint(delt/dtcldcr),1)
      dtcld = delt/loops
      if(delt.le.dtcldcr) dtcld = delt

      do loop = 1,loops




      WHERE (lmask)
        mstep(:) = 1
        flgcld(:) = .true.
      ENDWHERE


      do k = kts, kte
        CALL vsrec( tvec1(1), den(1,k), 64-1+1)
        WHERE( lmask )
          tvec1(:) = tvec1(:)*den0
        ENDWHERE
        CALL vssqrt( denfac(1,k), tvec1(1), 64-1+1)
      enddo





      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)

      do k = kts, kte
        WHERE( lmask )
          WHERE(t(:,k).lt.ttp) 
            xal(:) = xai
            xbl(:) = xbi
          ELSEWHERE
            xal(:) = xa
            xbl(:) = xb
          ENDWHERE
          tr_v=ttp/t(:,k)
          logtr_v=log(tr_v)
          qs(:,k,1)=psat*exp(logtr_v*(xa)+xb*(1.-tr_v))
          qs(:,k,1) = min(qs(:,k,1),0.99*p(:,k))
          qs(:,k,1) = ep2 * qs(:,k,1) / (p(:,k) - qs(:,k,1))
          qs(:,k,1) = max(qs(:,k,1),qmin)
          rh(:,k,1) = max(q(:,k) / qs(:,k,1),qmin)
          qs(:,k,2)=psat*exp(logtr_v*(xal(:))+xbl(:)*(1.-tr_v))
          qs(:,k,2) = min(qs(:,k,2),0.99*p(:,k))
          qs(:,k,2) = ep2 * qs(:,k,2) / (p(:,k) - qs(:,k,2))
          qs(:,k,2) = max(qs(:,k,2),qmin)
          rh(:,k,2) = max(q(:,k) / qs(:,k,2),qmin)
        ENDWHERE
      enddo






      do k = kts, kte
        do i = 1, min(irestrict,64)
          falk(i,k,1) = 0.
          falk(i,k,2) = 0.
          fall(i,k,1) = 0.
          fall(i,k,2) = 0.
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          xni(i,k) = 1.e3
        enddo
      enddo



      
      do k = kts, kte
        WHERE( lmask )
          tr_v = (den(:,k)*max(qci(:,k,2),qmin))
          tr_v = sqrt(sqrt(tr_v*tr_v*tr_v))
          xni(:,k) = min(max(5.38e7*tr_v,1.e3),1.e6)
        ENDWHERE
      enddo






      do k = kts, kte
        WHERE(lmask)



          WHERE (qrs(:,k,1).le.0.0)
            workrs(:,k) = 0.0
          ELSEWHERE (qrs(:,k,1).le.qcrmin)
            workrs(:,k) = pvtr*rsloperbmax*denfac(:,k)
          ELSEWHERE
          
            workrs(:,k) = pvtr*( exp(log( 1./sqrt(sqrt(pidn0r/((qrs(:,k,1))*(den(:,k))))) )*(bvtr)) )*denfac(:,k)
          ENDWHERE
        ENDWHERE
      enddo

      qrs(:,:,1) = den*qrs(:,:,1)
      call nislfv_rain_plm(idim,kdim,den,denfac,t,delz,workrs,qrs(:,:,1),  &
                           delqrs,dtcld,1,1,irestrict,lon,lat,.true.,1)
      fall(:,:,1) = qrs(:,:,1)*workrs/delz
      qrs(:,:,1) = max(qrs(:,:,1)/den,0.)
      fall(:,1,1) = delqrs/delz(:,1)/dtcld

      do k = kts, kte
        WHERE(lmask)
          WHERE (qrs(:,k,2).le.0.0)
            workrs(:,k) = 0.0
          ELSEWHERE (qrs(:,k,2).le.qcrmin)
            workrs(:,k) = pvts*rslopesbmax*denfac(:,k)
          ELSEWHERE
            workrs(:,k) = pvts*(exp(log(1./                                       &
      sqrt(sqrt(pidn0s*( max(min(exp(alpha*(t0c-t(:,k))),n0smax/n0s),1.))/((qrs(:,k,2))*(den(:,k)))))   &
              ) *(bvts)))*denfac(:,k)
          ENDWHERE
        ENDWHERE
      enddo
      qrs(:,:,2) = den*qrs(:,:,2)
      call nislfv_rain_plm(idim,kdim,den,denfac,t,delz,workrs,qrs(:,:,2),  &
                           delqrs,dtcld,2,1,irestrict,lon,lat,.false.,2)
      fall(:,:,2) = qrs(:,:,2)*workrs/delz
      qrs(:,:,2) = max(qrs(:,:,2)/den,0.)
      fall(:,1,2) = delqrs/delz(:,1)/dtcld

      
      xlf = xlf0
      do k = kte, kts, -1
        psmlt = 0.
        WHERE(lmask)
          supcol_v = t0c-t(:,k)
          n0sfac(:,k) = max(min(exp(alpha*supcol_v),n0smax/n0s),1.)
          WHERE(qrs(:,k,2).le.qcrmin)
            rslope_v(:,2) = rslopesmax
            rslopeb_v(:,2) = rslopesbmax
            rslope2_v(:,2) = rslopes2max
          ELSEWHERE
            rslope_v(:,2) = 1./sqrt(sqrt(pidn0s*(max(min(exp(alpha*(t0c-t(:,k))),n0smax/n0s),1.))/((qrs(:,k,2))*(den(:,k)))))
            rslopeb_v(:,2) = exp(log(rslope_v(:,2))*(bvts))
            rslope2_v(:,2) = rslope_v(:,2)*rslope_v(:,2)
          ENDWHERE
          WHERE (t(:,k).gt.t0c.and.qrs(:,k,2).gt.0.)




            work2(:,1)= (exp(log(((1.496e-6*((t(:,k))*sqrt(t(:,k)))        &
                        /((t(:,k))+120.)/(den(:,k)))/(8.794e-5             &
                        *exp(log(t(:,k))*(1.81))/p(:,k))))                 &
                        *((.3333333)))/sqrt((1.496e-6*((t(:,k))            &
                        *sqrt(t(:,k)))/((t(:,k))+120.)/(den(:,k))))        &
                        *sqrt(sqrt(den0/(den(:,k)))))
            tr_v = rslope2_v(:,2)*sqrt(rslope_v(:,2)*rslopeb_v(:,2))
            psmlt(:,1) = (1.414e3*(1.496e-6*((t(:,k))*sqrt(t(:,k)))        &
                        /((t(:,k))+120.)/(den(:,k)) )*(den(:,k)))          &
                        /xlf*(t0c-t(:,k))*pi/2.                            &
                        *n0sfac(:,k)*(precs1*rslope2_v(:,2)+precs2         &
                        *work2(:,1)*tr_v)
            psmlt(:,1) = min(max(psmlt(:,1)*dtcld/mstep(:),                &
                        -qrs(:,k,2)/mstep(:)),0.)
            qrs(:,k,2) = qrs(:,k,2) + psmlt(:,1)
            qrs(:,k,1) = qrs(:,k,1) - psmlt(:,1)
            t(:,k) = t(:,k) + xlf/cpm(:,k)*psmlt(:,1)
          ENDWHERE
        ENDWHERE
      enddo



      work1c = 0.
      WHERE (qci(:,:,2).gt.0.) 
        work1c =                                                            &
               1.49e4*exp(                                                  &
                 log(                                                       &
                   max(min(dicon * sqrt(den*qci(:,:,2)/xni),dimax), 1.e-25) &
                 )*(1.31)                                                   &
               )
      ENDWHERE
      denqci = den*qci(:,:,2)
      call nislfv_rain_plm(idim,kdim,den,denfac,t,delz,work1c,denqci,  &
                           delqrs,dtcld,1,0,irestrict,lon,lat,.false.,3)
      do k = kts, kte
        WHERE(lmask)
          qci(:,k,2) = max(denqci(:,k)/den(:,k),0.)
        ENDWHERE
      enddo
      WHERE(lmask)
        fallc(:,1) = delqrs(:)/delz(:,1)/dtcld
      ENDWHERE




        fallsum = fall(:,1,1)+fall(:,1,2)+fallc(:,1)
        fallsum_qsi = fall(:,1,2)+fallc(:,1)
        WHERE (lmask .and. fallsum.gt.0.)
          rainncv = fallsum*delz(:,1)/denr*dtcld*1000. + rainncv
          rain = fallsum*delz(:,1)/denr*dtcld*1000. + rain
        ENDWHERE
        WHERE (lmask .and. fallsum_qsi.gt.0.)
          tstepsnow   = fallsum_qsi*delz(:,kts)/denr*dtcld*1000.            &
                        +tstepsnow
          snowncv    = fallsum_qsi*delz(:,kts)/denr*dtcld*1000.            & 
                        +snowncv
          snow    = fallsum_qsi*delz(:,kts)/denr*dtcld*1000. + snow
        ENDWHERE
        WHERE (lmask.and.fallsum.gt.0.)sr=tstepsnow/(rainncv+1.e-12)





      do k = kts, kte
        WHERE(lmask)
          supcol_v = t0c-t(:,k)
          xlf_v = xlf0
          WHERE(supcol_v.ge.0.) xlf_v = xls-xl(:,k)
          WHERE(supcol_v.lt.0..and.qci(:,k,2).gt.0.) 
            qci(:,k,1) = qci(:,k,1) + qci(:,k,2)
            t(:,k) = t(:,k) - xlf_v/cpm(:,k)*qci(:,k,2)
            qci(:,k,2) = 0.
          ENDWHERE




          WHERE(supcol_v.gt.40..and.qci(:,k,1).gt.0.) 
            qci(:,k,2) = qci(:,k,2) + qci(:,k,1)
            t(:,k) = t(:,k) + xlf_v/cpm(:,k)*qci(:,k,1)
            qci(:,k,1) = 0.
          ENDWHERE




          
          WHERE(supcol_v.gt.0..and.qci(:,k,1).gt.0) 
            supcolt_v=min(supcol_v,50.)
            tr_v = min(pfrz1*(exp(pfrz2*supcolt_v)-1.)                        &
            *den(:,k)/denr/xncr*qci(:,k,1)*qci(:,k,1)*dtcld,qci(:,k,1))
            qci(:,k,2) = qci(:,k,2) + tr_v
            t(:,k) = t(:,k) + xlf_v/cpm(:,k)*tr_v
            qci(:,k,1) = qci(:,k,1)-tr_v
          ENDWHERE




          
          WHERE(supcol_v.gt.0..and.qrs(:,k,1).gt.0)
            supcolt_v=min(supcol_v,50.)
            WHERE ( qrs(:,k,1).le.qcrmin)
             temp_v = (rslopermax)
            ELSEWHERE
             temp_v = (1./sqrt(sqrt(pidn0r/((qrs(:,k,1))*(den(:,k))))))
            ENDWHERE
            temp_v = temp_v*temp_v*temp_v*temp_v*temp_v*temp_v*temp_v
            tr_v = min(20.*(pi*pi)*pfrz1*n0r*denr/den(:,k)                    &
                  *(exp(pfrz2*supcolt_v)-1.)*temp_v*dtcld,                    &
                  qrs(:,k,1))
            qrs(:,k,2) = qrs(:,k,2) + tr_v
            t(:,k) = t(:,k) + xlf_v/cpm(:,k)*tr_v
            qrs(:,k,1) = qrs(:,k,1)-tr_v
          ENDWHERE
        ENDWHERE
      enddo










      rdtcld = 1./dtcld

      do k = kts, kte
        do i = 1, min(irestrict,64)
          prevp(i,1) = 0.
          psdep(i,1) = 0.
          praut(i,1) = 0.
          psaut(i,1) = 0.
          pracw(i,1) = 0.
          psaci(i,1) = 0.
          psacw(i,1) = 0.
          pigen(i,1) = 0.
          pidep(i,1) = 0.
          pcond(i,1) = 0.
          psevp(i,1) = 0.
        enddo
        WHERE(lmask)
          work1(:,1,1) = (((((den(:,k))*(xl(:,k))*(xl(:,k)))*((t(:,k))+120.)        &      
                          *(den(:,k)))/(1.414e3*(1.496e-6*((t(:,k))*sqrt(t(:,k))))    &
                          *(den(:,k))*(rv*(t(:,k))*(t(:,k)))))                        &
                          +p(:,k)/((qs(:,k,1))*(8.794e-5*exp(log(t(:,k))*(1.81)))))
          work1(:,1,2) = ((((den(:,k))*(xls)*(xls))*((t(:,k))+120.)*(den(:,k)))&
                       /(1.414e3*(1.496e-6*((t(:,k))*sqrt(t(:,k))))*(den(:,k)) &
                       *(rv*(t(:,k))*(t(:,k))))                                &
                      + p(:,k)/(qs(:,k,2)*(8.794e-5*exp(log(t(:,k))*(1.81)))))
          work2(:,1) = (exp(.3333333*log(((1.496e-6 * ((t(:,k))*sqrt(t(:,k))))      &
                          *p(:,k))/(((t(:,k))+120.)*den(:,k)*(8.794e-5                &
                          *exp(log(t(:,k))*(1.81))))))*sqrt(sqrt(den0/(den(:,k)))))   &
                          /sqrt((1.496e-6*((t(:,k))*sqrt(t(:,k))))                    &
                          /(((t(:,k))+120.)*den(:,k)))
        ENDWHERE










        WHERE(lmask)
          WHERE(qrs(:,k,1).le.qcrmin)
            rslope_v(:,1) = rslopermax
            rslopeb_v(:,1) = rsloperbmax
            rslope2_v(:,1) = rsloper2max
            rslope3_v(:,1) = rsloper3max
          ELSEWHERE
            rslope_v(:,1) = 1./sqrt(sqrt(pidn0r/((qrs(:,k,1))*(den(:,k)))))
            rslopeb_v(:,1) = exp(log(rslope_v(:,1))*(bvtr))
            rslope2_v(:,1) = rslope_v(:,1)*rslope_v(:,1)
            rslope3_v(:,1) = rslope2_v(:,1)*rslope_v(:,1)
          ENDWHERE
          supsat_v = max(q(:,k),qmin)-qs(:,k,1)
          satdt_v = supsat_v/dtcld




          WHERE (qci(:,k,1).gt.qc0) 
            praut(:,1) = qck1*exp(log(qci(:,k,1))*((7./3.)))
            praut(:,1) = min(praut(:,1),qci(:,k,1)/dtcld)
          ENDWHERE




    
          WHERE (qrs(:,k,1).gt.qcrmin.and.qci(:,k,1).gt.qmin)
            pracw(:,1) = min(pacrr*rslope3_v(:,1)*rslopeb_v(:,1)               &
                         *qci(:,k,1)*denfac(:,k),qci(:,k,1)/dtcld)
          ENDWHERE




          WHERE(qrs(:,k,1).gt.0.)
            coeres_v = rslope2_v(:,1)*sqrt(rslope_v(:,1)*rslopeb_v(:,1))
            prevp(:,1) = (rh(:,k,1)-1.)*(precr1*rslope2_v(:,1)                 &
                         +precr2*work2(:,1)*coeres_v)/work1(:,1,1)
            WHERE(prevp(:,1).lt.0.)
              prevp(:,1) = max(prevp(:,1),-qrs(:,k,1)/dtcld)
              prevp(:,1) = max(prevp(:,1),satdt_v/2)
            ELSEWHERE
              prevp(:,1) = min(prevp(:,1),satdt_v/2)
            ENDWHERE
          ENDWHERE
        ENDWHERE

















        WHERE(lmask)



          WHERE(qrs(:,k,2).le.qcrmin)
            rslope_v(:,2) = rslopesmax
            rslopeb_v(:,2) = rslopesbmax
            rslope2_v(:,2) = rslopes2max
            rslope3_v(:,2) = rslopes3max
          ELSEWHERE
            rslope_v(:,2) = 1./sqrt(sqrt(pidn0s*(max(min(exp(alpha*(t0c-t(:,k))),n0smax/n0s),1.))/((qrs(:,k,2))*(den(:,k)))))
            rslopeb_v(:,2) = exp(log(rslope_v(:,2))*(bvts))
            rslope2_v(:,2) = rslope_v(:,2)*rslope_v(:,2)
            rslope3_v(:,2) = rslope2_v(:,2)*rslope_v(:,2)
          ENDWHERE
        ENDWHERE

        WHERE(lmask)
          supcol_v = t0c-t(:,k)
          n0sfac(:,1) = max(min(exp(alpha*supcol_v),n0smax/n0s),1.)
          supsat_v = max(q(:,k),qmin)-qs(:,k,2)
          satdt_v = supsat_v/dtcld
          ifsat_v(:) = 0



          temp_v = (den(:,k)*max(qci(:,k,2),qmin))
          temp_v = sqrt(sqrt(temp_v*temp_v*temp_v))
          xni(:,1) = min(max(5.38e7*temp_v,1.e3),1.e6)

          psaci(:,1) = 0.
          psacw(:,1) = 0.
          WHERE( supcol_v .gt. 0 )
            WHERE(qrs(:,k,2).gt.qcrmin.and.qci(:,k,2).gt.qmin)
              diameter_v  = min(dicon * sqrt(den(:,k)*qci(:,k,2)/xni(:,1)),dimax)




              psaci(:,1) = pi*qci(:,k,2)*(exp(0.07*(-supcol_v)))*n0s*n0sfac(:,1)                 &
                           *abs((pvts*rslopeb_v(:,2)*denfac(:,k))              &
                           -(1.49e4*diameter_v**1.31))*(2.*rslope3_v(:,2)+2.   &
                           *diameter_v*rslope2_v(:,2)            &
                           +diameter_v*diameter_v*rslope_v(:,2))/4.
            ENDWHERE
          ENDWHERE




          WHERE(qrs(:,k,2).gt.qcrmin.and.qci(:,k,1).gt.qmin)
            psacw(:,1) = min(pacrc*n0sfac(:,1)*rslope3_v(:,2)                  &
                           *rslopeb_v(:,2)*qci(:,k,1)*denfac(:,k)              &
                           ,qci(:,k,1)*rdtcld)
          ENDWHERE
          pidep(:,1) = 0.
          psdep(:,1) = 0.
          WHERE(supcol_v .gt. 0)




            WHERE(qci(:,k,2).gt.0.and.ifsat_v(:).ne.1)
              diameter_v = dicon * sqrt(den(:,k)*qci(:,k,2)/xni(:,1))
              pidep(:,1) = 4.*diameter_v*xni(:,1)*(rh(:,k,2)-1.)/work1(:,1,2)
              supice_v = satdt_v-prevp(:,1)
              WHERE(pidep(:,1).lt.0.)
                pidep(:,1) = max(max(pidep(:,1),satdt_v*.5),supice_v)
                pidep(:,1) = max(pidep(:,1),-qci(:,k,2)*rdtcld)
              ELSEWHERE
                pidep(:,1) = min(min(pidep(:,1),satdt_v*.5),supice_v)
              ENDWHERE
              WHERE(abs(prevp(:,1)+pidep(:,1)).ge.abs(satdt_v)) ifsat_v(:) = 1
            ENDWHERE




            WHERE(qrs(:,k,2).gt.0..and.ifsat_v(:).ne.1)
              psdep(:,1) = (rh(:,k,2)-1.)*n0sfac(:,1)                          &
                           *(precs1*rslope2_v(:,2)+precs2                      &
                           *work2(:,1)*(rslope2_v(:,2)*sqrt(rslope_v(:,2)*rslopeb_v(:,2))))/work1(:,1,2)
              supice_v = satdt_v-prevp(:,1)-pidep(:,1)
              WHERE(psdep(:,1).lt.0.)
                psdep(:,1) = max(psdep(:,1),-qrs(:,k,2)*rdtcld)
                psdep(:,1) = max(max(psdep(:,1),satdt_v*.5),supice_v)
              ELSEWHERE
                psdep(:,1) = min(min(psdep(:,1),satdt_v*.5),supice_v)
              ENDWHERE
              WHERE(abs(prevp(:,1)+pidep(:,1)+psdep(:,1)).ge.abs(satdt_v))          &
                ifsat_v(:) = 1
            ENDWHERE




            pigen(:,1) = 0.
            WHERE(supsat_v.gt.0.and.ifsat_v(:).ne.1)
              supice_v = satdt_v-prevp(:,1)-pidep(:,1)-psdep(:,1)
              pigen(:,1) = max(0.,((4.92e-11*exp(log(1.e3*exp(0.1*supcol_v)) &
                          *(1.33)))/den(:,k)-max(qci(:,k,2),0.))          &
                          *rdtcld)
              pigen(:,1) = min(min(pigen(:,1),satdt_v),supice_v)
            ENDWHERE





            psaut(:,1) = 0.
            WHERE(qci(:,k,2).gt.0.)
              psaut(:,1) = max(0.,(qci(:,k,2)-(roqimax/den(:,k)))*rdtcld)
            ENDWHERE
          ENDWHERE




          psevp(:,1) = 0.
          WHERE(supcol_v.lt.0.)
            WHERE(qrs(:,k,2).gt.0..and.rh(:,k,1).lt.1.)
              psevp(:,1) = psdep(:,1)*work1(:,1,2)/work1(:,1,1)
            ENDWHERE
            psevp(:,1) = min(max(psevp(:,1),-qrs(:,k,2)*rdtcld),0.)
          ENDWHERE
        ENDWHERE










        do i = 1, min(irestrict,64)
          if(t(i,k).le.t0c) then



            value = max(qmin,qci(i,k,1))
            source = (praut(i,1)+pracw(i,1)+psacw(i,1))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,1) = praut(i,1)*factor
              pracw(i,1) = pracw(i,1)*factor
              psacw(i,1) = psacw(i,1)*factor
            endif



            value = max(qmin,qci(i,k,2))
            source = (psaut(i,1)+psaci(i,1)-pigen(i,1)-pidep(i,1))*dtcld
            if (source.gt.value) then
              factor = value/source
              psaut(i,1) = psaut(i,1)*factor
              psaci(i,1) = psaci(i,1)*factor
              pigen(i,1) = pigen(i,1)*factor
              pidep(i,1) = pidep(i,1)*factor
            endif




            value = max(qmin,qrs(i,k,1))
            source = (-praut(i,1)-pracw(i,1)-prevp(i,1))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,1) = praut(i,1)*factor
              pracw(i,1) = pracw(i,1)*factor
              prevp(i,1) = prevp(i,1)*factor
            endif



            value = max(qmin,qrs(i,k,2))
            source = (-psdep(i,1)-psaut(i,1)-psaci(i,1)-psacw(i,1))*dtcld  
            if (source.gt.value) then
              factor = value/source
              psdep(i,1) = psdep(i,1)*factor
              psaut(i,1) = psaut(i,1)*factor
              psaci(i,1) = psaci(i,1)*factor
              psacw(i,1) = psacw(i,1)*factor
            endif

            work2(i,1)=-(prevp(i,1)+psdep(i,1)+pigen(i,1)+pidep(i,1))

            q(i,k) = q(i,k)+work2(i,1)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,1)+pracw(i,1)                 &
                        +psacw(i,1))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,1)+pracw(i,1)                 &
                        +prevp(i,1))*dtcld,0.)
            qci(i,k,2) = max(qci(i,k,2)-(psaut(i,1)+psaci(i,1)                 &
                        -pigen(i,1)-pidep(i,1))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psdep(i,1)+psaut(i,1)                 &
                        +psaci(i,1)+psacw(i,1))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xls*(psdep(i,1)+pidep(i,1)+pigen(i,1))                  &
                      -xl(i,k)*prevp(i,1)-xlf*psacw(i,1)
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          else



            value = max(qmin,qci(i,k,1))
            source=(praut(i,1)+pracw(i,1)+psacw(i,1))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,1) = praut(i,1)*factor
              pracw(i,1) = pracw(i,1)*factor
              psacw(i,1) = psacw(i,1)*factor
            endif



            value = max(qmin,qrs(i,k,1))
            source = (-praut(i,1)-pracw(i,1)-prevp(i,1)-psacw(i,1))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,1) = praut(i,1)*factor
              pracw(i,1) = pracw(i,1)*factor
              prevp(i,1) = prevp(i,1)*factor
              psacw(i,1) = psacw(i,1)*factor
            endif  



            value = max(qcrmin,qrs(i,k,2))
            source=(-psevp(i,1))*dtcld
            if (source.gt.value) then
              factor = value/source
              psevp(i,1) = psevp(i,1)*factor
            endif
            work2(i,1)=-(prevp(i,1)+psevp(i,1))

            q(i,k) = q(i,k)+work2(i,1)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,1)+pracw(i,1)                 &
                        +psacw(i,1))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,1)+pracw(i,1)                 &
                        +prevp(i,1) +psacw(i,1))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+psevp(i,1)*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xl(i,k)*(prevp(i,1)+psevp(i,1))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld

          endif
        enddo
      enddo   




      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = kts, kte
        WHERE(lmask)
          tr_v=ttp/t(:,k)
          logtr_v = log(tr_v)
          qs(:,k,1)=psat*exp(logtr_v*(xa)+xb*(1.-tr_v))
          qs(:,k,1) = min(qs(:,k,1),0.99*p(:,k))
          qs(:,k,1) = ep2 * qs(:,k,1) / (p(:,k) - qs(:,k,1))
          qs(:,k,1) = max(qs(:,k,1),qmin)
        ENDWHERE
      enddo






      do k = kts, kte
        do i = 1, min(irestrict,64)
          work1(i,1,1) = ((max(q(i,k),qmin)-(qs(i,k,1)))/(1.+(xl(i,k))         &   
                        *(xl(i,k))/(rv*(cpm(i,k)))*(qs(i,k,1))                 &
                        /((t(i,k))*(t(i,k)))))
          pcond(i,1) = min(max(work1(i,1,1)/dtcld,0.),max(q(i,k),0.)/dtcld)
          if(qci(i,k,1).gt.0..and.work1(i,1,1).lt.0.)                          &
            pcond(i,1) = max(work1(i,1,1),-qci(i,k,1))/dtcld
          q(i,k) = q(i,k)-pcond(i,1)*dtcld
          qci(i,k,1) = max(qci(i,k,1)+pcond(i,1)*dtcld,0.)
          t(i,k) = t(i,k)+pcond(i,1)*xl(i,k)/cpm(i,k)*dtcld
        enddo
      enddo





      do k = kts, kte
        do i = 1, min(irestrict,64)
          if(qci(i,k,1).le.qmin) qci(i,k,1) = 0.0
          if(qci(i,k,2).le.qmin) qci(i,k,2) = 0.0
        enddo
      enddo
      enddo                  

  END SUBROUTINE wsm52d


  SUBROUTINE slope_rain(qrs,den,denfac,t,rslope,rslopeb,                   &
                            vt,irestrict,kts,kte,lmask) 
  IMPLICIT NONE
  INTEGER       ::               irestrict,kts,kte
  REAL, DIMENSION( 1:64 , kts:kte) ::                                       &
                                                                          qrs, &
                                                                           vt, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, DIMENSION( 1:64          ) ::                                       &
                                                                       rslope, &
                                                                      rslopeb
  REAL, PARAMETER  :: t0c = 273.15
  REAL       ::  lamdar, x, y, z
  LOGICAL*4 :: lmask(1:64)
  integer :: i, j, k



      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      

      do k = kts, kte
        WHERE(lmask)
          WHERE(qrs(:,k).le.qcrmin)
            rslope(:) = rslopermax
            rslopeb(:) = rsloperbmax
          ELSEWHERE
            rslope(:) = 1./sqrt(sqrt(pidn0r/((qrs(:,k))*(den(:,k)))))
            rslopeb(:) = rslope(:)**bvtr

          ENDWHERE
          WHERE(qrs(:,k).le.0.0) 
            vt(:,k) = 0.0
          ELSEWHERE
            vt(:,k) = pvtr*rslopeb(:)*denfac(:,k)
          ENDWHERE
        ENDWHERE
      enddo
  END SUBROUTINE slope_rain

  SUBROUTINE slope_snow(qrs,den,denfac,t,rslope,rslopeb,                   &
                        vt,irestrict,kts,kte,lmask)
  IMPLICIT NONE
  INTEGER       ::               irestrict,kts,kte
  REAL, DIMENSION( 1:64 , kts:kte) ::                                       &
                                                                          qrs, &
                                                                           vt, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, DIMENSION( 1:64          ) ::                                       &
                                                                       rslope, &
                                                                      rslopeb
  REAL, PARAMETER  :: t0c = 273.15
  integer :: i, j, k
  LOGICAL*4 lmask(64)




      do k = kts, kte
        WHERE(lmask)



          WHERE(qrs(:,k).le.qcrmin)
            rslope(:) = rslopesmax
            rslopeb(:) = rslopesbmax
          ELSEWHERE
            rslope(:) = 1./sqrt(sqrt(pidn0s*((max(min(exp(alpha*(t0c-t(:,k))),n0smax/n0s),1.)))/((qrs(:,k))*(den(:,k)))))
            rslopeb(:) = rslope(:)**bvts

          ENDWHERE
          WHERE(qrs(:,k).le.0.0)
            vt(:,k) = 0.0
          ELSEWHERE
            vt(:,k) = pvts*rslopeb(:)*denfac(:,k)
          ENDWHERE
        ENDWHERE
      enddo
  END SUBROUTINE slope_snow


  SUBROUTINE nislfv_rain_plm(im0,km,den,denfac,tk,dz,ww0,qq0,precip0,dt,id,iter,irestrict,lon,lat,doit,call)






















      implicit none
      integer,intent(in ) ::  im0,km,id,irestrict
      real, intent(in   ) ::  den(64,km)
      real, intent(in   ) ::  denfac(64,km)
      real, intent(in   ) ::  tk(64,km)
      real, intent(in   ) ::  dz(64,km)
      real, intent(in   ) ::  ww0(64,km)
      real, intent(  out) ::  precip0(64)
      real, intent(inout) ::  qq0(64,km)
      real, intent(in   ) ::  dt
      logical, intent(in) :: doit
      integer :: call

      integer  i,k,m,kk,iter
      integer n
      real  dim(64),dip(64),c1,con1,fa1,fa2
      real  allold(64), allnew(64), zz(64), dzamin(64), cflmax(64), decfl(64)

      real  qq(64,km), ww(64,km),precip(64)
      real  wd(64,km), wa(64,km), was(64,km)
      real  wi(64,km+1), zi(64,km+1), za(64,km+1)
      real  qn(64,km),qr(64,km),tmp(64,km),tmp1(64,km),tmp2(64,km),tmp3(64,km)
      real  dza(64,km+1), qa(64,km+1), qmi(64,km+1), qpi(64,km+1)
      logical*4 lmask(64)

      INTEGER minkb, minkt
      LOGICAL, DIMENSION(64) :: intp_mask, tmask
      INTEGER, DIMENSION(64) :: kb, kt
      REAL,    DIMENSION(64) :: tl,tl2,th,th2,qqd,qqh,qql,zsum,qsum,dql,dqh
      REAL,    DIMENSION(64) :: za_gath_t,za_gath_b
      REAL,    DIMENSION(64) :: qa_gath_b
      REAL,    DIMENSION(64) :: dza_gath_t,dza_gath_b
      REAL,    DIMENSION(64) :: qpi_gath_t,qpi_gath_b
      REAL,    DIMENSION(64) :: qmi_gath_t,qmi_gath_b

      integer, intent(in) ::  lat,lon

      precip = 0.0


       qq = qq0
       ww = ww0

      lmask = .false.
      do i = 1,min( irestrict, 64 )
        lmask(i) = .true.
      enddo
      allold = 0.0
      do k=1,km
        where(lmask)allold = allold + qq(:,k)
      enddo
      lmask = lmask .and. ( allold .gt. 0.0 ) 
      IF ( .NOT. ANY(lmask) ) THEN
        precip0 = precip
        RETURN
      ENDIF



      zi(:,1)=0.0
      do k=1,km
        where(lmask) zi(:,k+1) = zi(:,k)+dz(:,k)
      enddo


     wd = ww
     DO n = 0, iter
      where(lmask)


        wi(:,1) = ww(:,1)
        wi(:,km+1) = ww(:,km)
      endwhere
      do k=2,km
        where(lmask)wi(:,k) = (ww(:,k)*dz(:,k-1)+ww(:,k-1)*dz(:,k))/(dz(:,k-1)+dz(:,k))
      enddo


      fa1 = 9./16.
      fa2 = 1./16.
      where(lmask)
        wi(:,1) = ww(:,1)
        wi(:,2) = 0.5*(ww(:,2)+ww(:,1))
      endwhere
      do k=3,km-1
        where(lmask)wi(:,k) = fa1*(ww(:,k)+ww(:,k-1))-fa2*(ww(:,k+1)+ww(:,k-2))
      enddo
      where(lmask)
        wi(:,km) = 0.5*(ww(:,km)+ww(:,km-1))
        wi(:,km+1) = ww(:,km)
      endwhere


      do k=2,km
        where(lmask .and. ww(:,k).eq.0.0 ) wi(:,k)=ww(:,k-1)
      enddo


      con1 = 0.05
      do k=km,1,-1
        where (lmask) 
          decfl = (wi(:,k+1)-wi(:,k))*dt/dz(:,k)
        elsewhere
          decfl = 0.
        endwhere
        where (lmask ) 
          where (decfl .gt. con1 )
            wi(:,k) = wi(:,k+1) - con1*dz(:,k)/dt
          endwhere
        endwhere
      enddo


      do k=1,km+1
        where (lmask) za(:,k) = zi(:,k) - wi(:,k)*dt
      enddo

      do k=1,km
        where (lmask) dza(:,k) = za(:,k+1)-za(:,k)
      enddo
      where (lmask) dza(:,km+1) = zi(:,km+1) - za(:,km+1)



      do k=1,km
        where (lmask) qa(:,k) = qq(:,k)*dz(:,k)/dza(:,k)
        where (lmask) qr(:,k) = qa(:,k)/den(:,k)
      enddo
      where (lmask) qa(:,km+1) = 0.0


      if( n.le.iter-1 ) then
        if (id.eq.1) then
          call slope_rain(qr,den,denfac,tk,tmp,tmp1,wa,irestrict,1,km,lmask)
        else
          call slope_snow(qr,den,denfac,tk,tmp,tmp1,wa,irestrict,1,km,lmask)
        endif 
        do k=1,km
          if( n.ge.1 ) where (lmask) wa(:,k)=0.5*(wa(:,k)+was(:,k))
          where (lmask ) ww(:,k) = 0.5* ( wd(:,k)+wa(:,k) )
        enddo
        was = wa
      endif
     ENDDO  


      do k=2,km
        where (lmask )
          dip=(qa(:,k+1)-qa(:,k))/(dza(:,k+1)+dza(:,k))
          dim=(qa(:,k)-qa(:,k-1))/(dza(:,k-1)+dza(:,k))
          where( dip*dim.le.0.0 )
            qmi(:,k)=qa(:,k)
            qpi(:,k)=qa(:,k)
          elsewhere
            qpi(:,k)=qa(:,k)+0.5*(dip+dim)*dza(:,k)
            qmi(:,k)=2.0*qa(:,k)-qpi(:,k)
            where( qpi(:,k).lt.0.0 .or. qmi(:,k).lt.0.0 )
              qpi(:,k) = qa(:,k)
              qmi(:,k) = qa(:,k)
            endwhere
          endwhere
        endwhere
      enddo
      where (lmask )
        qpi(:,1)=qa(:,1)
        qmi(:,1)=qa(:,1)
        qmi(:,km+1)=qa(:,km+1)
        qpi(:,km+1)=qa(:,km+1)
      endwhere


      qn = 0.0
      kb=1  
      kt=1  
      INTP : do k=1,km
             kb=max(kb-1,1)
             kt=max(kt-1,1)

             intp_mask = ( zi(:,k).lt.za(:,km+1) .AND. lmask )
             tmask = intp_mask
             minkb = 999
             minkt = 999
             DO i=1,64
               IF ( tmask(i) .AND. kb(i) .lt. minkb ) minkb = kb(i)
               IF ( tmask(i) .AND. kt(i) .lt. minkt ) minkt = kt(i)
             ENDDO
             find_kb : do kk=minkb,km
               WHERE ( tmask .AND. zi(:,k).le.za(:,kk+1) )
                 kb = kk
                 tmask = .FALSE.
               END WHERE
             enddo find_kb

             tmask = intp_mask
             find_kt : do kk=minkt,km
               WHERE ( tmask .AND. zi(:,k+1).le.za(:,kk) )
                 kt = kk
                 tmask = .FALSE.
               END WHERE
             enddo find_kt
             kt = max(kt - 1,1)


!DEC$ SIMD
             DO i = 1, 64
               qa_gath_b(i) = qa((i+(kb(i)-1)*64),1)
               za_gath_b(i) = za((i+(kb(i)-1)*64),1)
               dza_gath_b(i) = dza((i+(kb(i)-1)*64),1)
               qpi_gath_b(i) = qpi((i+(kb(i)-1)*64),1)
               qmi_gath_b(i) = qmi((i+(kb(i)-1)*64),1)
             ENDDO
!DEC$ SIMD
             DO i = 1, 64
               za_gath_t(i) = za((i+(kt(i)-1)*64),1)
               dza_gath_t(i) = dza((i+(kt(i)-1)*64),1)
               qpi_gath_t(i) = qpi((i+(kt(i)-1)*64),1)
               qmi_gath_t(i) = qmi((i+(kt(i)-1)*64),1)
             ENDDO

             WHERE ( kt .eq. kb .AND. intp_mask ) 
               tl=(zi(:,k)-za_gath_b)/dza_gath_b
               th=(zi(:,k+1)-za_gath_b)/dza_gath_b
               tl2 = tl*tl
               th2 = th*th
               qqd=0.5*(qpi_gath_b-qmi_gath_b)
               qqh=qqd*th2+qmi_gath_b*th
               qql=qqd*tl2+qmi_gath_b*tl
               qn(:,k) = (qqh-qql)/(th-tl)
             ELSE WHERE ( kt .gt. kb .AND. intp_mask ) 
               tl=(zi(:,k)-za_gath_b)/dza_gath_b
               tl2=tl*tl
               qqd=0.5*(qpi_gath_b-qmi_gath_b)
               qql=qqd*tl2+qmi_gath_b*tl
               dql = qa_gath_b-qql
               zsum  = (1.-tl)*dza_gath_b
               qsum  = dql*dza_gath_b
             END WHERE
             DO i = 1, 64
               if( kt(i)-kb(i).gt.1 .AND. intp_mask(i) ) then
                 do m=kb(i)+1,kt(i)-1
                     zsum(i) = zsum(i) + dza(i,m)
                     qsum(i) = qsum(i) + qa(i,m) * dza(i,m)
                 enddo
               endif
             enddo
             WHERE ( kt .gt. kb .AND. intp_mask )
               th=(zi(:,k+1)-za_gath_t)/dza_gath_t
               th2 = th*th
               qqd=0.5*(qpi_gath_t-qmi_gath_t)
               dqh=qqd*th2+qmi_gath_t*th
               zsum  = zsum + th*dza_gath_t
               qsum  = qsum + dqh*dza_gath_t
               qn(:,k) = qsum/zsum
             END WHERE
       ENDDO intp


      intp_mask = lmask
      sum_precip: do k=1,km
             WHERE ( za(:,k).lt.0.0 .and. za(:,k+1).lt.0.0 .AND. intp_mask )
               precip = precip + qa(:,k)*dza(:,k)
             ELSE WHERE ( za(:,k).lt.0.0 .and. za(:,k+1).ge.0.0 .AND. intp_mask)
               precip = precip + qa(:,k)*(0.0-za(:,k))
               intp_mask = .FALSE.
             END WHERE
      enddo sum_precip






      do k = 1,km
        where(lmask) qq0(:,k) = qn(:,k)
      enddo
      precip0 = 0.
      where (lmask) precip0 = precip






  END SUBROUTINE nislfv_rain_plm








      subroutine refl10cm_wsm5 (qv1d, qr1d, qs1d,                       &
                       t1d, p1d, dBZ, kts, kte, ii, jj)

      IMPLICIT NONE


      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                      qv1d, qr1d, qs1d, t1d, p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ


      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho
      REAL, DIMENSION(kts:kte):: rr, rs
      REAL:: temp_C

      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilams
      DOUBLE PRECISION, DIMENSION(kts:kte):: N0_r, N0_s
      DOUBLE PRECISION:: lamr, lams
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs

      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow
      DOUBLE PRECISION:: fmelt_s

      INTEGER:: i, k, k_0, kbot, n
      LOGICAL:: melti

      DOUBLE PRECISION:: cback, x, eta, f_d
      REAL, PARAMETER:: R=287.



      do k = kts, kte
         dBZ(k) = -35.0
      enddo




      do k = kts, kte
         temp(k) = t1d(k)
         temp_C = min(-0.001, temp(K)-273.15)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))

         if (qr1d(k) .gt. 1.E-9) then
            rr(k) = qr1d(k)*rho(k)
            N0_r(k) = n0r
            lamr = (xam_r*xcrg(3)*N0_r(k)/rr(k))**(1./xcre(1))
            ilamr(k) = 1./lamr
            L_qr(k) = .true.
         else
            rr(k) = 1.E-12
            L_qr(k) = .false.
         endif

         if (qs1d(k) .gt. 1.E-9) then
            rs(k) = qs1d(k)*rho(k)
            N0_s(k) = min(n0smax, n0s*exp(-alpha*temp_C))
            lams = (xam_s*xcsg(3)*N0_s(k)/rs(k))**(1./xcse(1))
            ilams(k) = 1./lams
            L_qs(k) = .true.
         else
            rs(k) = 1.E-12
            L_qs(k) = .false.
         endif
      enddo




      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt.273.15) .and. L_qr(k) .and. L_qs(k+1) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue







      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (xam_s/900.0)*(xam_s/900.0)          &
                                 * N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
      enddo










      if (melti .and. k_0.ge.kts+1) then
       do k = k_0-1, kts, -1


          if (L_qs(k) .and. L_qs(k_0) ) then
           fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           lams = 1./ilams(k)
           do n = 1, nrbins
              x = xam_s * xxDs(n)**xbm_s
              call rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif
       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k))*1.d18)
      enddo

      end subroutine refl10cm_wsm5


  SUBROUTINE wsm5init(den0,denr,dens,cl,cpv,allowed_to_read)

  IMPLICIT NONE


   REAL, INTENT(IN) :: den0,denr,dens,cl,cpv
   LOGICAL, INTENT(IN) :: allowed_to_read

   pi = 4.*atan(1.)
   xlv1 = cl-cpv

   qc0  = 4./3.*pi*denr*r0**3*xncr/den0  
   qck1 = .104*9.8*peaut/(xncr*denr)**(1./3.)/xmyu*den0**(4./3.) 
   pidnc = pi*denr/6.        

   bvtr1 = 1.+bvtr
   bvtr2 = 2.5+.5*bvtr
   bvtr3 = 3.+bvtr
   bvtr4 = 4.+bvtr
   g1pbr = rgmma(bvtr1)
   g3pbr = rgmma(bvtr3)
   g4pbr = rgmma(bvtr4)            
   g5pbro2 = rgmma(bvtr2)          
   pvtr = avtr*g4pbr/6.
   eacrr = 1.0
   pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
   precr1 = 2.*pi*n0r*.78
   precr2 = 2.*pi*n0r*.31*avtr**.5*g5pbro2
   xmmax = (dimax/dicon)**2
   roqimax = 2.08e22*dimax**8

   bvts1 = 1.+bvts
   bvts2 = 2.5+.5*bvts
   bvts3 = 3.+bvts
   bvts4 = 4.+bvts
   g1pbs = rgmma(bvts1)    
   g3pbs = rgmma(bvts3)
   g4pbs = rgmma(bvts4)    
   g5pbso2 = rgmma(bvts2)
   pvts = avts*g4pbs/6.
   pacrs = pi*n0s*avts*g3pbs*.25
   precs1 = 4.*n0s*.65
   precs2 = 4.*n0s*.44*avts**.5*g5pbso2
   pidn0r =  pi*denr*n0r
   pidn0s =  pi*dens*n0s
   pacrc = pi*n0s*avts*g3pbs*.25*eacrc

   rslopermax = 1./lamdarmax
   rslopesmax = 1./lamdasmax
   rsloperbmax = rslopermax ** bvtr
   rslopesbmax = rslopesmax ** bvts
   rsloper2max = rslopermax * rslopermax
   rslopes2max = rslopesmax * rslopesmax
   rsloper3max = rsloper2max * rslopermax
   rslopes3max = rslopes2max * rslopesmax





   xam_r = PI*denr/6.
   xbm_r = 3.
   xmu_r = 0.
   xam_s = PI*dens/6.
   xbm_s = 3.
   xmu_s = 0.
   xam_g = PI*dens/6.      
   xbm_g = 3.
   xmu_g = 0.

   call radar_init


  END SUBROUTINE wsm5init

      REAL FUNCTION rgmma(x)

  IMPLICIT NONE


      REAL :: euler
      PARAMETER (euler=0.577215664901532)
      REAL :: x, y
      INTEGER :: i
      if(x.eq.1.)then
        rgmma=0.
          else
        rgmma=x*exp(euler*x)
        do i=1,10000
          y=float(i)
          rgmma=rgmma*(1.000+x/y)*exp(-x/y)
        enddo
        rgmma=1./rgmma
      endif
      END FUNCTION rgmma


     subroutine effectRad_wsm5 (t, qc, qi, qs, rho, qmin, t0c,        &
                                re_qc, re_qi, re_qs, kts, kte, ii, jj)










      implicit none


      integer, intent(in) :: kts, kte, ii, jj
      real, intent(in) :: qmin
      real, intent(in) :: t0c
      real, dimension( kts:kte ), intent(in)::  t
      real, dimension( kts:kte ), intent(in)::  qc
      real, dimension( kts:kte ), intent(in)::  qi
      real, dimension( kts:kte ), intent(in)::  qs
      real, dimension( kts:kte ), intent(in)::  rho
      real, dimension( kts:kte ), intent(inout):: re_qc
      real, dimension( kts:kte ), intent(inout):: re_qi
      real, dimension( kts:kte ), intent(inout):: re_qs

      integer:: i,k
      integer :: inu_c
      real, dimension( kts:kte ):: ni
      real, dimension( kts:kte ):: rqc
      real, dimension( kts:kte ):: rqi
      real, dimension( kts:kte ):: rni
      real, dimension( kts:kte ):: rqs
      real :: temp
      real :: lamdac
      real :: supcol, n0sfac, lamdas
      real :: diai      
      logical :: has_qc, has_qi, has_qs

      real, parameter :: R1 = 1.E-12
      real, parameter :: R2 = 1.E-6

      real, parameter :: bm_r = 3.0
      real, parameter :: obmr = 1.0/bm_r
      real, parameter :: nc0  = 3.E8

      has_qc = .false.
      has_qi = .false.
      has_qs = .false.

      do k = kts, kte
        
        rqc(k) = max(R1, qc(k)*rho(k))
        if (rqc(k).gt.R1) has_qc = .true.
        
        rqi(k) = max(R1, qi(k)*rho(k))
        temp = (rho(k)*max(qi(k),qmin))
        temp = sqrt(sqrt(temp*temp*temp))
        ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
        rni(k)= max(R2, ni(k)*rho(k))
        if (rqi(k).gt.R1 .and. rni(k).gt.R2) has_qi = .true.
        
        rqs(k) = max(R1, qs(k)*rho(k))
        if (rqs(k).gt.R1) has_qs = .true.
      enddo

      if (has_qc) then
        do k=kts,kte
          if (rqc(k).le.R1) CYCLE
          lamdac   = (pidnc*nc0/rqc(k))**obmr
          re_qc(k) =  max(2.51E-6,min(1.5*(1.0/lamdac),50.E-6))
        enddo
      endif

     if (has_qi) then
        do k=kts,kte
          if (rqi(k).le.R1 .or. rni(k).le.R2) CYCLE
          diai = 11.9*sqrt(rqi(k)/ni(k))
          re_qi(k) = max(10.01E-6,min(0.75*0.163*diai,125.E-6))
        enddo
      endif

      if (has_qs) then
        do k=kts,kte
          if (rqs(k).le.R1) CYCLE
          supcol = t0c-t(k)
          n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          lamdas = sqrt(sqrt(pidn0s*n0sfac/rqs(k)))
          re_qs(k) = max(25.E-6,min(0.5*(1./lamdas), 999.E-6))
        enddo
      endif

      end subroutine effectRad_wsm5

                                                          
END MODULE module_mp_wsm5
