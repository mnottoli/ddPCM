      subroutine gradR(ISph,VPlm,VCos,VSin,BasLoc,dBsLoc,G,Y,FX)
      use ddcosmo
      implicit none
C
C     Compute the gradient of ddPCM R and contract it
C     < Y, grad R (PhiE - Phi) >
C
C     Physical quantities
      real*8 G(nbasis,*), Y(ngrid,*), FX(*)
C     Various scratch arrays
      Real*8 VIK(3), SIK(3), VKI(3), SKI(3), VKJ(3), SKJ(3), VJI(3),
     $  SJI(3), VA(3), VB(3), A(3), D(3)
C     Jacobian matrix
      Real*8 SJac(3,3)
C     Other scratch arrays
      Real*8 VPlm(*), VCos(*), VSin(*), BasLoc(*), dBsLoc(3,*)
      integer isph
C     indexes and scalars
      integer its, ik, ksph, l, m, ind, jsph, icomp, jcomp
      real*8 cx, cy, cz, vvki, tki, gg, fl, fac, vvkj, tkj
      real*8 tt, fcl, dij, fjj, gi, fii, vvji, tji, qji
      real*8 b, vvik, tik, qik, tlow, thigh, duj
C
      tlow  = one - pt5*(one - se)*eta
      thigh = one + pt5*(one + se)*eta
C
C     First set of contributions:
C     diagonal block, Kc and part of Kb
C
      Do 10 ITs = 1, ngrid
C
C       Sum over KSph in neighbors of ISph
        Do 20 IK = INL(ISph), INL(ISph+1) - 1
          KSph = NL(IK)
C
C         Build geometrical quantities
          cx = CSph(1,KSph) + RSph(KSph)*grid(1,ITs)
          cy = CSph(2,KSph) + RSph(KSph)*grid(2,ITs)
          cz = CSph(3,KSph) + RSph(KSph)*grid(3,ITs)
          VKI(1) = Cx - CSph(1,ISph)
          VKI(2) = Cy - CSph(2,ISph)
          VKI(3) = Cz - CSph(3,ISph)
          VVKI = Sqrt(VKI(1)*VKI(1) + VKI(2)*VKI(2) +
     $      VKI(3)*VKI(3))
          TKI  = VVKI/RSph(ISph)
C
C         Contributions involving grad I of UK come from the switching
C         region.
C         Note: ui avoids contributions from points that are in the
C         switching between ISph and KSph but are buried in a third
C         sphere.
C
          If ((TKI.gt.TLow).and.(TKI.lt.THigh) .and.
     $      ui(ITs,KSph).gt.Zero) Then
C
C           Other geometrical quantities
            SKI = VKI/VVKI
CCC            write(6,*) VKI(1), VKI(2), VKI(3), VVKI, RSph(ISph), TKI,
CCC     $        SKI(1), SKI(2), SKI(3)
C
C           Diagonal block KK contribution, with K in N(I)
            GG = Zero
            Do 21 L = 0, lmax
              Ind = L*L + L + 1
              fL = Dble(L)
              Fac = Two*Pi/(Two*fL + One)
              Do 21 M = -L, L 
                GG = GG + Fac*basis(Ind+M,ITs)*G(Ind+M,KSph)
CCC               Write(6,*) Fac, basis(Ind+M,ITs), G(Ind+M,KSph)
  21            Continue
CCC           Write(6,*) GG
C
C           Kc contribution
            Do 120 JSph = 1, NSph
              If (JSph.ne.KSph .and. JSph.ne.ISph) Then 
                VKJ(1) = Cx - CSph(1,JSph)
                VKJ(2) = Cy - CSph(2,JSph)
                VKJ(3) = Cz - CSph(3,JSph)
                VVKJ = Sqrt(VKJ(1)*VKJ(1) + VKJ(2)*VKJ(2) +
     $            VKJ(3)*VKJ(3))
                TKJ  = VVKJ/RSph(JSph)
                SKJ  = VKJ/VVKJ
                Call YlmBas(SKJ,BasLoc,VPlm,VCos,VSin)
                TT = One/TKJ
                Do 121 L = 0, lmax
                  Ind = L*L + L + 1
                  FcL = - Four*PI*Dble(L)/(Two*Dble(L)+One)*TT
                  Do 122 M = -L, L
                    GG = GG + FcL*G(Ind+M,JSph)*BasLoc(Ind+M)
  122               Continue
                  TT = TT/TKJ
  121             Continue
                End If
  120         Continue
CCC            Write(6,*) GG
C
C
C           Part of Kb contribution
            Call YlmBas(SKI,BasLoc,VPlm,VCos,VSin)
            TT = One/TKI
            Do 251 L = 0, lmax
              Ind = L*L + L + 1
              FcL = - Four*Pi*Dble(L)/(Two*Dble(L)+One)*TT
              Do 252 M = -L, L
                GG = GG + FcL*G(Ind+M,ISph)*BasLoc(Ind+M)
  252           Continue
              TT = TT/TKI
  251         Continue
CCC           Write(6,*) GG
C
C
C           Common step, product with grad I UJ
            dUJ = dFSW(TKI,se,eta)/RSph(ISph)
            FJJ = dUj*w(ITs)*GG*Y(ITs,KSph)
            write(6,'(7f10.5)') dUj, w(ITs), GG, Y(ITs,KSph), SKI(1),
     $      SKI(2), SKI(3)
            FX(1) = FX(1) - FJJ*SKI(1)
            FX(2) = FX(2) - FJJ*SKI(2)
            FX(3) = FX(3) - FJJ*SKI(3)
            End If
   20     Continue
C
C       Diagonal block II contribution
        If (ui(ITs,ISph).gt.Zero.and.ui(ITs,ISph).lt.One) Then
          Gi = Zero
          Do 30 L = 0, lmax
            Ind = L*L + L + 1
            fL = Dble(L)
            Fac = Two*Pi/(Two*fL + One)
            Do 30 M = -L, L 
              Gi = Gi + Fac*basis(Ind+M,ITs)*G(Ind+M,ISph)
   30         Continue
          FII = w(ITs)*Gi*Y(ITs,ISph)
          FX(1) = FX(1) + FII*zi(1,ITs,ISph)
          FX(2) = FX(2) + FII*zi(2,ITs,ISph)
          FX(3) = FX(3) + FII*zi(3,ITs,ISph)
          End If
   10   Continue
C
C     Second set of contributions:
C     part of Kb and Ka
C
      Do 200 ITs = 1, ngrid
C
C       Run over all the spheres except ISph 
        Do 210 JSph = 1, NSph
          If (ui(ITs,JSPh).gt.Zero .and. JSph.ne.ISph) Then
C
C           Build geometrical quantities
            Cx = CSph(1,JSph) + RSph(JSph)*grid(1,ITs)
            Cy = CSph(2,JSph) + RSph(JSph)*grid(2,ITs)
            Cz = CSph(3,JSph) + RSph(JSph)*grid(3,ITs)
            VJI(1) = Cx - CSph(1,ISph)
            VJI(2) = Cy - CSph(2,ISph)
            VJI(3) = Cz - CSph(3,ISph)
            VVJI = Sqrt(VJI(1)*VJI(1) + VJI(2)*VJI(2) +
     $        VJI(3)*VJI(3))
            TJI = VVJI/RSph(ISph)
            QJI = One/VVJI
            SJI = VJI/VVJI
C
C           Build the jacobian of SJI
            SJac = Zero
            SJac(1,1) = - One
            SJac(2,2) = - One
            SJac(3,3) = - One
            Do 211 IComp = 1, 3
              Do 211 JComp = 1, 3
                SJac(IComp,JComp) = QJI*(SJac(IComp,JComp)
     $            + SJI(IComp)*SJI(JComp))
  211           Continue
C
C           Assemble the local basis and its gradient
            Call DBasis(SJI,BasLoc,dBsLoc,VPlm,VCos,VSin)
C
C           Assemble the contribution
            A = Zero
            TT = One/(TJI)
            Do 220 L = 0, lmax
              Ind = L*L + L + 1
              fL = Dble(L)
              FcL = - TT*Four*Pi*fL/(Two*fL + One)
              Do 221 M = -L, L
                Fac = FcL*G(Ind+M,ISph)
                B = (fL + One)*BasLoc(Ind+M)/(RSph(ISph)*TJI)
C
C               Apply the jacobian to grad Y
                VA(1) = SJac(1,1)*dBsLoc(1,Ind+M) +
     $           SJac(1,2)*dBsLoc(2,Ind+M) + SJac(1,3)*dBsLoc(3,Ind+M)
                VA(2) = SJac(2,1)*dBsLoc(1,Ind+M) +
     $           SJac(2,2)*dBsLoc(2,Ind+M) + SJac(2,3)*dBsLoc(3,Ind+M)
                VA(3) = SJac(3,1)*dBsLoc(1,Ind+M) +
     $           SJac(3,2)*dBsLoc(2,Ind+M) + SJac(3,3)*dBsLoc(3,Ind+M)
                A(1) = A(1) + Fac*(SJI(1)*B + VA(1))
                A(2) = A(2) + Fac*(SJI(2)*B + VA(2))
                A(3) = A(3) + Fac*(SJI(3)*B + VA(3))
  221           Continue
              TT = TT/TJI
  220         Continue
            Fac = ui(ITs,JSph)*w(ITs)*Y(ITs,JSph)
            FX(1) = FX(1) - Fac*A(1)
            FX(2) = FX(2) - Fac*A(2)
            FX(3) = FX(3) - Fac*A(3)
            End If
  210     Continue
  200   Continue
C
C     Ka contribution
C
      Do 300 ITs = 1, ngrid
        Cx = CSph(1,ISph) + RSph(ISph)*grid(1,ITs)
        Cy = CSph(2,ISph) + RSph(ISph)*grid(2,ITs)
        Cz = CSph(3,ISph) + RSph(ISph)*grid(3,ITs)
        A = Zero
C
C       Iterate on all the spheres except ISph
        Do 310 KSph = 1, NSph
          If (ui(ITs,ISph).gt.Zero .and. KSph.ne.ISph) Then
C
C           Geometrical stuff
            VIK(1) = Cx - CSph(1,KSph)
            VIK(2) = Cy - CSph(2,KSph)
            VIK(3) = Cz - CSph(3,KSph)
            VVIK = Sqrt(VIK(1)*VIK(1) + VIK(2)*VIK(2) + 
     $        VIK(3)*VIK(3))
            TIK = VVIK/RSph(KSph)
            QIK = One/VVIK
            SIK = VIK/VVIK
C
C           Build the jacobian of SIK
            SJac = Zero
            SJac(1,1) = One
            SJac(2,2) = One
            SJac(3,3) = One
            Do 311 IComp = 1, 3
              Do 311 JComp = 1, 3
                SJac(IComp,JComp) = QIK*(SJac(IComp,JComp)
     $            - SIK(IComp)*SIK(JComp))
  311           Continue
C
C           If we are in the switching region, recover grad_i U_i
C           Not the best way to do this
            VB = Zero
            If (ui(ITs,ISph).lt.One) Then
              VB(1) = zi(1,ITs,ISph)
              VB(2) = zi(2,ITs,ISph)
              VB(3) = zi(3,ITs,ISph)
            End If
C
C           Assemble the local basis and its gradient
            Call DBasis(SIK,BasLoc,dBsLoc,VPlm,VCos,VSin)
C
C           Assemble the contribution
            TT = One/(TIK)
            Do 320 L = 0, lmax
              Ind = L*L + L + 1
              fL = DBle(L)
              FcL = - TT*Four*Pi*fL/(Two*fL + One)
              Do 321 M = -L, L
                Fac = FcL*G(Ind+M,KSph)
                Fac = - Fac*BasLoc(Ind+M)
                A(1) = A(1) + Fac*VB(1)
                A(2) = A(2) + Fac*VB(2) 
                A(3) = A(3) + Fac*VB(3)
C
                Fac = ui(ITs,ISph)*FcL*G(Ind+M,KSph)
                B = - (fL + One)*BasLoc(Ind+M)/(RSph(KSph)*TIK)
C
C               Apply the jacobian to grad Y
                VA(1) = SJac(1,1)*dBsLoc(1,Ind+M) +
     $           SJac(1,2)*dBsLoc(2,Ind+M) + SJac(1,3)*dBsLoc(3,Ind+M)
                VA(2) = SJac(2,1)*dBsLoc(1,Ind+M) +
     $           SJac(2,2)*dBsLoc(2,Ind+M) + SJac(2,3)*dBsLoc(3,Ind+M)
                VA(3) = SJac(3,1)*dBsLoc(1,Ind+M) +
     $           SJac(3,2)*dBsLoc(2,Ind+M) + SJac(3,3)*dBsLoc(3,Ind+M)
                A(1) = A(1) + Fac*(SIK(1)*B + VA(1))
                A(2) = A(2) + Fac*(SIK(2)*B + VA(2))
                A(3) = A(3) + Fac*(SIK(3)*B + VA(3))
  321           Continue 
              TT = TT/TIK
  320         Continue
            End If
  310     Continue
        Fac = w(ITs)*Y(ITs,ISph)
        FX(1) = FX(1) - Fac*A(1)
        FX(2) = FX(2) - Fac*A(2)
        FX(3) = FX(3) - Fac*A(3)
  300   Continue
C
      Return
      End
