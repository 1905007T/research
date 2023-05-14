!----------------------------------------------------------------------
!        共通変数の定義
!----------------------------------------------------------------------
      module fdtd_variable

        integer::order,lpml                    !PMLの次数，層数
        real::rmax                             !要求精度〔dB〕
        parameter(lpml=64,order=4,rmax=-120.0)
  !解析領域の分割数，計算領域の分割数，時間ステップ総数
        integer::nz0,nz,nstep
        parameter(nz0=2000, nstep=5500)
        parameter(nz=nz0+2*lpml)
  !電磁界の配列，係数の配列，媒質定数の配列
        real(kind=8)::ex(0:nz),ey(0:nz),ez(0:nz),hy(0:nz)
        real(kind=8)::jx(0:nz),jy(0:nz),jz(0:nz),je(0:nz)
        real(kind=8)::ae(0:nz),be(0:nz),am(0:nz),bm(0:nz)
        real::dz,dt,t
        !parameter(dz=75e-6,dt=0.125e-12)
        real::epsd(-1:nz),sgmd(-1:nz),mud(-1:nz),msgmd(-1:nz)
  !定数
        real::eps0,mu0,z0,c,pi
        parameter(eps0=8.8541878e-12,mu0=1.2566371e-6)
        parameter(z0=376.73031)
        parameter(c=2.9979246e8)
        parameter(pi=3.1415926535)
        real::aj,wj,wl
        parameter(aj=5e4)
        real::zp,a                        !励振パルスのパラメータ
        real::rd,rr
        real,parameter::duration=0.1e-8,t0=4.0*duration  !�p���X���C�s�[�N���� 
        real(kind=8),parameter::qe=1.602176462d-19,mel = 9.10938188d-31
        !real(kind=8),parameter::wpe=(2*pi)*28.7e9,wce_t=(2*pi)*28.7e9,nue=20e8
        real(kind=8),parameter::wpe=(2*pi)*1e9,wce_t=(2*pi)*1e8,nue=(2*pi)*1e8
        real(kind=8) :: wce(3) !cyclotron frequency
        real(kind=8) :: imat(3,3),omat_e(3,3)
        real(kind=8) :: sa_e(3,3),sb_e(3,3),sab_e(3,3),tc_e(3,3),tu_e(3,3)
        real(kind=8) :: vx_e(0:nz),vy_e(0:nz),vz_e(0:nz)
        real(kind=8) :: flag(0:nz,3),nd_e(0:nz)
        real(kind=8) :: jx_e(0:nz),jy_e(0:nz),jz_e(0:nz)

        end module fdtd_variable
  !----------------------------------------------------------------------
  !        メインプログラム
  !----------------------------------------------------------------------
        program fdtd_1d
        use fdtd_variable                    !共通変数の引用
        integer::n,pls,damp
        character(50)::fname
        real(kind=8)::ee,he,hin,hre
        ee=0
        he=0
        hin=0
        hre=0
        pls=1                                !0:プラズマなし　1:非磁化プラズマ　2:磁化プラズマ
        damp=1
        if(damp==0) then
            rd=0.0
            rr=0.0
            write(*,*)'nodamp'
        else
            rd=1.0
            rr=1.0
            write(*,*)'damp'
        endif


        open(02,file='damp.csv')
        write(02,*)'t',',','hy'
        open(04,file='timevariation.csv')
        write(04,*)'t',',','ex'
        open(03,file='e1500.csv')
        write(03,*)'t',',','ex'
        open(06,file='damping.csv')
        write(06,*)'t,',',','damping'
        if(pls==0) then
            open(05,file='nopls.csv')
            write(*,*)'vacuum'
        else if(pls==1) then
            open(05,file='pls.csv')
            write(*,*)'non-magnetized plasma'
        else
            open(05,file='mpls.csv')
            write(*,*)'magnetized plasma'
        end if
        write(05,*)'t',',','ex'

        call setup                           !初期設定
        write(*,*)'dt=',dt,'wl/dz=',wl/dz,'duration=', duration
        !write(*,*)(2-nue*dt)/(2+nue*dt),2*eps0*(wpe**2)*dt/(2+nue*dt)
        t=dt                                 !時間
        do n=1,nstep
          call e_cal                         !電界の計算
          write(03,*)n,',',ex(1600)
          !call feed()
          !if(pls==0) then
          !  write(05,*)t,',',ex(400)
          !else
          !  write(05,*)t,',',ex(40)
          !end if
          if(mod(n,30)==0) then
              write(fname,"(i0,'.csv')")n
              open(n,file='./graph/'//fname)
              write(n,*)'k',',','ex'
              do k=0,nz
                  write(n,*)k,',',ex(k)
              end do
              close(n)
          end if
          do k=nz/2,nz
            ee=ee+ex(k)**2
          end do
          write(04,*)n,',',ee
          ee=0
          if(n==2248) then
            do k=nz/2,nz
                hin=hin+hy(k)**2
            end do
          end if
          if(n==4496) then
            do k=nz/2,nz
                hre=hre+hy(k)**2
            end do
          end if
          t=t+0.5*dt                         !時間を半ステップ前進
          call h_cal                         !磁界の計算
          if(pls==1) then
            call current
          else if(pls==2) then
            call currentb(n)
          end if
          t=t+0.5*dt                         !時間を半ステップ前進
        end do
        !write(*,*)wl,',',hre/hin
        close(02)
        close(03)
        close(04)
        close(05)
        end program fdtd_1d
  !----------------------------------------------------------------------
  !        解析モデルの設定と係数の計算
  !----------------------------------------------------------------------
        subroutine setup()
        use fdtd_variable                    !共通変数の引用
        real::eps,sgm,mu,msgm,a0,b
        
        d=0.1                                !誘電体層の厚さ[m]
        nd=50                                !誘電体層の分割数
        dz=d/nd                              !セルサイズ
        wl=dz*15
        wj=c/wl
        epsr=1.0                             !誘電体層の比誘電率
        kd=nz/2                              !誘電体層の左側の位置
        do k=-1,nz
           if(k<kd.or.k>=kd+nd) then         !真空領域
              epsd(k)=1.0                    !真空の比誘電率と導電率
              sgmd(k)=0.0
              mud(k)=1.0
              msgmd(k)=0.0
           else                              !誘電体層の領域
              epsd(k)=epsr                   !誘電体の比誘電率と導電率
              sgmd(k)=0.0
              mud(k)=1.0
              msgmd(k)=0.0
           end if
        end do
       
        !dt=0.99999*dz/c                       !時間ステップ
        dt=0.5000*dz/c
        do k=0,nz-1
           eps=0.5*(epsd(k)+epsd(k-1))*eps0   !セル表面の媒質定数は両側の平均
           sgm=0.5*(sgmd(k)+sgmd(k-1))
           b=dt/eps
           a0=0.5*sgm*b
           ae(k)=(1.0-a0)/(1.0+a0)            !電界の係数
           be(k)=b/(1.0+a0)/dz
           mu=mud(k)*mu0
           msgm=msgmd(k)
           b=dt/mu
           a0=0.5*msgm*b
           am(k)=(1.0-a0)/(1.0+a0)            !磁界の係数
           bm(k)=b/(1.0+a0)/dz
        end do

        zp=(lpml+100)*dz                      !パルスの中心位置
        a=100.0*dz                             !パルス幅
        !do k=0,nz                             !電界の初期値
           !z=k*dz
           !ex(k)=pulse(z,0.0)
        !end do
        !do k=0,nz-1                           !磁界の初期値
           !z=(k+0.5)*dz
           !hy(k)=pulse(z,0.5*dt)/z0
        !end do
        end subroutine setup
  !----------------------------------------------------------------------
  !        入射パルス  (z:位置    tm:時刻)
  !----------------------------------------------------------------------
        real function pulse(z,tm)
        use fdtd_variable
        
        !pulse=exp(-((z-zp-c*tm)/a)**2)       !ガウスパルス
        pulse=-0.000001*((z-zp-c*tm)/(a**3))*exp(-((z-zp-c*tm)**2/(2*(a**2))))
        end function pulse
!-----------------------------------------------------------------------
!       電流密度の計算
!-----------------------------------------------------------------------
        subroutine current()
            use fdtd_variable
            integer::k
            do k=0,nz
                jx(k)=(2-nue*dt)*jx(k)/(2+nue*dt)&
                &+2*eps0*(wpe**2)*dt*ex(k)/(2+nue*dt)
            end do
        end subroutine current
!-----------------------------------------------------------------------
!       電流密度の計算
!-----------------------------------------------------------------------
        subroutine currentb(n)
        use fdtd_variable
        integer::k,i,j
        integer,intent(in)::n
        integer :: n0,lwork,infom
        integer :: ipiv(3),work(192)
        real(kind=8) :: tmp(3,3)
        real(kind=8)::vxtmp(0:nz),vytmp(0:nz),vztmp(0:nz)
        n0 = 3
        lwork = n0*64

        if(n<2) then

        wce(1) = wce_t
        wce(2) = wce_t
        wce(3) = 0

        imat(1,1) = 1.0d0
        imat(2,2) = 1.0d0
        imat(3,3) = 1.0d0

        omat_e(1,1) = nue
        omat_e(1,2) = wce(3)
        omat_e(1,3) =-wce(2)

        omat_e(2,1) =-wce(3)
        omat_e(2,2) = nue 
        omat_e(2,3) = wce(1)

        omat_e(3,1) = wce(2)
        omat_e(3,2) =-wce(1)
        omat_e(3,3) = nue 

        omat_e = dt/2.0d0*omat_e
        do i=1,3
            do j=1,3
                sa_e(i,j) = imat(i,j) + omat_e(i,j)
            end do
        end do
        !sa_e = adding_mat(imat(1:3,1:3),omat_e(1:3,1:3))
        !sa_e = inversing_mat(sa_e(1:3,1:3))
        tmp = sa_e
        call dgetrf(n0,n0,sa_e,n0,ipiv,infom)
        call dgetri(n0,sa_e,n0,ipiv,work,lwork,infom)
        tmp = matmul(tmp,sa_e)
        !sa_e = adding_mat(imat(:,:),omat_e(:,:))
        !sa_e = inversing_mat(sa_e(:,:))


        do i=1,3
            do j=1,3
                sb_e(i,j) = imat(i,j) - omat_e(i,j)
            end do
        end do
        !sb_e = adding_mat(imat(1:3,1:3),-omat_e(1:3,1:3))

        sab_e = matmul(sa_e,sb_e)
        tc_e = qe*dt/mel*sa_e

        do k=0,nz
            nd_e(k) = mel*eps0*(wpe**2.0d0)/(qe**2.0d0)
            flag(k,1)=1.0d0
            flag(k,2)=1.0d0
            flag(k,3)=1.0d0
        end do

    endif

        do k=0,nz
            vxtmp(k)=vx_e(k)
            vytmp(k)=vy_e(k)
            vztmp(k)=vz_e(k)
        end do

        do k=0,nz
            vx_e(k) = sab_e(1,1)*vxtmp(k)+sab_e(1,2)*vytmp(k)+sab_e(1,3)*vztmp(k)-tc_e(1,1)*ex(k)-tc_e(1,2)*ey(k)-tc_e(1,3)*ez(k)
            vy_e(k) = sab_e(2,1)*vxtmp(k)+sab_e(2,2)*vytmp(k)+sab_e(2,3)*vztmp(k)-tc_e(2,1)*ex(k)-tc_e(2,2)*ey(k)-tc_e(2,3)*ez(k)
            vz_e(k) = sab_e(3,1)*vxtmp(k)+sab_e(3,2)*vytmp(k)+sab_e(3,3)*vztmp(k)-tc_e(3,1)*ex(k)-tc_e(3,2)*ey(k)-tc_e(3,3)*ez(k)
            vx_e(k) = vx_e(k)*flag(k,1)
            vy_e(k) = vy_e(k)*flag(k,2)
            vz_e(k) = vz_e(k)*flag(k,3)
        enddo

        do k=0,nz
            jx_e(k) = -qe*nd_e(k)*vx_e(k)*flag(k,1)
            jy_e(k) = -qe*nd_e(k)*vy_e(k)*flag(k,2)
            jz_e(k) = -qe*nd_e(k)*vz_e(k)*flag(k,3)

            jx(k) = jx_e(k)
            jy(k) = jy_e(k)
            jz(k) = jz_e(k)
        enddo

        contains
            function inversing_mat(x) result(ans)
                integer :: n,lwork,ifail,infom
                integer :: ipiv(3),work(192)
                integer :: i,j
                real(kind=8) :: ans(3,3),tmp(3,3)
                real(kind=8),intent(inout) :: x(3,3)
                n = 3
                lwork = n*64
                tmp = x
                call dgetrf(n,n,x,n,ipiv,infom)
                call dgetri(n,x,n,ipiv,work,lwork,infom)
                tmp = matmul(tmp,x)
                ans = x
            end function
        end subroutine currentb
!-----------------------------------------------------------------------
!       �d�����Ƃ��̓_�̓d�E
!-----------------------------------------------------------------------      
      subroutine feed()
      use fdtd_variable
      !real::iz
      real,parameter::dl=0.001

      !iz=exp(-((t-0.5*dt-t0)/duration)**2)              !�K�E�X�p���X
      !ex(ifed)=exp(-((t-t0)/duration)**2)/dl
      !ez(ifed,jfed)=ez(ifed,jfed)-befed*iz/(dx*dy)      !�d�E
      !ez(ifed,jfed)=exp(-((t-t0)/duration)**2)/dl
      ex(100)=(sin(wpe*t))**2
      end subroutine feed
  !----------------------------------------------------------------------
  !        電界の計算
  !----------------------------------------------------------------------
        subroutine e_cal()  
        use fdtd_variable
        integer::k

        je(nz/2)=-aj*exp(-((t-t0)/duration)**2)*sin(wj*t)/dt

        do k=1,nz-1
           ex(k)=fm(real(k),rd)*(ae(k)*ex(k)&
           &-fm(real(k),rr)*be(k)*(hy(k)-hy(k-1))-fm(real(k),rr)*dt*jx(k)/eps0)&
           &-fm(real(k),rd)*fm(real(k),rr)*dt*je(k)
        end do
        end subroutine e_cal
  !----------------------------------------------------------------------
  !        磁界の計算
  !----------------------------------------------------------------------
        subroutine h_cal() 
        use fdtd_variable      
        integer::k
  
        do k=0,nz-1
           hy(k)=fm(real(k+0.5),rd)*(am(k)*hy(k)-fm(real(k+0.5),rr)*bm(k)*(ex(k+1)-ex(k)))
        end do
        end subroutine h_cal
!----------------------------------------------------------------------
!        ?��?��?��?��?��̈�̌v?��Z
!----------------------------------------------------------------------
      real function fm(k,r)
        use fdtd_variable
        real::l,k,r

        l=k-(nz0/2+lpml)

        if(nz0/2<l) then
            fm=1-(r*(l-nz0/2)/lpml)**2
            !fm=1.0
        else if(l<-nz0/2) then
            fm=1-(r*(-l-nz0/2)/lpml)**2
            !fm=1.0
        else
            fm=1.0
        end if

        return

      end function fm

    function adding_mat(x,y) result(ans)
        use fdtd_variable
        implicit none 
        integer :: i,j
        real(kind=8),intent(in) :: x(3,3),y(3,3)
        real(kind=8):: ans(3,3)
        
        do i=1,3
            do j=1,3
                ans(i,j) = x(i,j) + y(i,j)
            enddo
        enddo
        
    end function 

    function inversing_mat(x) result(ans)
        use fdtd_variable
        implicit none
        real(kind=8) :: ans(3,3),tmp(3,3)
        real(kind=8),intent(inout) :: x(3,3)
        tmp = x
        tmp = matmul(tmp,x)
        ans = x
    end function