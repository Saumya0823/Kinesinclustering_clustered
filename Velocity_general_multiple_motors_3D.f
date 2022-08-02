	Program cargonunderexternalforce
	parameter (N=10)
	real Lo
        dimension tau_ran(3),angle(3),copy_beadmx(N),
     1          copy_beadmz(N),copy_headmx(N),copy_headmy(N),
     1          copy_headmz(N),copy_beadmy(N),x(10),y(10),
     1		sort(10),store_dx(4)
	dimension beadmotx(N),beadmoty(N),beadmotz(N),
     1		headmotx(N),headmoty(N),headmotz(N),
     1          dup_motbx(N),dup_motby(N),dup_motbz(N),
     1          dup_mothx(N),dup_mothy(N),dup_mothz(N),
     1          rel_motbx(N),rel_motby(N),rel_motbz(N),
     1          rel_mothx(N),rel_mothy(N),rel_mothz(N),
     1          rel_range(N),stalk_motl(N),range(N),
     1		dup_stalk(N),sum_torq(200),del(200,200)
c     1          ang_mot(N,3),Fall_mot(N),Fall_motx(N),
c     1          Fall_moty(N),Fall_motz(N),taual1(N),taual2(N),
c     1          taual3(N)
	dimension thetamot(N),phimot(N),StepR(N),DetR(N),AttR(N),
     1		F_motx(N),F_moty(N),F_motz(N),F_mot(N),Fonmot(N),
     1		u_motx(N),u_moty(N),u_motz(N),bead_MTdist(N),
     1		tau(200,200),ps(200),pd(200),pa(200),
     1		anchorx(N),anchory(N),anchorz(N)
        parameter (Fs=5.,iTemp=300,springk= 320,Bindingrate=5,w=2.,
     1          d=0.008,Vo=1.0,Rso=Vo/d,Rdo=1,Fd=3.0,Lo=0.110,
     1          boltzmann=1.3806e-5,pi=4*atan(1.),
     1          beadRadius=0.250,Rcube= beadRadius**3,
     1          twokbt= 2*boltzmann*iTemp,viscosity=1e-3,twopi=2*pi,
     1          gamnew= 6*pi*viscosity*beadRadius,
     1          gamma= 8*pi*viscosity*Rcube, del_t=1e-5, maxconf=1000)

	open(unit=21,file="Rate_V_clus_N10.dat")
        open(unit=23,file="V_clus_N10.dat")

c        SA= 4*pi*(beadRadius**2)
c        avgden= 9.7444*1e-5
c        N= SA*avgden
c	iZ = (N-1)/2
        variate= del_t/gamma
        vari2= del_t/gamnew

	do imtd= 0,0

	average_m=0
        do iF= 0,10,1								! force loop starts
          F_ext= -float(iF)
          Fx_ext= F_ext
	  Fy_ext=0; Fz_ext=0
	  print*, Fx_ext, Fy_ext, Fz_ext

	sigconf2=0
        sigconf=0; sum_act=0; sum_t=0
        sum_det=0; sum_att=0; sum_step=0
        sumsq_mot=0; tsq_sum=0
        detsq_sum=0; attsq_sum=0; stepsq_sum=0
	vel_sum=0; vel_sqsum=0
        print*, "conf loop starts"

        do 								!iconf=1,maxconf                                                            ! conf loop starts
	  print*, "conf=", sigconf2, beadRadius, N, imtd, iF
c	  write(10,*)"new configuration starts"
	  beadinitx= 0
	  beadinity= 0
	  beadinitz= beadRadius + imtd
c	  write(10,*) "cargo=", beadinitx, beadinity, beadinitz

	  do i=1,N
	   beadmotx(i)=0; beadmoty(i)=0; beadmotz(i)=0
	   headmotx(i)=0; headmoty(i)=0; headmotz(i)=0
	   thetamot(i)=0; phimot(i)=0
	   stalk_motl(i)=0
	  end do

c         do i=1,N
c1          phimot(i)= (pi/18)*rand()+ ((16*pi)/18)
c           thetamot(i)=  (pi/18)*rand()+((7*pi)/18)
c           phimot(i)= (pi/6)*rand()+ ((4*pi)/6)
c           thetamot(i)=  (pi/6)*rand()+((7*pi)/6)
c           caviar1= beadRadius*cos(thetamot(i))*sin(phimot(i))
c           beadmotx(i)= beadinitx + caviar1
c           caviar2= beadRadius*sin(thetamot(i))*sin(phimot(i))
c           beadmoty(i)= beadinity + caviar2
c           beadmotz(i)= beadinitz + beadRadius*cos(phimot(i))
c           term= sqrt((beadmotx(i)-beadinitx)**2+
c     1          (beadmoty(i)-beadinity)**2+
c     1          (beadmotz(i)-beadinitz)**2)
c           term_diff= abs(term-beadRadius)
c          if(term_diff .gt. 1e-4)then
c                go to 1
c           end if
c           write(*,*)beadmotx(i),beadmoty(i),beadmotz(i)
c           write(11,*)beadmotx(i),beadmoty(i),beadmotz(i)
c         end do
c
c         do i=iZ+1,N
c600        phimot(i)= (pi/18)*rand()+ ((16*pi)/18)
c           thetamot(i)=  (pi/18)*rand()+((16*pi)/18)
c600       phimot(i)= (pi/6)*rand()+ ((4*pi)/6)
c          thetamot(i)=  (pi/6)*rand()+((10*pi)/6)
c           caviar1= beadRadius*cos(thetamot(i))*sin(phimot(i))
c           beadmotx(i)= beadinitx + caviar1
c           caviar2= beadRadius*sin(thetamot(i))*sin(phimot(i))
c           beadmoty(i)= beadinity + caviar2
c           beadmotz(i)= beadinitz + beadRadius*cos(phimot(i))
c           term= sqrt((beadmotx(i)-beadinitx)**2+
c     1          (beadmoty(i)-beadinity)**2+
c     1          (beadmotz(i)-beadinitz)**2)
c           term_diff= abs(term-beadRadius)
c          if(term_diff .gt. 1e-4)then
c               go to 600
c           end if
c         end do

         do i=1,N
1          phimot(i)= acos(2*rand()-1)
           thetamot(i)= twopi*rand()
           caviar1= beadRadius*cos(thetamot(i))*sin(phimot(i))
           beadmotx(i)= beadinitx + caviar1
           caviar2= beadRadius*sin(thetamot(i))*sin(phimot(i))
           beadmoty(i)= beadinity + caviar2
           beadmotz(i)= beadinitz + beadRadius*cos(phimot(i))
           term= sqrt((beadmotx(i)-beadinitx)**2+
     1          (beadmoty(i)-beadinity)**2+
     1          (beadmotz(i)-beadinitz)**2)
           term_diff= abs(term-beadRadius)
           if(term_diff .gt. 1e-4)then
                go to 1
           end if
c          b_dis= sqrt((beadmoty(i)**2)+(beadmotz(i)**2))
c          if(b_dis .gt. Lo)go to 1
          b_dis= 1000*sqrt((beadmoty(i)**2)+(beadmotz(i)**2))
          if(b_dis .gt. 6.125)go to 1
c          write(10,*)beadmotx(i),beadmoty(i),beadmotz(i)
          end do

          do i=1,N
           range(i)=0
          end do

        do i=1,N                                                                        !Lattice Microtubule
           bead_MTdist(i)= sqrt((beadmotz(i)**2))
c           print*, "dist=", bead_MTdist(i)
           if(bead_MTdist(i) .le. Lo)then
             range(i)=1
            if(abs(beadmoty(i)) .gt. 0.012)then
                i_mot= int((beadmotx(i)*1000)/8)
                icheck= i_mot*8
                icheck1= (i_mot+1)*8
                xcheck= sqrt((abs(beadmoty(i))-0.012)**2+
     1          (beadmotx(i)-(0.001*icheck))**2+beadmotz(i)**2)
                xcheck1= sqrt((abs(beadmoty(i))-0.012)**2+
     1          (beadmotx(i)-(0.001*icheck1))**2+beadmotz(i)**2)
                if(xcheck .lt. xcheck1)then
                 if(xcheck .le. Lo)then
                  if(beadmoty(i) .lt. 0)headmoty(i)=-0.012
                  if(beadmoty(i) .gt. 0)headmoty(i)=0.012
                  headmotx(i)= (0.001*icheck); headmotz(i)=0
                 else
                  headmotx(i)=beadmotx(i)
                  headmoty(i)=beadmoty(i)
                  headmotz(i)=beadmotz(i)
                 end if
                else
                 if(xcheck1 .le. Lo)then
                  if(beadmoty(i) .lt. 0)headmoty(i)=-0.012
                  if(beadmoty(i) .gt. 0)headmoty(i)=0.012
                  headmotx(i)= (0.001*icheck1); headmotz(i)=0
                 else
                  headmotx(i)=beadmotx(i)
                  headmoty(i)=beadmoty(i)
                  headmotz(i)=beadmotz(i)
                 end if
                end if
            else if(abs(beadmoty(i)) .lt. 0.012)then
                i_mot= int((beadmotx(i)*1000)/8)
                i_moty= int((beadmoty(i)*1000)/4)
                icheck = i_mot*8
                icheck1= (i_mot+1)*8
                ichecky= i_moty*4
                ichecky1= (i_moty+1)*4
                sort(1)= sqrt((beadmoty(i)-(0.001*ichecky))**2+
     1          (beadmotx(i)-(0.001*icheck))**2+beadmotz(i)**2)
                x(1)=icheck; y(1)=ichecky
                sort(2)= sqrt((beadmoty(i)-(0.001*ichecky1))**2+
     1          (beadmotx(i)-(0.001*icheck))**2+beadmotz(i)**2)
                x(2)= icheck; y(2)= ichecky1
                sort(3)= sqrt((beadmoty(i)-(0.001*ichecky))**2+
     1          (beadmotx(i)-(0.001*icheck1))**2+beadmotz(i)**2)
                x(3)= icheck1; y(3)= ichecky
                sort(4)= sqrt((beadmoty(i)-(0.001*ichecky1))**2+
     1          (beadmotx(i)-(0.001*icheck1))**2+beadmotz(i)**2)
                x(4)= icheck1; y(4)= ichecky1
                do l=1,4
                 do j=1,l
                   if(sort(l) .lt. sort(j))then
                    tempor= sort(l)
                    sort(l)= sort(j)
                    sort(j)= tempor
                    temporx= x(l)
                    x(l)=x(j)
                    x(j)=temporx
                    tempory= y(l)
                    y(l)=y(j)
                    y(j)=tempory
                   end if
                 end do
                end do
                if(sort(1) .le. Lo)then
                  headmoty(i)= y(1)*0.001
                  headmotx(i)= x(1)*0.001; headmotz(i)=0
                else
                  headmotx(i)=beadmotx(i)
                  headmoty(i)=beadmoty(i)
                  headmotz(i)=beadmotz(i)
                end if
             end if
            else
              headmotx(i)=beadmotx(i)
              headmoty(i)=beadmoty(i)
              headmotz(i)=beadmotz(i)
            end if
        end do

          do i=1,N
           do j=1,i
           if((headmotx(i)==headmotx(j)).and.(i.ne.j))then
	     if(headmoty(i)==headmoty(j))then
             headmotx(i)=beadmotx(i)
             headmoty(i)=beadmoty(i)
             headmotz(i)=beadmotz(i)
	     end if
           end if
           end do
           s1= (headmotx(i)-beadmotx(i))**2+
     1          (headmoty(i)-beadmoty(i))**2+
     1          (headmotz(i)-beadmotz(i))**2
           stalk_motl(i)= sqrt(s1)
c           if(sigconf==3)print*, i, stalk_motl(i)
          end do
 
	  do i=1,N
	    rel_motbx(i)=0; rel_motby(i)=0; rel_motbz(i)=0
	    rel_mothx(i)=0; rel_mothy(i)=0; rel_mothz(i)=0
	    rel_range(i)=0
	    dup_motbx(i)=0; dup_motby(i)=0; dup_motbz(i)=0
	    dup_mothx(i)=0; dup_mothy(i)=0; dup_mothz(i)=0
	    dup_stalk(i)=0 
	  end do

	  k=0
	  do i=1,N
	    if(range(i).eq.1 .or. stalk_motl(i).gt.0)then
	      k=k+1
	      rel_motbx(k)= beadmotx(i)
	      rel_motby(k)= beadmoty(i)
	      rel_motbz(k)= beadmotz(i)
	      rel_mothx(k)= headmotx(i)
	      rel_mothy(k)= headmoty(i)
	      rel_mothz(k)= headmotz(i)
	      rel_range(k)= range(i)
	      dup_stalk(k)= stalk_motl(i)
c	      write(10,*) i, k, rel_range(k)
	     end if
	  end do
	  m=k

c	  write(10,*)"initial motor orientation"
	  do i=1,m
	   dup_motbx(i)= rel_motbx(i)
	   dup_motby(i)= rel_motby(i)
	   dup_motbz(i)= rel_motbz(i)
	   dup_mothx(i)= rel_mothx(i)
	   dup_mothy(i)= rel_mothy(i)
	   dup_mothz(i)= rel_mothz(i)
c	   write(10,*)dup_mothx(i)-dup_motbx(i),dup_mothx(i)
	  end do

	  v=0
	  do i=1,N
	    if(stalk_motl(i).ne.0)v=v+1
	  end do
	  average_m= average_m+ v
	  print*, "no of attached motors=", m

c	  write(10,*)"no of attached motors=",m
	  if(m==0)then
		go to 300
	  else
		sigconf= sigconf +1
	  end if 
	  print*, "m=", m
c------------------------------------------------------------------
	  final_cargox= beadinitx
	  final_cargoy= beadinity
	  final_cargoz= beadinitz
          copy_cargox= final_cargox
          copy_cargoy= final_cargoy
          copy_cargoz= final_cargoz

	  do i=1,N								! flushing of copy and final angular position arrays of motors' attachment points
	   copy_beadmx(i)=0; copy_beadmy(i)=0
	   copy_beadmz(i)=0; copy_headmx(i)=0
	   copy_headmy(i)=0; copy_headmz(i)=0
	  end do
	  do i=1,4
	    store_dx(i)=0
	  end do

          do i=1,N                                                              ! flushing of copy and final angular position arrays of motors' attachment points
           copy_beadmx(i)= beadmotx(i) 
	   copy_beadmy(i)= beadmoty(i)
           copy_beadmz(i)= beadmotz(i)
	   copy_headmx(i)= headmotx(i)
           copy_headmy(i)= headmoty(i)
	   copy_headmz(i)= headmotz(i)
          end do

	  print*, "time loop starts",sigconf2
	  Rs=0;mcount=0
	  attach=0; detach=0; step=0
	  time=0; act_mot=m
	  do 									! time loop starts
	    copy_cargox= final_cargox
	    copy_cargoy= final_cargoy
	    copy_cargoz= final_cargoz
c	    write(10,*)"---------------------------------"
c	    write(10,*) "new time configuration-"
c	    write(10,*)"time", time, "cargo pos=", copy_cargox,
c     1		copy_cargoy, copy_cargoz	
c	    print*, time

	    do j=1,m
	     dup_motbx(j)= rel_motbx(j)
	     dup_motby(j)= rel_motby(j)
	     dup_motbz(j)= rel_motbz(j)
	     dup_mothx(j)= rel_mothx(j)
	     dup_mothy(j)= rel_mothy(j)
	     dup_mothz(j)= rel_mothz(j)
	    end do

c	    write(10,*)"time=",time
            do i=1,N
             copy_beadmx(i)= beadmotx(i)
             copy_beadmy(i)= beadmoty(i)
             copy_beadmz(i)= beadmotz(i)
             copy_headmx(i)= headmotx(i)
             copy_headmy(i)= headmoty(i)
             copy_headmz(i)= headmotz(i)
c            write(10,*)copy_beadmx(i),copy_beadmy(i),copy_beadmz(i)
c            write(10,*)copy_headmx(i),copy_headmy(i),copy_headmz(i)
            end do
            
c	    write(10,*)"stalk length of relevant motors="
c            do i=1,m
c	      s1= (dup_motbx(i)-dup_mothx(i))**2+
c     1			(dup_motby(i)-dup_mothy(i))**2+
c     1			(dup_motbz(i)-dup_mothz(i))**2
c              dup_stalk(i)= sqrt(s1)
c	      write(10,*)i, sigconf2,time, dup_stalk(i)
c            end do
c
c	    write(10,*)"stalk length of all motors="
c	    do i=1,N
c		s1= (copy_beadmx(i)-copy_headmx(i))**2+
c     1			(copy_beadmy(i)-copy_headmy(i))**2+
c     1			(copy_beadmz(i)-copy_headmz(i))**2
c		stalk_motl(i)= sqrt(s1)
c		print*, "here",i, stalk_motl(i)
c	        write(10,*)i, stalk_motl(i)
c	    end do

            do i=1,N
                F_mot(i)=0; Fonmot(i)=0
                F_motx(i)=0;F_moty(i)=0
                F_motz(i)=0
            end do

c	    write(10,*)"force by relevant motors="
            do i=1,m
             if (dup_stalk(i) .le. Lo)then                                            ! case of detached/unattached motor
               F_mot(i)=0
             else
               F_mot(i)= (dup_stalk(i)-Lo)*springk                            	       ! magnitude of force exerted by attached and extended motors
	     end if
	     if((dup_mothx(i)-dup_motbx(i)).lt.0)then
		  Fonmot(i)= F_mot(i)
	     end if
	     if((dup_mothx(i)-dup_motbx(i)).ge.0)then
	  	  Fonmot(i)= -F_mot(i)
             end if
c	     write(10,*)i, dup_mothx(i)-dup_motbx(i)
c	     write(10,*) i, F_mot(i), Fonmot(i)
            end do

c            do i=1,m
c                del(i,1)=0; del(i,2)=0; del(i,3)=0
c                tau(i,1)=0; tau(i,2)=0; tau(i,3)=0
c                F_motx(i)=0; F_moty(i)=0; F_motz(i)=0
c                u_motx(i)=0; u_moty(i)=0; u_motz(i)=0
c            end do
c
c           write(10,*)"force and its component by motors="
c            sumF_x=0; sumF_y=0; sumF_z=0
c            sum_torq(1)=0; sum_torq(2)=0; sum_torq(3)=0
c
c           write(10,*)"sum of force by motors="
c            do i=1,m
c             if (F_mot(i).eq. 0)then                                            ! case of detached/unattached motor
c                F_motx(i)=0; F_moty(i)=0; F_motz(i)=0
c                del(i,1)=0; del(i,2)=0; del(i,3)=0
c                u_motx(i)=0; u_moty(i)=0; u_motz(i)=0
c                tau(i,1)=0; tau(i,2)=0; tau(i,3)=0
c             else if(F_mot(i).ne.0)then
c               F_mot(i)= (dup_stalk(i)-Lo)*springk                            ! magnitude of force exerted by attached and extended motors
c                del(i,1)= dup_motbx(i)- copy_cargox
c                del(i,2)= dup_motby(i)- copy_cargoy
c                del(i,3)= dup_motbz(i)- copy_cargoz
c                u_motx(i)= (dup_mothx(i)-dup_motbx(i))/dup_stalk(i)                 ! unit vector of x,y,z line connecting head and tail of motors
c                u_moty(i)= (dup_mothy(i)-dup_motby(i))/dup_stalk(i)
c                u_motz(i)= (dup_mothz(i)-dup_motbz(i))/dup_stalk(i)
c                F_motx(i)= F_mot(i)*u_motx(i)
c                F_moty(i)= F_mot(i)*u_moty(i)
c                F_motz(i)= F_mot(i)*u_motz(i)
c		tau(i,1)= del(i,2)*F_motz(i)-del(i,3)*F_moty(i)
c	        tau(i,2)= -del(i,1)*F_motz(i)+del(i,3)*F_motx(i)
c                tau(i,3)= del(i,1)*F_moty(i)-del(i,2)*F_motx(i)
c             end if
c
c             sumF_x= sumF_x + F_motx(i)
c             sumF_y= sumF_y + F_moty(i)
c             sumF_z= sumF_z + F_motz(i)
c             sum_torq(1)= sum_torq(1) + tau(i,1)
c             sum_torq(2)= sum_torq(2) + tau(i,2)
c             sum_torq(3)= sum_torq(3) + tau(i,3)
c            end do
c           write(10,*) sumF_x, sumF_y, sumF_z
c           write(10,*) sum_torq(1), sum_torq(2), sum_torq(3)
c
c3           r1= 2*rand()-1
c            r2= 2*rand()-1
c            gauss= r1**2 + r2**2
c            if(gauss.ge.1 .or. gauss.eq.0)go to 3
c            var= r1*sqrt(-2*log(gauss)/gauss)
c            x_random= var*sqrt(twokbt*vari2)
c
c9           b1= 2*rand()-1
c            b2= 2*rand()-1
c            gem= b1**2 + b2**2
c            if(gem.ge.1 .or. gem.eq.0)go to 9
c            vr1= b1*sqrt(-2*log(gem)/gem)
c            y_random= vr1*sqrt(twokbt*vari2)
c
c8           c1= 2*rand()-1
c            c2= 2*rand()-1
c            gas= c1**2 + c2**2
c            if(gas.ge.1 .or. gas.eq.0)go to 8
c            var2= c1*sqrt(-2*log(gas)/gas)
c            z_random= var2*sqrt(twokbt*vari2)
c	    write(10,*)"random force="
c           write(10,*) Fran_x, Fran_y, Fran_z
c
c2           d1= 2*rand()-1
c            d2= 2*rand()-1
c            fim= d1**2+d2**2
c            if(fim.ge.1 .or. fim.eq.0)go to 2
c            var3= d1*sqrt(-2*log(fim)/fim)
c            tau_ran(1)= var3*sqrt(twokbt*variate)
c
c6           e1= 2*rand()-1
c            e2= 2*rand()-1
c            tim= e1**2+ e2**2
c            if(tim.ge.1 .or. tim.eq.0)go to 6
c            var4= e1*sqrt(-2*log(tim)/tim)
c            tau_ran(2)= var4*sqrt(twokbt*variate)
c
c4           f1= 2*rand()-1
c            f2= 2*rand()-1
c            bim= f1**2+f2**2
c            if(bim.ge.1 .or. bim.eq.0)go to 4
c            var5= f1*sqrt(-2*log(bim)/bim)
c            tau_ran(3)= var5*sqrt(twokbt*variate)
c
c            write(10,*)"random torque="
c            write(10,*)tau_ran(1), tau_ran(2), tau_ran(3)
c
c            C_MT_dist= sqrt(copy_cargoy**2+copy_cargoz**2)                      ! when cargo volume gets cut by MT (steric spring force acting here)
c	    if(C_MT_dist .gt. beadRadius)then
c		cross_dist=0
c		Fsteric_x=0; Fsteric_y=0; Fsteric_z=0
c            else if(C_MT_dist .le. beadRadius)then                
c		cross_dist= abs(beadRadius-C_MT_dist)
c                Fsteric= -sterick*(C_MT_dist)
c                unit_dispx= 0
c                unit_dispy=-(copy_cargoy*cross_dist)/C_MT_dist
c                unit_dispz=-(copy_cargoz*cross_dist)/C_MT_dist
c		sub_dist= sqrt(unit_dispy**2+unit_dispz**2)
c                Fsteric_x= 0
c                Fsteric_y= Fsteric*(unit_dispy/sub_dist)
c                Fsteric_z= Fsteric*(unit_dispz/sub_dist)
c		write(10,*)unit_dispy, unit_dispz, Fsteric
c		write(10,*)unit_dispz/sub_dist
c		write(10,*)Fsteric_x,Fsteric_y,Fsteric_z
c                write(10,*)"reading this loop"
c		angle(1)=0; angle(2)=0; angle(3)=0
c		go to 500
c            end if
c
c            C_MT_dist= sqrt(copy_cargoy**2+copy_cargoz**2)                      ! when cargo volume gets cut by MT (steric spring force acting here)
c            C_MT_dist= sqrt(copy_cargoz**2)
c            if(C_MT_dist .lt. (beadRadius))then
c                cross_dist= abs(beadRadius-C_MT_dist)
c                unit_dispy=(copy_cargoy*cross_dist)/C_MT_dist
c                unit_dispz=(copy_cargoz*cross_dist)/C_MT_dist
c		dx= 0
c		dy= 0
c		dy= unit_dispy
c		dz= unit_dispz
c                span=110
c               angle(1)=0; angle(2)=0; angle(3)=0
c               go to 500
c            end if

            do j=1,m								! force on each motor along x-axis
c	      write(10,*)"Fonmot=", Fonmot(j)
              StepR(j)=0; DetR(j)=0; AttR(j)=0
	      ps(j)=0; pd(j)=0; pa(j)=0
            end do

c	    write(10,*) "rates of relevant motors="
	    do i=1,m								! calculation of step/attach/detach rates acc to force acting on motors
	      if(dup_stalk(i).ne.0)then
	       AttR(i)=0
	       if(Fonmot(i).le.0)then						! backward force
	        if(abs(Fonmot(i)).le.abs(Fs))then
	           stepcons= 1-(abs(Fonmot(i)/Fs)**w)
		   detachcons= abs(Fonmot(i)/Fd)
	  	   StepR(i)=Rso*stepcons
	           DetR(i)=Rdo*exp(detachcons)
	        else 
	          StepR(i)=0
                  DetR(i)= 1.07+ abs(0.186*Fonmot(i))
	        end if
	       else 									       	     !write(10,*)'backward motor'	
	        StepR(i)=Rso
	 	if(abs(Fonmot(i)).le.abs(Fs))then
		  detachcons= abs(Fonmot(i)/Fd)
	 	  DetR(i)= Rdo*exp(detachcons)
		else 
                  DetR(i)= 1.07+ abs(0.186*Fonmot(i))
	 	end if
	       end if
	      else if(dup_stalk(i)== 0)then
	       AttR(i)= Bindingrate
	       DetR(i)=0; StepR(i)=0
	      end if
c             write(10,*)"Rates=",i,AttR(i),DetR(i),StepR(i)
	    end do

c	    write(10,*)"probabilities="
	    do i=1,m
		ps(i)= StepR(i)*del_t
		pd(i)= DetR(i)*del_t
		pa(i)= AttR(i)*del_t
c            write(10,*)ps(i), pd(i), pa(i)
	    end do
	    
c	    write(10,*)"probabilities and decision="
	    do i=1,m
             r=rand()
	     if(r .le. ps(i))then
	      	dup_mothx(i)= rel_mothx(i) + d
	        if(time .ge. 0.5)step= step +1
c		step= step +1
                do k=1,m
                 if(dup_mothx(k)==dup_mothx(i) .and. (k.ne.i))then
                   dup_mothx(i)= dup_mothx(i)-d
c		   step= step -1
                   if(time .ge. 0.5)step= step - 1
                 end if
                end do
c               write(10,*)"stepped"
	     end if
	     if(r.gt.ps(i) .and. r .le.(ps(i)+pd(i)))then
		dup_mothx(i)= rel_motbx(i)
		dup_mothy(i)= rel_motby(i)
		dup_mothz(i)= rel_motbz(i)			
		dup_stalk(i)=0
c                detach= detach + 1
		if(time .ge. 0.5)detach= detach + 1
	     end if
	     if(r.gt.(ps(i)+pd(i)).and.r.le.(ps(i)+pd(i)+pa(i)))then
             if(rel_range(i)==1)then
              if(abs(dup_motby(i)) .gt. 0.012)then
                i_mot= int((dup_motbx(i)*1000)/8)
                icheck= i_mot*8
                icheck1= (i_mot+1)*8
                sort(1)= sqrt((abs(dup_motby(i))-0.012)**2+
     1          (dup_motbx(i)-(0.001*icheck))**2+dup_motbz(i)**2)
                x(1)= icheck; y(1)=12
                sort(2)= sqrt((abs(dup_motby(i))-0.012)**2+
     1          (abs(dup_motbx(i))-(0.001*icheck1))**2+dup_motbz(i)**2)
                x(2)= icheck1; y(2)=12
                 if(sort(2) .lt. sort(1))then
                    tempor= sort(2)
                    sort(2)= sort(1)
                    sort(1)= tempor
                    temporx= x(2)
                    x(2)=x(1)
                    x(1)=temporx
                 end if
                 if(sort(1) .le. Lo)then
                  if(dup_motby(i).lt.0)dup_mothy(i)=-(y(1)*0.001)
                  if(dup_motby(i).gt.0)dup_mothy(i)=(y(1)*0.001)
                  dup_mothx(i)= (x(1)*0.001); dup_mothz(i)=0
                  if(time .ge. 0.5)attach= attach + 1
c 		  attach= attach + 1
                 else
                  if(sort(2) .le. Lo)then
                   if(dup_motby(i).lt.0)dup_mothy(i)=-(y(2)*0.001)
                   if(dup_motby(i).gt.0)dup_mothy(i)=(y(2)*0.001)
                   dup_mothx(i)= (x(2)*0.001); dup_mothz(i)=0
                   if(time .ge. 0.5)attach= attach + 1
c		   attach= attach + 1
                  end if
                 end if
              else if(abs(dup_motby(i)) .lt. 0.012)then
                i_mot= int((dup_motbx(i)*1000)/8)
                i_moty= int((dup_motby(i)*1000)/4)
                icheck = i_mot*8
                icheck1= (i_mot+1)*8
                ichecky= i_moty*4
                ichecky1= (i_moty+1)*4
                sort(1)= sqrt((abs(dup_motby(i))-(0.001*ichecky))**2+
     1          (abs(dup_motbx(i))-(0.001*icheck))**2+dup_motbz(i)**2)
                x(1)=icheck; y(1)=ichecky
                sort(2)= sqrt((abs(dup_motby(i))-(0.001*ichecky1))**2+
     1          (abs(dup_motbx(i))-(0.001*icheck))**2+dup_motbz(i)**2)
                x(2)= icheck; y(2)= ichecky1
                sort(3)= sqrt((abs(dup_motby(i))-(0.001*ichecky))**2+
     1          (abs(dup_motbx(i))-(0.001*icheck1))**2+dup_motbz(i)**2)
                x(3)= icheck1; y(3)= ichecky
                sort(4)= sqrt((abs(dup_motby(i))-(0.001*ichecky1))**2+
     1          (abs(dup_motbx(i))-(0.001*icheck1))**2+dup_motbz(i)**2)
                x(4)= icheck1; y(4)= ichecky1
                do l=1,4
                 do j=1,l
                   if(sort(l) .lt. sort(j))then
                    tempor= sort(l)
                    sort(l)= sort(j)
                    sort(j)= tempor
                    temporx= x(l)
                    x(l)=x(j)
                    x(j)=temporx
                    tempory= y(l)
                    y(l)=y(j)
                    y(j)=tempory
                   end if
                 end do
                end do
                do j=1,4
                if(sort(j) .le. Lo)then
                  if(dup_motby(i).lt.0)dup_mothy(i)=-(y(j)*0.001)
                  if(dup_motby(i).gt.0)dup_mothy(i)=(y(j)*0.001)
                  dup_mothx(i)= (x(j)*0.001); dup_mothz(i)=0
                  if(time .ge. 0.5)attach= attach + 1
c		  attach= attach + 1
                  exit
                end if
                end do
              end if
             end if
            end if
            end do
            time= time +1e-5
	    itime= int(time*1e5)
            act_mot= act_mot +m
	
	    do i=1,m
	     do k=1,i
	     if(dup_mothx(i)==dup_mothx(k).and.(i.ne.k))then
		if(dup_mothy(i)==dup_mothy(k))then
		dup_mothx(i)=dup_motbx(i)
		dup_mothy(i)=dup_motby(i)
		dup_mothz(i)=dup_motbz(i)
		dup_stalk(i)=0
		end if
	     end if
	     end do
	    end do

	    do i=1,N
	     do k=1,m
		if(copy_beadmx(i)==dup_motbx(k).and.
     1		copy_beadmy(i)==dup_motby(k) .and.
     1		copy_beadmz(i)==dup_motbz(k))then
                  copy_headmx(i)= dup_mothx(k)
                  copy_headmy(i)= dup_mothy(k)
                  copy_headmz(i)= dup_mothz(k)
		end if
	     end do
                s1= (copy_beadmx(i)-copy_headmx(i))**2+
     1          (copy_beadmy(i)-copy_headmy(i))**2+
     1          (copy_beadmz(i)-copy_headmz(i))**2
                stalk_motl(i)= sqrt(s1)
c		print*, i, stalk_motl(i)
	    end do

	    icount=0	    
	    do i=1,N							! update the length and position of motor due to its step/detach/attachment
	      if(stalk_motl(i).ne.0)then
		icount= icount+1
	      end if
	    end do
c	    write(10,*)"no of attached motors=",icount
	
	   if(icount==0 .and. itime .lt. 1e5)then
		go to 300
	   end if
           if(itime== 1e5)then
                if(icount==0)then
                  go to 300
                end if
           end if
           if(itime== 0.5e5)then
                if(icount==0)then
                  go to 300
                else
                 bead_first = copy_cargox
                end if
           end if
           if(itime == 1e5)then
		bead_final= copy_cargox
		Rs= bead_final-bead_first
c                if(Rs .lt. 0 .and. abs(Rs) .gt. 1e5)then
c                  sigconf= sigconf-1
c                  go to 300
c                end if
c                if(Rs .gt. 1e6)then
c                  sigconf= sigconf-1
c                  go to 300
c                end if
		act_mot= ((act_mot/time)*1e-5)
		sigconf2= sigconf2 +1
		tk_time= time-0.5
                detach= detach/act_mot
                attach= attach/tk_time
                step= step/tk_time
		print*, "here",act_mot, detach, attach, step
		velocity= Rs/tk_time

                sum_act= sum_act + act_mot
                sum_t= sum_t + tk_time
                sum_det= sum_det + detach
                sum_att= sum_att + attach
                sum_step= sum_step + step
		vel_sum= vel_sum + velocity

                sumsq_mot= sumsq_mot + (act_mot**2)
                tsq_sum= tsq_sum + (tk_time**2)
                detsq_sum= detsq_sum + (detach**2)
                attsq_sum= attsq_sum + (attach**2)
                stepsq_sum= stepsq_sum + (step**2)
		vel_sqsum= vel_sqsum + (velocity**2)
                print*, "step=",step, "RL=",Rs,"time=",
     1          tk_time,"V=",velocity,"T=",tk_time
                go to 300                 
          end if

            do i=1,N
                F_mot(i)=0
                F_motx(i)=0;F_moty(i)=0
                F_motz(i)=0
            end do

c           write(10,*)"force by relevant motors="
            do i=1,m
             if (dup_stalk(i) .le. Lo)then                                            ! case of detached/unattached motor
               F_mot(i)=0
             else
               F_mot(i)= (dup_stalk(i)-Lo)*springk                                     ! magnitude of force exerted by attached and extended motors
             end if
c             if((dup_mothx(i)-dup_motbx(i)).lt.0)then
c                  Fonmot(i)= F_mot(i)
c             end if
c             if((dup_mothx(i)-dup_motbx(i)).ge.0)then
c                  Fonmot(i)= -F_mot(i)
c             end if
c            write(10,*)i, dup_mothx(i)-dup_motbx(i)
c            write(10,*) i, F_mot(i), Fonmot(i)
            end do

            do i=1,m
                del(i,1)=0; del(i,2)=0; del(i,3)=0
                tau(i,1)=0; tau(i,2)=0; tau(i,3)=0
                F_motx(i)=0; F_moty(i)=0; F_motz(i)=0
                u_motx(i)=0; u_moty(i)=0; u_motz(i)=0
            end do

c           write(10,*)"force and its component by motors="
            sumF_x=0; sumF_y=0; sumF_z=0
            sum_torq(1)=0; sum_torq(2)=0; sum_torq(3)=0

            do i=1,m
             if (F_mot(i).eq. 0)then                                            ! case of detached/unattached motor
                F_motx(i)=0; F_moty(i)=0; F_motz(i)=0
                del(i,1)=0; del(i,2)=0; del(i,3)=0
                u_motx(i)=0; u_moty(i)=0; u_motz(i)=0
                tau(i,1)=0; tau(i,2)=0; tau(i,3)=0
             else if(F_mot(i).ne.0)then
c               F_mot(i)= (dup_stalk(i)-Lo)*springk                            ! magnitude of force exerted by attached and extended motors
                del(i,1)= dup_motbx(i)- copy_cargox
                del(i,2)= dup_motby(i)- copy_cargoy
                del(i,3)= dup_motbz(i)- copy_cargoz
                u_motx(i)= (dup_mothx(i)-dup_motbx(i))/dup_stalk(i)                 ! unit vector of x,y,z line connecting head and tail of motors
                u_moty(i)= (dup_mothy(i)-dup_motby(i))/dup_stalk(i)
                u_motz(i)= (dup_mothz(i)-dup_motbz(i))/dup_stalk(i)
                F_motx(i)= F_mot(i)*u_motx(i)
                F_moty(i)= F_mot(i)*u_moty(i)
                F_motz(i)= F_mot(i)*u_motz(i)
                tau(i,1)= del(i,2)*F_motz(i)-del(i,3)*F_moty(i)
                tau(i,2)= -del(i,1)*F_motz(i)+del(i,3)*F_motx(i)
                tau(i,3)= del(i,1)*F_moty(i)-del(i,2)*F_motx(i)
             end if

             sumF_x= sumF_x + F_motx(i)
             sumF_y= sumF_y + F_moty(i)
             sumF_z= sumF_z + F_motz(i)
             sum_torq(1)= sum_torq(1) + tau(i,1)
             sum_torq(2)= sum_torq(2) + tau(i,2)
             sum_torq(3)= sum_torq(3) + tau(i,3)
            end do
c           write(10,*) sumF_x, sumF_y, sumF_z
c           write(10,*) sum_torq(1), sum_torq(2), sum_torq(3)

3           r1= 2*rand()-1
            r2= 2*rand()-1
            gauss= r1**2 + r2**2
            if(gauss.ge.1 .or. gauss.eq.0)go to 3
            var= r1*sqrt(-2*log(gauss)/gauss)
            x_random= var*sqrt(twokbt*vari2)

9           b1= 2*rand()-1
            b2= 2*rand()-1
            gem= b1**2 + b2**2
            if(gem.ge.1 .or. gem.eq.0)go to 9
            vr1= b1*sqrt(-2*log(gem)/gem)
            y_random= vr1*sqrt(twokbt*vari2)

8           c1= 2*rand()-1
            c2= 2*rand()-1
            gas= c1**2 + c2**2
            if(gas.ge.1 .or. gas.eq.0)go to 8
            var2= c1*sqrt(-2*log(gas)/gas)
            z_random= var2*sqrt(twokbt*vari2)
c           write(10,*)"random force="
c           write(10,*) Fran_x, Fran_y, Fran_z

2           d1= 2*rand()-1
            d2= 2*rand()-1
            fim= d1**2+d2**2
            if(fim.ge.1 .or. fim.eq.0)go to 2
            var3= d1*sqrt(-2*log(fim)/fim)
            tau_ran(1)= var3*sqrt(twokbt*variate)

6           e1= 2*rand()-1
            e2= 2*rand()-1
            tim= e1**2+ e2**2
            if(tim.ge.1 .or. tim.eq.0)go to 6
            var4= e1*sqrt(-2*log(tim)/tim)
            tau_ran(2)= var4*sqrt(twokbt*variate)

4           f1= 2*rand()-1
            f2= 2*rand()-1
            bim= f1**2+f2**2
            if(bim.ge.1 .or. bim.eq.0)go to 4
            var5= f1*sqrt(-2*log(bim)/bim)
            tau_ran(3)= var5*sqrt(twokbt*variate)

c            C_MT_dist= sqrt(copy_cargoy**2+copy_cargoz**2)                      ! when cargo volume gets cut by MT (steric spring force acting here)
            C_MT_dist= sqrt(copy_cargoz**2)
            if(C_MT_dist .lt. (beadRadius))then
                cross_dist= abs(beadRadius-C_MT_dist)
c                unit_dispy=(copy_cargoy*cross_dist)/C_MT_dist
                unit_dispz=(copy_cargoz*cross_dist)/C_MT_dist
                dx= 0
                dy= 0
c               dy= unit_dispy
                dz= unit_dispz
c                span=110
               angle(1)=0; angle(2)=0; angle(3)=0
               go to 500
            end if

c           x_random=0; y_random=0; z_random=0
c          tau_ran(1)=0; tau_ran(2)=0; tau_ran(3)=0
c          sum_torq(1)=0; sum_torq(2)=0; sum_torq(3)=0
          do i=1,3
                angle(i)= (sum_torq(i)*variate)+ tau_ran(i)
c               write(10,*) "angle=",angle(i)
           end do
           dx= ((sumF_x + Fx_ext)*vari2) + x_random
           dy= (sumF_y*vari2) + y_random
           dz= (sumF_z*vari2) + z_random
c            write(10,*) "dx=", dx, dy, dz

            if(time .le. 4e-5)then
                k=int(time*1e5)
                store_dx(k)= dx
            else if(time .gt. 4e-5)then
                store_dx(1)=store_dx(2)
                store_dx(2)=store_dx(3)
                store_dx(3)=store_dx(4)
                store_dx(4)=dx
            end if
            st_dx_4=abs(store_dx(4))-abs(store_dx(3))
            st_dx_3=abs(store_dx(3))-abs(store_dx(2))
            st_dx_2=abs(store_dx(2))-abs(store_dx(1))

            if(st_dx_4 .gt. 5 .and. st_dx_3 .gt. 5)then
             if(st_dx_2 .gt. 5)then
                mcount= mcount +1
c               print*, "arrested case"
c                write(10,*) "arrested case"
             end if
            end if
            icount=0
            do i=1,m
                if(dup_stalk(i) . gt. 0)icount= icount +1
            end do

            if(icount .gt. 1 .and. mcount .ge. 4)then
             do i=2,m
                if(dup_stalk(i) .gt. 0)then
                do j=1,i-1
                  if(dup_mothx(i) .lt. dup_mothx(j) .and.
     1                  dup_stalk(j) .gt. 0)then
                   xdown = dup_mothx(i)
                   ydown= dup_mothy(i)
                  end if
                end do
                end if
             end do
             do i=2,m
                if(dup_stalk(i) .gt. 0)then
                do j=1,i-1
                  if(dup_mothx(i) .gt. dup_mothx(j) .and.
     1                  dup_stalk(j) .gt. 0)then
                   xup = dup_mothx(i)
                   yup= dup_mothy(i)
                  end if
                end do
                end if
             end do
             dx= ((xup+xdown)/2)-copy_cargox
             if(yup .gt. ydown)then
                dy= ((yup+ydown)/2)-copy_cargoy
             else
                dy= ((ydown+yup)/2)-copy_cargoy
             end if
             dz= beadRadius-copy_cargoz
             angle(1)=0; angle(2)=0; angle(3)=0
           end if

c500         t11= cos(angle(2))*cos(angle(3))                                            !clockwise
c            t12= cos(angle(2))*sin(angle(3)); t13= -sin(angle(2))
c            t21= sin(angle(1))*sin(angle(2))*cos(angle(3))-
c     1          sin(angle(3))*cos(angle(1))
c            t22= sin(angle(1))*sin(angle(2))*sin(angle(3))+
c     1          cos(angle(3))*cos(angle(1))
c            t23= sin(angle(1))*cos(angle(2))
c            t31= cos(angle(1))*sin(angle(2))*cos(angle(3))+
c     1          sin(angle(3))*sin(angle(1))
c            t32= cos(angle(1))*sin(angle(2))*sin(angle(3))-
c     1          cos(angle(3))*sin(angle(1))
c            t33= cos(angle(1))*cos(angle(2))

500         t11= cos(angle(2))*cos(angle(3))                                            !anti-clockwise
            t12= -cos(angle(2))*sin(angle(3)); t13= sin(angle(2))
            t21= sin(angle(1))*sin(angle(2))*cos(angle(3))+
     1          sin(angle(3))*cos(angle(1))
            t22= -sin(angle(1))*sin(angle(2))*sin(angle(3))+
     1          cos(angle(3))*cos(angle(1))
            t23= -sin(angle(1))*cos(angle(2))
            t31= -cos(angle(1))*sin(angle(2))*cos(angle(3))+
     1          sin(angle(3))*sin(angle(1))
            t32= cos(angle(1))*sin(angle(2))*sin(angle(3))+
     1          cos(angle(3))*sin(angle(1))
            t33= cos(angle(1))*cos(angle(2))

	   do i=1,N
		anchorx(i)=0; anchory(i)=0; anchorz(i)=0
	   end do

	    temp1=0; temp2=0; temp3=0
            do i=1,N                                                                      !ROTATION OF ATTACHMENT POINTS ON BEAD SURFACE ABOUT X-, Y- AND Z- AXES
             temp1= t11*(copy_beadmx(i)-copy_cargox)+
     1       t12*(copy_beadmy(i)-copy_cargoy)+
     1		t13*(copy_beadmz(i)-copy_cargoz)
             anchorx(i)= temp1+ copy_cargox

             temp2= t21*(copy_beadmx(i)-copy_cargox)+
     1       t22*(copy_beadmy(i)-copy_cargoy)+
     1		t23*(copy_beadmz(i)-copy_cargoz)
             anchory(i)= temp2+ copy_cargoy

             temp3= t31*(copy_beadmx(i)-copy_cargox)+
     1       t32*(copy_beadmy(i)-copy_cargoy)+
     1		t33*(copy_beadmz(i)-copy_cargoz)
             anchorz(i)= temp3+ copy_cargoz
	    end do

c           if(span == 110)then
	    final_cargox= copy_cargox + dx
	    final_cargoy= copy_cargoy + dy
	    final_cargoz= copy_cargoz + dz
c           write(10,*)"f_C=",final_cargox,
c     1         final_cargoy,final_cargoz

	    do i=1,N
		copy_beadmx(i)= dx + anchorx(i)   						!+ copy_beadmx(i)
		copy_beadmy(i)= dy + anchory(i)							!+ copy_beadmy(i)
		copy_beadmz(i)= dz + anchorz(i)  						!+ copy_beadmz(i)
             pseudo_rad=sqrt((copy_beadmx(i)-final_cargox)**2+
     1          (copy_beadmy(i)-final_cargoy)**2+
     1          (copy_beadmz(i)-final_cargoz)**2)
	     delta= pseudo_rad-beadRadius
	     if(abs(delta).gt.1e-6)then
		ux= (copy_beadmx(i)-final_cargox)/beadRadius
		uy= (copy_beadmy(i)-final_cargoy)/beadRadius
		uz= (copy_beadmz(i)-final_cargoz)/beadRadius
		copy_beadmx(i)= copy_beadmx(i)- (ux*delta)
		copy_beadmy(i)= copy_beadmy(i)- (uy*delta)
		copy_beadmz(i)= copy_beadmz(i)- (uz*delta)
	     end if
             if(stalk_motl(i)==0)then
                copy_headmx(i)= copy_beadmx(i)
                copy_headmy(i)= copy_beadmy(i)
                copy_headmz(i)= copy_beadmz(i)
             end if
c             write(10,*)copy_beadmx(i),copy_beadmy(i),copy_beadmz(i)
c             write(10,*)copy_headmx(i),copy_headmy(i),copy_headmz(i)
	    end do
c            print*, final_cargox, final_cargoy,final_cargoz
c           else
c            final_cargox= xdown + dx
c            final_cargoy= ydown + dy
c            final_cargoz= beadRadius
c            write(*,*)copy_beadmx(i),copy_beadmy(i),copy_beadmz(i)
c            write(*,*)copy_headmx(i),copy_headmy(i),copy_headmz(i)
c           print*, final_cargox, final_cargoy,final_cargoz                                                     !copy_cargox
c           end if

           do i=1,N
           bead_MTdist(i)= sqrt((copy_beadmz(i)**2))
c           print*, "dist=", bead_MTdist(i)
           if(bead_MTdist(i) .le. Lo)then
            if(abs(copy_beadmy(i)) .gt. 0.012)then
                i_mot= int((copy_beadmx(i)*1000)/8)
                icheck= i_mot*8
                icheck1= (i_mot+1)*8
                xcheck= sqrt((abs(copy_beadmy(i))-0.012)**2+
     1          (copy_beadmx(i)-(0.001*icheck))**2+copy_beadmz(i)**2)
                xcheck1= sqrt((abs(copy_beadmy(i))-0.012)**2+
     1          (copy_beadmx(i)-(0.001*icheck1))**2+copy_beadmz(i)**2)
                if(xcheck .le. Lo .or. xcheck1 .le. Lo)then
                  range(i)=1
                 else
                  range(i)=0
                end if
            else if(abs(copy_beadmy(i)) .lt. 0.012)then
                i_mot= int((copy_beadmx(i)*1000)/8)
                i_moty= int((copy_beadmy(i)*1000)/4)
                icheck = i_mot*8
                icheck1= (i_mot+1)*8
                ichecky= i_moty*4
                ichecky1= (i_moty+1)*4
                sort(1)= sqrt((copy_beadmy(i)-(0.001*ichecky))**2+
     1          (copy_beadmx(i)-(0.001*icheck))**2+copy_beadmz(i)**2)
                x(1)=icheck; y(1)=ichecky
                sort(2)= sqrt((copy_beadmy(i)-(0.001*ichecky1))**2+
     1          (copy_beadmx(i)-(0.001*icheck))**2+copy_beadmz(i)**2)
                x(2)= icheck; y(2)= ichecky1
                sort(3)= sqrt((copy_beadmy(i)-(0.001*ichecky))**2+
     1          (copy_beadmx(i)-(0.001*icheck1))**2+copy_beadmz(i)**2)
                x(3)= icheck1; y(3)= ichecky
                sort(4)= sqrt((copy_beadmy(i)-(0.001*ichecky1))**2+
     1          (copy_beadmx(i)-(0.001*icheck1))**2+copy_beadmz(i)**2)
                x(4)= icheck1; y(4)= ichecky1
                if(sort(1) .le. Lo .or. sort(2) .le. Lo .or.
     1          sort(3) .le. Lo .or. sort(4) .le. Lo) then
                 range(i)=1
                else
                 range(i)=0
                end if
             end if
            else
               range(i)=0
            end if
            s1=((copy_beadmx(i)-copy_headmx(i))**2+
     1         (copy_beadmy(i)-copy_headmy(i))**2+
     1         (copy_beadmz(i)-copy_headmz(i))**2)
             stalk_motl(i)= sqrt(s1)
           end do

c            do i=1,N
c             bead_MTdist(i)=sqrt((copy_beadmy(i)**2)+
c     1          (copy_beadmz(i)**2))
c             if(bead_MTdist(i) .le. Lo)then
c              icheck= int(copy_beadmx(i)/8)
c              x_mot= icheck*8
c              x_mot1= (icheck+1)*8
c              xcheck=sqrt((copy_beadmx(i)-x_mot)**2+
c     1          copy_beadmy(i)**2+copy_beadmz(i)**2)
c              xcheck1=sqrt((copy_beadmx(i)-x_mot1)**2+
c     1          copy_beadmy(i)**2+ copy_beadmz(i)**2)
c              if(xcheck .le. Lo .or. xcheck1 .le. Lo)then
c                range(i)=1
c              else
c                range(i)=0
c              end if
c             else
c                range(i)=0
c             end if
c             s1=(copy_beadmx(i)-copy_headmx(i))**2+
c     1         (copy_beadmy(i)-copy_headmy(i))**2+
c     1         (copy_beadmz(i)-copy_headmz(i))**2
c             stalk_motl(i)= sqrt(s1)
c            write(10,*)i ,stalk_motl(i),range(i)
c            end do

	    k=0
	    do i=1,N
		if(range(i).eq.1 .or. stalk_motl(i).gt.0)then
		k=k+1
		dup_motbx(k)=copy_beadmx(i)
		dup_motby(k)=copy_beadmy(i)
		dup_motbz(k)=copy_beadmz(i)
		dup_mothx(k)=copy_headmx(i)
		dup_mothy(k)=copy_headmy(i)
		dup_mothz(k)=copy_headmz(i)
		rel_range(k)=range(i)
		dup_stalk(k)=stalk_motl(i)
		end if
c	     write(10,*)copy_beadmx(i),copy_beadmy(i),copy_beadmz(i)
c             write(10,*)copy_headmx(i),copy_headmy(i),copy_headmz(i)
	    end do
	    m=k	
		
	    do i=1,m
	      rel_motbx(i)= dup_motbx(i)
	      rel_motby(i)= dup_motby(i)
	      rel_motbz(i)= dup_motbz(i)
	      rel_mothx(i)= dup_mothx(i)
	      rel_mothy(i)= dup_mothy(i)
	      rel_mothz(i)= dup_mothz(i)
c	     write(10,*)dup_motbx(i),dup_motby(i),dup_motbz(i)
c            write(10,*)dup_mothx(i),dup_mothy(i),dup_mothz(i)
c	     write(10,*)dup_stalk(i)
	    end do

c            do i=1,N
c             if(stalk_motl(i) .ge. 1e7)then
c                sigconf=sigconf-1
c                go to 300
c             end if
c            end do

c            write(10,*)"anchor points and cargo center="
            do i=1,N                                                            ! updating original array of motor head and tail (bead) positions
	     pseudo_rad=0
             beadmotx(i)= copy_beadmx(i)
             beadmoty(i)= copy_beadmy(i)
             beadmotz(i)= copy_beadmz(i)
             headmotx(i)= copy_headmx(i)
             headmoty(i)= copy_headmy(i)
	     headmotz(i)= copy_headmz(i)
c            write(10,*)beadmotx(i),beadmoty(i),beadmotz(i)
c            write(10,*)headmotx(i),headmoty(i),headmotz(i)
	    end do
            copy_cargox= final_cargox
            copy_cargoy= final_cargoy
            copy_cargoz= final_cargoz

	  end do
300       continue
c	  print*, "significant conf=", sigconf
	  if (sigconf2==maxconf)exit  
	 end do
         avgact= sum_act/sigconf2
         avgt= sum_t/sigconf2
         detavg= sum_det/(avgt*sigconf2)
         attavg= sum_att/sigconf2
         stpavg= sum_step/sigconf2

         Actsq_avg= sumsq_mot/sigconf2
         active_SD= sqrt(Actsq_avg-(avgact**2))
         active_SEM= active_SD/sqrt(sigconf2)

c         tavg_sq= tsq_sum/sigconf2
c         tavg_SD= sqrt(tavg_sq-(avgt**2))
c         tavg_SEM= tavg_SD/sqrt(sigconf2)

         detsq_avg= detsq_sum/(avgt*sigconf2)
         det_SD= sqrt(detsq_avg-(detavg**2))
         det_SEM= det_SD/sqrt(sigconf2)

         attsq_avg= attsq_sum/sigconf2
         att_SD= sqrt(attsq_avg-(attavg**2))
         att_SEM= att_SD/sqrt(sigconf2)

         stpsq_avg= stepsq_sum/sigconf2
         step_SD= sqrt(stpsq_avg-(stpavg**2))
         step_SEM= step_SD/sqrt(sigconf2)

         avgvel= vel_sum/sigconf2
         vel_sqavg= vel_sqsum/sigconf2
         vel_SD= sqrt(vel_sqavg- (avgvel**2))
         vel_SEM= vel_SD/sqrt(sigconf2)

         print*, iF, "detach= ",detavg," attach= ",attavg," step= ",
     1          stpavg, sigconf," active motors=", avgact,
     1		"time=", avgt,"Vel=", avgvel

         write(21,*) iF,N,imtd,beadRadius,detavg,det_SD,det_SEM,
     1          attavg,att_SD,att_SEM,stpavg,step_SD,step_SEM

         write(23,*)iF,N,imtd,beadRadius,avgvel,vel_SD,vel_SEM,
     1		avgact,active_SD,active_SEM
c
c         write(24,*) iF, N, imtd,beadRadius,detsq_avg,attsq_avg,
c     1          stpsq_avg,Xavgsq,tavg_sq,vel_sqavg,Actsq_avg

	end do
	end do
c	close(10)
	close(21)
	close(23)
c	close(24)
	stop
	end
