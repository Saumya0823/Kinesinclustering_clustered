	Program cargorotationunderexternalforce
	parameter (N=5)
	real Lo
	dimension tau_ran(3),angle(3),copy_beadmx(N),
     1		copy_beadmz(N),copy_headmx(N),copy_headmy(N),
     1		copy_headmz(N),copy_beadmy(N),sort(10),x(10),y(10)
	dimension beadmotx(N),beadmoty(N),beadmotz(N),store_dx(4),
     1		headmotx(N),headmoty(N),headmotz(N),stalk_motl(200),
     1          dup_motbx(N),dup_motby(N),dup_motbz(N),
     1          dup_mothx(N),dup_mothy(N),dup_mothz(N),
     1          rel_motbx(N),rel_motby(N),rel_motbz(N),
     1          rel_mothx(N),rel_mothy(N),rel_mothz(N),
     1          rel_range(N),StepR(N),DetR(N),AttR(N),
     1		dup_stalk(N),sum_torq(3),del(200,200),
     1		ang_mot(N,3),Fall_mot(N),Fall_motx(N),
     1		Fall_moty(N),Fall_motz(N),taual1(N),taual2(N),
     1		taual3(N),sum_motin(1000),iterm(1000),agl(N),
     1          avg_area(1000),sm_area(1000),poly_x(N),poly_y(N)
	dimension thetamot(N),phimot(N),range(N),avg_motin(1000),
     1		F_motx(N),F_moty(N),F_motz(N),F_mot(N),
     1		Fonmot(N),u_motx(N),u_moty(N),u_motz(N),
     1		tau(200,200),ps(N),pd(N),pa(N),t_array(N),
     1		anchorx(N),anchory(N),anchorz(N),bead_MTdist(N)
        parameter (Fs=5.,iTemp=300,springk= 320,Bindingrate=5,w=2,
     1          d=0.008,Vo=1.0,Rso=Vo/d,Rdo=1,Fd=3.,Lo=0.110,
     1          boltzmann=1.3806e-5,pi=4*atan(1.),
     1          beadRadius=0.250,Rcube= beadRadius**3,
     1          twokbt= 2*boltzmann*iTemp,viscosity=1e-3,
     1          gamnew= 6*pi*viscosity*beadRadius,v= 180/pi,
     1          gamma= 8*pi*viscosity*Rcube, del_t=1e-5, maxconf=1000)

c        open(unit=10,file="Torque_clusmot_N5_D001.dat")
c        open(unit=11,file="angle_clusmot_N5_D001.dat")
c        open(unit=22,file="Rate_clusmot_N5_D001.dat")
c        open(unit=23,file="RLV_clusmot_N5_D001.dat")
        open(unit=20, file="motor_tails_clusmot_N5_D001_F0.dat")
        open(unit=21, file="time_points_clusmot_N5_D001_F0.dat")

c        SA= 4*pi*(beadRadius**2)
c        avgden= 9.7444*1e-5
c        N= SA*avgden
c	N= 1
c	iZ= N/2									!Total number of motors
c        imtd= 30
        variate= del_t/gamma
        vari2= del_t/gamnew
	vmot= (0.01/beadRadius)							!motor diffusion constant
	vgamm= (twokbt/vmot)							!motor gamma
	cn= del_t/vgamm
        
	do imtd=0,0						!60,10

	average_m=0
        do iF= 0,0								! force loop starts
          F_ext= -float(iF)
          Fz_ext= 0
          Fy_ext= 0
          Fx_ext= F_ext
	  print*, Fx_ext, Fy_ext, Fz_ext

        sigconf=0; xsum=0; sum_act=0; sum_t=0
        sum_det=0; sum_att=0; sum_step=0
        xsqsum=0; sumsq_mot=0; tsq_sum=0
        detsq_sum=0; attsq_sum=0; stepsq_sum=0
        vel_sum=0; vel_sqsum=0
        torq_avg=0; avg_angle=0
        torq_avg1=0; avg_angle1=0
        torq_avg2=0; avg_angle2=0
	Favg=0; Fxavg=0; Fyavg=0; Fzavg=0
        t_mot_avg=0
        print*, "conf loop starts"
c        do i=1,1000
c         avg_motin(i)=0; sum_motin(i)=0; iterm(i)=0
c         avg_area(i)=0; sm_area(i)=0
c        end do

        do                        		                                 ! conf loop starts
	  print*, "conf=", sigconf, beadRadius,'imtd=',imtd,N,iF
	  beadinitx= 0
	  beadinity= 0
	  beadinitz= beadRadius + imtd

	  do i=1,N
	   beadmotx(i)=0; beadmoty(i)=0; beadmotz(i)=0
	   headmotx(i)=0; headmoty(i)=0; headmotz(i)=0
	   thetamot(i)=0; phimot(i)=0
	   stalk_motl(i)=0
	  end do

c         do i=1,N
c1          phimot(i)= (pi/18)*rand()+ ((16*pi)/18)
c           thetamot(i)=  (pi/18)*rand()+((7*pi)/18)
c1          phimot(i)= (pi/6)*rand()+ ((4*pi)/6)
c           thetamot(i)=  (pi/6)*rand()+((7*pi)/6)
c           caviar1= beadRadius*cos(thetamot(i))*sin(phimot(i))
c           beadmotx(i)= beadinitx + caviar1
c           caviar2= beadRadius*sin(thetamot(i))*sin(phimot(i))
c          beadmoty(i)= beadinity + caviar2
c           beadmotz(i)= beadinitz + beadRadius*cos(phimot(i))
c           term= sqrt((beadmotx(i)-beadinitx)**2+
c     1          (beadmoty(i)-beadinity)**2+
c     1          (beadmotz(i)-beadinitz)**2)
c           term_diff= abs(term-beadRadius)
c           if(term_diff .gt. 1e-4)then
c                go to 1
c           end if
c           write(*,*)beadmotx(i),beadmoty(i),beadmotz(i)
c           write(11,*)beadmotx(i),beadmoty(i),beadmotz(i)
c        end do
c
c        do i=iZ+1,N
c600        phimot(i)= (pi/18)*rand()+ ((16*pi)/18)
c           thetamot(i)=  (pi/18)*rand()+((16*pi)/18)
c600        phimot(i)= (pi/6)*rand()+ ((4*pi)/6)
c           thetamot(i)=  (pi/6)*rand()+((10*pi)/6)
c           caviar1= beadRadius*cos(thetamot(i))*sin(phimot(i))
c           beadmotx(i)= beadinitx + caviar1
c           caviar2= beadRadius*sin(thetamot(i))*sin(phimot(i))
c           beadmoty(i)= beadinity + caviar2
c           beadmotz(i)= beadinitz + beadRadius*cos(phimot(i))
c           term= sqrt((beadmotx(i)-beadinitx)**2+
c     1          (beadmoty(i)-beadinity)**2+
c     1          (beadmotz(i)-beadinitz)**2)
c           term_diff= abs(term-beadRadius)
c           if(term_diff .gt. 1e-4)then
c               go to 600
c           end if
c        end do
c
        do i=1,N
1         phimot(i)= acos(2*rand()-1)
          thetamot(i)= twopi*rand()
          caviar1= beadRadius*cos(thetamot(i))*sin(phimot(i))
          beadmotx(i)= beadinitx + caviar1
          caviar2= beadRadius*sin(thetamot(i))*sin(phimot(i))
          beadmoty(i)= beadinity + caviar2
          beadmotz(i)= beadinitz + beadRadius*cos(phimot(i))
          term= sqrt((beadmotx(i)-beadinitx)**2+
     1         (beadmoty(i)-beadinity)**2+
     1         (beadmotz(i)-beadinitz)**2)
          term_diff= abs(term-beadRadius)*0.001
          if(term_diff .gt. 1e-4)then
               go to 1
          end if
          b_dis= 1000*sqrt((beadmoty(i)**2)+(beadmotz(i)**2))
          if(b_dis .gt. 6.125)go to 1
c          write(*,*)beadmotx(i),beadmoty(i),beadmotz(i)
         end do

	  do i=1,N
	   range(i)=0
	  end do

           do i=1,N
           bead_MTdist(i)= sqrt((beadmotz(i)**2))
c           print*, "dist=", bead_MTdist(i)
           if(bead_MTdist(i) .le. Lo)then
	     range(i)=1
	    if(abs(beadmoty(i)) .gt. 0.012)then
		i_mot= int((beadmotx(i)*1000)/8)
		icheck= i_mot*8
		icheck1= (i_mot+1)*8
		xcheck= sqrt((abs(beadmoty(i))-0.012)**2+
     1		(beadmotx(i)-(0.001*icheck))**2+beadmotz(i)**2)
                xcheck1= sqrt((abs(beadmoty(i))-0.012)**2+
     1          (beadmotx(i)-(0.001*icheck1))**2+beadmotz(i)**2)
	        if(xcheck .lt. xcheck1)then
		 if(xcheck .le. Lo)then
		  if(beadmoty(i) .lt. 0)headmoty(i)=-0.012
		  if(beadmoty(i) .gt. 0)headmoty(i)=0.012
		  headmotx(i)= icheck*0.001; headmotz(i)=0
		 else
                  headmotx(i)=beadmotx(i)
                  headmoty(i)=beadmoty(i)
                  headmotz(i)=beadmotz(i)
		 end if
		else
                 if(xcheck1 .le. Lo)then
                  if(beadmoty(i) .lt. 0)headmoty(i)=-0.012
                  if(beadmoty(i) .gt. 0)headmoty(i)=0.012
                  headmotx(i)= icheck1*0.001; headmotz(i)=0
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
                sort(1)= sqrt((beadmoty(i)-(ichecky*0.001))**2+
     1          (beadmotx(i)-(icheck*0.001))**2+beadmotz(i)**2)
		x(1)=icheck; y(1)=ichecky
                sort(2)= sqrt((beadmoty(i)-(ichecky1*0.001))**2+
     1          (beadmotx(i)-(icheck*0.001))**2+beadmotz(i)**2)
		x(2)= icheck; y(2)= ichecky1
                sort(3)= sqrt((beadmoty(i)-(ichecky*0.001))**2+
     1          (beadmotx(i)-(icheck1*0.001))**2+beadmotz(i)**2)
		x(3)= icheck1; y(3)= ichecky
                sort(4)= sqrt((beadmoty(i)-(ichecky1*0.001))**2+
     1          (beadmotx(i)-(icheck1*0.001))**2+beadmotz(i)**2)
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
             headmotx(i)=beadmotx(i)
             headmoty(i)=beadmoty(i)
             headmotz(i)=beadmotz(i)
           end if
           end do
           s1= (headmotx(i)-beadmotx(i))**2+
     1		(headmoty(i)-beadmoty(i))**2+
     1          (headmotz(i)-beadmotz(i))**2
           stalk_motl(i)= sqrt(s1)
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
	     end if
	  end do
	  m=k

	  do i=1,m
	   dup_motbx(i)= rel_motbx(i)
	   dup_motby(i)= rel_motby(i)
	   dup_motbz(i)= rel_motbz(i)
	   dup_mothx(i)= rel_mothx(i)
	   dup_mothy(i)= rel_mothy(i)
	   dup_mothz(i)= rel_mothz(i)
	  end do

          iv=0
          do i=1,N
            if(stalk_motl(i).ne.0)iv=iv+1
          end do
          average_m= average_m+ iv
          print*, "no of attached motors=", iv
          print*, "no of relevant motors=", m

          if(iv==0)then
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
	   t_array(i)=0
	  end do
          do i=1,4
            store_dx(i)=0
          end do

	  print*, "time loop starts"
	  Rs=0; torq_sum=0; angle_sum=0
	  angle_sum1=0; angle_sum2=0
	  torq_sum1=0; torq_sum2=0
	  attach=0; detach=0; step=0
	  time=1e-5; act_mot=m; Fsum=0
	  Fxsum=0; Fysum=0; Fzsum=0
	  mcount=0; itcount=0
          counterprev=0
          do 									! time loop starts
	    copy_cargox= final_cargox
	    copy_cargoy= final_cargoy
	    copy_cargoz= final_cargoz
c	    print*, time, copy_cargox

	    do j=1,m
	     dup_motbx(j)= rel_motbx(j)
	     dup_motby(j)= rel_motby(j)
	     dup_motbz(j)= rel_motbz(j)
	     dup_mothx(j)= rel_mothx(j)
	     dup_mothy(j)= rel_mothy(j)
	     dup_mothz(j)= rel_mothz(j)
	    end do

            do i=1,N
             copy_beadmx(i)= beadmotx(i)
             copy_beadmy(i)= beadmoty(i)
             copy_beadmz(i)= beadmotz(i)
             copy_headmx(i)= headmotx(i)
             copy_headmy(i)= headmoty(i)
             copy_headmz(i)= headmotz(i)
            end do
            
            do i=1,m
                F_mot(i)=0; Fonmot(i)=0
                F_motx(i)=0;F_moty(i)=0
                F_motz(i)=0
            end do

            do i=1,m
             if (dup_stalk(i) .le. Lo)then 					      ! case of detached/unattached motor
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
c             Fsum= Fsum + Fonmot(i) 
            end do

            do j=1,m								! force on each motor along x-axis
              StepR(j)=0; DetR(j)=0; AttR(j)=0
              ps(j)=0; pd(j)=0; pa(j)=0
            end do

	    do i=1,m								! calculation of step/attach/detach rates acc to force acting on motors
	      if(dup_stalk(i).ne.0)then
	       AttR(i)=0
	       if(Fonmot(i).le.0)then						! backward force
	        if(abs(Fonmot(i)).lt.abs(Fs))then
	           stepcons= 1-abs((Fonmot(i)/Fs)**w)
		   detachcons= abs(Fonmot(i)/Fd)
	  	   StepR(i)=Rso*stepcons
c		   StepR(i)= Rso
	           DetR(i)=Rdo*exp(detachcons)
	        else
	          StepR(i)=0
		  DetR(i)= 1.07+ abs(0.186*Fonmot(i))  
	        end if
	       else 							! forward force
	        StepR(i)=Rso
c		DetR(i)=0
	 	if(abs(Fonmot(i)).lt.abs(Fs))then
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
	    end do

	    do i=1,m
		ps(i)= StepR(i)*del_t
		pd(i)= DetR(i)*del_t
		pa(i)= AttR(i)*del_t
	    end do

c           C_MT_dist= sqrt(copy_cargoy**2+copy_cargoz**2)                      ! when cargo volume gets cut by MT (steric spring force acting here)
            C_MT_dist= sqrt(copy_cargoz**2)
c           write(10,*) "CMTD=", C_MT_dist
            if(C_MT_dist .lt. beadRadius)then
                cross_dist= abs(beadRadius-C_MT_dist)
c                unit_dispy=(copy_cargoy*cross_dist)/C_MT_dist
                unit_dispz=(copy_cargoz*cross_dist)/C_MT_dist
                dx= 0
                dy=0
c               dy= unit_dispy
                dz= unit_dispz
               angle(1)=0; angle(2)=0; angle(3)=0
               go to 500
            end if
	    
	    do i=1,m
             r=rand()
	     if(r .le. ps(i))then
	      	dup_mothx(i)= rel_mothx(i) + d
		step= step + 1
		do k=1,m
		 if(dup_mothx(k)==dup_mothx(i) .and. (k.ne.i))then
		   dup_mothx(i)= dup_mothx(i)-d
		   step= step - 1
		 end if
		end do
	     end if
	     if(r.gt.ps(i) .and. r .le.(ps(i)+pd(i)))then
		dup_mothx(i)= rel_motbx(i)
		dup_mothy(i)= rel_motby(i)
		dup_mothz(i)= rel_motbz(i)	
		detach= detach + 1
		dup_stalk(i)=0
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
		 else
		  if(sort(2) .le. Lo)then
                   if(dup_motby(i).lt.0)dup_mothy(i)=-(y(2)*0.001)
                   if(dup_motby(i).gt.0)dup_mothy(i)=(y(2)*0.001)
                   dup_mothx(i)= (x(2)*0.001); dup_mothz(i)=0
		  end if
		 end if
              else if(abs(dup_motby(i)) .lt. 0.012)then
                i_mot= int((dup_motbx(i)*1000)/8)
                i_moty= int((dup_motby(i)*1000)/4)
                icheck = i_mot*8
                icheck1= (i_mot+1)*8
                ichecky= i_moty*4
                ichecky1= (i_moty+1)*4
                sort(1)= sqrt((abs(dup_motby(i))-(ichecky*0.001))**2+
     1          (abs(dup_motbx(i))-(icheck*0.001))**2+dup_motbz(i)**2)
                x(1)=icheck; y(1)=ichecky
                sort(2)= sqrt((abs(dup_motby(i))-(ichecky1*0.001))**2+
     1          (abs(dup_motbx(i))-(icheck*0.001))**2+dup_motbz(i)**2)
                x(2)= icheck; y(2)= ichecky1
                sort(3)= sqrt((abs(dup_motby(i))-(ichecky*0.001))**2+
     1          (abs(dup_motbx(i))-(icheck1*0.001))**2+dup_motbz(i)**2)
                x(3)= icheck1; y(3)= ichecky
                sort(4)= sqrt((abs(dup_motby(i))-(ichecky1*0.001))**2+
     1          (abs(dup_motbx(i))-(icheck1*0.001))**2+dup_motbz(i)**2)
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
		  exit
		end if
                end do
              end if
             end if
	    end if
	    end do
            time= time +1e-5
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
     1          copy_beadmy(i)==dup_motby(k) .and.
     1          copy_beadmz(i)==dup_motbz(k))then
                  copy_headmx(i)= dup_mothx(k)
                  copy_headmy(i)= dup_mothy(k)
                  copy_headmz(i)= dup_mothz(k)
                end if
             end do
             s1= (copy_beadmx(i)-copy_headmx(i))**2+
     1          (copy_beadmy(i)-copy_headmy(i))**2+
     1          (copy_beadmz(i)-copy_headmz(i))**2
             stalk_motl(i)= sqrt(s1)
            end do

            do i=1,N
             if(stalk_motl(i) .gt. 0)then
               t_array(i)= t_array(i) + 1e-5
             end if
            end do

	    icount=0	    
	    do i=1,N								! update the length and position of motor due to its step/detach/attachment
	      if(stalk_motl(i).ne.0)then
		icount= icount+1
	      end if
	    end do
	
c            print*, time, copy_cargox
	  if(icount==0)then
		Rs= copy_cargox
		torq_avg= torq_avg + ((torq_sum/time)*1e-5)
                torq_avg1= torq_avg1 + ((torq_sum1/time)*1e-5)
                torq_avg2= torq_avg2 + ((torq_sum2/time)*1e-5)
		avg_angle= avg_angle + ((angle_sum*180)/pi)
                avg_angle1= avg_angle1 + ((angle_sum1*180)/pi)
                avg_angle2= avg_angle2 + ((angle_sum2*180)/pi)
c               write(10,*) ((torq_sum/time)*1e-5),
c     1         ((torq_sum1/time)*1e-5),((torq_sum2/time)*1e-5)
c               write(11,*)(angle_sum*180)/pi,(angle_sum1*180)/pi,
c     1         (angle_sum2*180)/pi
                act_mot= ((act_mot/time)*1e-5)
                tk_time= time
                detach= detach/act_mot
                attach= attach/(act_mot*tk_time)
                step= step/(act_mot*tk_time)
                Fsum= (Fsum/(time*act_mot))*1e-5
                Fxsum= (Fxsum/(time*act_mot))*1e-5
                Fysum= (Fysum/(time*act_mot))*1e-5
                Fzsum= (Fzsum/(time*act_mot))*1e-5
c                do i=1,N
c                  t_array(i)= t_array(i)*del_t
c                end do
                print*, act_mot, detach, attach, step
                velocity= Rs/tk_time

                xsum= xsum + Rs
                sum_act= sum_act + act_mot
                sum_t= sum_t + tk_time
                sum_det= sum_det + detach
                sum_att= sum_att + attach
                sum_step= sum_step + step
                Favg= Favg + Fsum
                Fxavg= Fxavg + Fxsum
                Fyavg= Fyavg + Fysum
                Fzavg= Fzavg + Fzsum
                vel_sum= vel_sum + velocity
                do i=1,N
                 t_mot_avg= t_mot_avg + t_array(i)
                end do

                xsqsum= xsqsum + (Rs**2)
                sumsq_mot= sumsq_mot + (act_mot**2)
                tsq_sum= tsq_sum + (tk_time**2)
                detsq_sum= detsq_sum + ((detach/tk_time)**2)
                attsq_sum= attsq_sum + (attach**2)
                stepsq_sum= stepsq_sum + (step**2)
                vel_sqsum= vel_sqsum + (velocity**2)
                print*, "step=",step, "RL=",Rs,"time=",
     1          tk_time,"V=",velocity,"T=",time
c		print*, bead_final,(time*del_t),velocity
                if(sigconf== maxconf)then
                 do i=1,1000
                   avg_motin(i)= sum_motin(i)/iterm(i)
                   avg_area(i)= sm_area(i)/iterm(i)
                   write(21,*)-F_ext,vmot,N,i*0.2,iterm(i),
     1                  avg_motin(i), avg_area(i)
                 end do
                end if
                go to 300
c	  else
c		time= time + 1
c		act_mot = act_mot + m
          end if

            do i=1,m
                F_mot(i)=0; Fonmot(i)=0
                F_motx(i)=0;F_moty(i)=0
                F_motz(i)=0
            end do

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
c             Fsum= Fsum + Fonmot(i)
            end do

            do i=1,m
                del(i,1)=0; del(i,2)=0; del(i,3)=0
                tau(i,1)=0; tau(i,2)=0; tau(i,3)=0
                F_motx(i)=0; F_moty(i)=0; F_motz(i)=0
                u_motx(i)=0; u_moty(i)=0; u_motz(i)=0
            end do

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
            torq_sum= torq_sum + sum_torq(1)
            torq_sum1= torq_sum1 + sum_torq(2)
            torq_sum2= torq_sum2 + sum_torq(3)

3           r1= 2*rand()-1
            r2= 2*rand()-1
            gauss= r1**2+ r2**2
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

c          x_random=0; y_random=0; z_random=0
           dx= ((sumF_x + Fx_ext)*vari2) + x_random
           dy= (sumF_y*vari2) + y_random
           dz= (sumF_z*vari2) + z_random

c          tau_ran(1)=0; tau_ran(2)=0; tau_ran(3)=0
c          sum_torq(1)=0; sum_torq(2)=0; sum_torq(3)=0
           do i=1,3
                angle(i)= (sum_torq(i)*variate)+ tau_ran(i)
c               write(10,*) "angle=",angle(i)
           end do

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
c              print*, "arrested case"
c               if(sigconf==125)write(10,*) "arrested case"
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

          angle_sum= angle_sum+ angle(1)
          angle_sum1= angle_sum1+ angle(2)
          angle_sum2= angle_sum2+ angle(3)

          for_var= sqrt(sumF_x**2 +sumF_y**2 +sumF_z**2)
          Fsum= Fsum + for_var
          Fxsum= Fxsum + sumF_x
          Fysum= Fysum + sumF_y
          Fzsum= Fzsum + sumF_z

          do i=1,N
            ang_mot(i,1)=0; ang_mot(i,2)=0; ang_mot(i,3)=0
          end do

	   tau_mot1=0; tau_mot2=0; tau_mot3=0
           do i=1,N
            if (stalk_motl(i) .le. Lo)then                                            ! case of detached/unattached motor
              Fall_mot(i)=0
            else
              Fall_mot(i)= (stalk_motl(i)-Lo)*springk  			     ! magnitude of force exerted by attached and extended motors
            end if
            if (Fall_mot(i).eq. 0)then                                            ! case of detached/unattached motor
              taual1(i)=0; taual2(i)=0; taual3(i)=0
            else if(Fall_mot(i).ne.0)then
              dl_all1= copy_beadmx(i)- copy_cargox
              dl_all2= copy_beadmy(i)- copy_cargoy
              dl_all3= copy_beadmz(i)- copy_cargoz
              uall1= (copy_headmx(i)-copy_beadmx(i))/stalk_motl(i)                 ! unit vector of x,y,z line connecting head and tail of motors
              uall2= (copy_headmy(i)-copy_beadmy(i))/stalk_motl(i)
              uall3= (copy_headmz(i)-copy_beadmz(i))/stalk_motl(i)
              Fall_motx(i)= Fall_mot(i)*uall1
              Fall_moty(i)= Fall_mot(i)*uall2
              Fall_motz(i)= Fall_mot(i)*uall3
              taual1(i)= dl_all2*Fall_motz(i)-dl_all3*Fall_moty(i)
              taual2(i)= -dl_all1*Fall_motz(i)+dl_all3*Fall_motx(i)
              taual3(i)= dl_all1*Fall_moty(i)-dl_all2*Fall_motx(i)
            end if

c100         r_mo1= 2*rand()-1
c            r_mo2= 2*rand()-1
c            grad= r_mo1**2+r_mo2**2
c            if(grad.ge.1 .or. grad.eq.0)go to 100
c            def= r_mo1*sqrt(-2*log(grad)/grad)
c            tau_mot1= def*sqrt(2*vmot*del_t)
c
c101         s_mot1= 2*rand()-1
c            s_mot2= 2*rand()-1
c            hrad= s_mot1**2+ s_mot2**2
c            if(hrad.ge.1 .or. hrad.eq.0)go to 101
c            dfmot= s_mot1*sqrt(-2*log(hrad)/hrad)
c            tau_mot2= dfmot*sqrt(2*vmot*del_t)
c
c102         t_mot1= 2*rand()-1
c            t_mot2= 2*rand()-1
c            frad= t_mot1**2+t_mot2**2
c            if(frad.ge.1 .or. frad.eq.0)go to 102
c            demot= t_mot1*sqrt(-2*log(frad)/frad)
c            tau_mot3= demot*sqrt(2*vmot*del_t)

            ang_mot(i,1)= (taual1(i)*cn)+tau_mot1
            ang_mot(i,2)= (taual2(i)*cn)+tau_mot2
            ang_mot(i,3)= (taual3(i)*cn)+tau_mot3         
          end do

          do i=1,N
             anchorx(i)=0; anchory(i)=0; anchorz(i)=0
          end do

          temp1=0; temp2=0; temp3=0; delta=0
	  do i=1,N
           eul11= cos(ang_mot(i,2))*cos(ang_mot(i,3))                                            !anti-clockwise
           eul12= -cos(ang_mot(i,2))*sin(ang_mot(i,3))
	   eul13= sin(ang_mot(i,2))
          eul21= sin(ang_mot(i,1))*sin(ang_mot(i,2))*cos(ang_mot(i,3))+
     1          sin(ang_mot(i,3))*cos(ang_mot(i,1))
         eul22= -sin(ang_mot(i,1))*sin(ang_mot(i,2))*sin(ang_mot(i,3))+
     1          cos(ang_mot(i,3))*cos(ang_mot(i,1))
          eul23= -sin(ang_mot(i,1))*cos(ang_mot(i,2))
         eul31= -cos(ang_mot(i,1))*sin(ang_mot(i,2))*cos(ang_mot(i,3))+
     1          sin(ang_mot(i,3))*sin(ang_mot(i,1))
          eul32= cos(ang_mot(i,1))*sin(ang_mot(i,2))*sin(ang_mot(i,3))+
     1          cos(ang_mot(i,3))*sin(ang_mot(i,1))
           eul33= cos(ang_mot(i,1))*cos(ang_mot(i,2))

           temp1= eul11*(copy_beadmx(i)-copy_cargox)+
     1       eul12*(copy_beadmy(i)-copy_cargoy)+
     1          eul13*(copy_beadmz(i)-copy_cargoz)
           anchorx(i)= temp1+ copy_cargox

           temp2= eul21*(copy_beadmx(i)-copy_cargox)+
     1       eul22*(copy_beadmy(i)-copy_cargoy)+
     1          eul23*(copy_beadmz(i)-copy_cargoz)
           anchory(i)= temp2+ copy_cargoy

           temp3= eul31*(copy_beadmx(i)-copy_cargox)+
     1       eul32*(copy_beadmy(i)-copy_cargoy)+
     1          eul33*(copy_beadmz(i)-copy_cargoz)
           anchorz(i)= temp3+ copy_cargoz

	   copy_beadmx(i)= anchorx(i)
	   copy_beadmy(i)= anchory(i)
	   copy_beadmz(i)= anchorz(i)
           pserad=sqrt((copy_beadmx(i)-copy_cargox)**2+
     1          (copy_beadmy(i)-copy_cargoy)**2+
     1          (copy_beadmz(i)-copy_cargoz)**2)
           delta= pserad-beadRadius
           if(abs(delta) .gt. 1e-6)then
                ux= (copy_beadmx(i)-copy_cargox)/beadRadius
                uy= (copy_beadmy(i)-copy_cargoy)/beadRadius
                uz= (copy_beadmz(i)-copy_cargoz)/beadRadius
                copy_beadmx(i)= copy_beadmx(i)- (ux*delta)
                copy_beadmy(i)= copy_beadmy(i)- (uy*delta)
                copy_beadmz(i)= copy_beadmz(i)- (uz*delta)
             end if
             if(stalk_motl(i)==0)then
                copy_headmx(i)= copy_beadmx(i)
                copy_headmy(i)= copy_beadmy(i)
                copy_headmz(i)= copy_beadmz(i)
             end if
          end do

         do i=1,N
            s1=((copy_beadmx(i)-copy_headmx(i))**2+
     1         (copy_beadmy(i)-copy_headmy(i))**2+
     1         (copy_beadmz(i)-copy_headmz(i))**2)
             stalk_motl(i)= sqrt(s1)
	 end do

500         t11= cos(angle(2))*cos(angle(3))						!anti-clockwise
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

c500	   t11= cos(angle(3))*cos(angle(2))
c	   t12= sin(angle(3))*cos(angle(1))+
c     1		cos(angle(3))*sin(angle(2))*sin(angle(1))
c	   t13= sin(angle(3))*sin(angle(1))-
c     1		cos(angle(3))*sin(angle(2))*cos(angle(1))
c	   t21= sin(angle(3))*cos(angle(2))
c	   t22= cos(angle(3))*cos(angle(1))-
c     1		sin(angle(3))*sin(angle(2))*sin(angle(1))
c	   t23= cos(angle(3))*sin(angle(1))+
c     1		sin(angle(3))*sin(angle(2))*cos(angle(1))
c	   t31= sin(angle(2))
c	   t32= -cos(angle(2))*sin(angle(1))
c	   t33= cos(angle(2))*cos(angle(1))

	   do i=1,N
		anchorx(i)=0; anchory(i)=0; anchorz(i)=0
	   end do

	   temp1=0; temp2=0; temp3=0; delta=0
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

	    final_cargox= copy_cargox + dx
	    final_cargoy= copy_cargoy + dy
	    final_cargoz= copy_cargoz + dz
c           C_MT_dist=sqrt(final_cargoy**2+final_cargoz**2)

	    do i=1,N
		copy_beadmx(i)= dx + anchorx(i)							!+ copy_beadmx(i)					
		copy_beadmy(i)= dy + anchory(i)							!+ copy_beadmy(i)					
		copy_beadmz(i)= dz + anchorz(i)							!+ copy_beadmz(i)					
             pseudo_rad=sqrt((copy_beadmx(i)-final_cargox)**2+
     1          (copy_beadmy(i)-final_cargoy)**2+
     1          (copy_beadmz(i)-final_cargoz)**2)
	     delta= pseudo_rad-beadRadius
	     if(abs(delta) .gt. 1e-6)then
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
	    end do

           do i=1,N
           bead_MTdist(i)= sqrt((copy_beadmz(i)**2))
c           print*, "dist=", bead_MTdist(i)
           if(bead_MTdist(i) .le. Lo)then
            if(abs(copy_beadmy(i)) .gt. 0.012)then
                i_mot= int((copy_beadmx(i)*1000)/8)
                icheck= i_mot*8
                icheck1= (i_mot+1)*8
                xcheck= sqrt((abs(copy_beadmy(i))-0.012)**2+
     1          (copy_beadmx(i)-(icheck*0.001))**2+copy_beadmz(i)**2)
                xcheck1= sqrt((abs(copy_beadmy(i))-0.012)**2+
     1          (copy_beadmx(i)-(icheck1*0.001))**2+copy_beadmz(i)**2)
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
                sort(1)= sqrt((copy_beadmy(i)-(ichecky*0.001))**2+
     1          (copy_beadmx(i)-(icheck*0.001))**2+copy_beadmz(i)**2)
                x(1)=icheck; y(1)=ichecky
                sort(2)= sqrt((copy_beadmy(i)-(ichecky1*0.001))**2+
     1          (copy_beadmx(i)-(icheck*0.001))**2+copy_beadmz(i)**2)
                x(2)= icheck; y(2)= ichecky1
                sort(3)= sqrt((copy_beadmy(i)-(ichecky*0.001))**2+
     1          (copy_beadmx(i)-(icheck1*0.001))**2+copy_beadmz(i)**2)
                x(3)= icheck1; y(3)= ichecky
                sort(4)= sqrt((copy_beadmy(i)-(ichecky1*0.001))**2+
     1          (copy_beadmx(i)-(icheck1*0.001))**2+copy_beadmz(i)**2)
                x(4)= icheck1; y(4)= ichecky1
                if(sort(1) .le. Lo .or. sort(2) .le. Lo .or. 
     1		sort(3) .le. Lo .or. sort(4) .le. Lo) then
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
		dup_stalk(k)= stalk_motl(i)
		end if
	    end do
	    m=k	
		
	    do i=1,m
	      rel_motbx(i)= dup_motbx(i)
	      rel_motby(i)= dup_motby(i)
	      rel_motbz(i)= dup_motbz(i)
	      rel_mothx(i)= dup_mothx(i)
	      rel_mothy(i)= dup_mothy(i)
	      rel_mothz(i)= dup_mothz(i)
	    end do

            do i=1,N                                                            ! updating original array of motor head and tail (bead) positions
	     pseudo_rad=0
             beadmotx(i)= copy_beadmx(i)
             beadmoty(i)= copy_beadmy(i)
             beadmotz(i)= copy_beadmz(i)
             headmotx(i)= copy_headmx(i)
             headmoty(i)= copy_headmy(i)
	     headmotz(i)= copy_headmz(i)
	    end do

            itime= int(time*1e5)
            counter= mod(itime,20000)

            if(counter==0)then
            counterr1= itime
            if(counterr1 .ne. counterprev)then
            itcount= itcount +1
            iterm(itcount)= iterm(itcount)+1
            sum_x=0

            do i=1,N
             do j=1,N
               conss= (copy_beadmx(i)-copy_beadmx(j))**2+
     1          (copy_beadmy(i)-copy_beadmy(j))**2+
     1          (copy_beadmz(i)-copy_beadmz(j))**2
                sum_x= sum_x + sqrt(conss)
             end do
             if(i==N)then
               if(N .gt. 2)then
                sum_x= sum_x/(N*(N-1))
               else
                sum_x= sum_x/2
               end if
               sum_motin(itcount)= sum_motin(itcount)+ sum_x
             end if
c             write(20,*)-F_ext,sigconf,time,i,copy_beadmx(i),
c     1          copy_beadmy(i),copy_beadmz(i)
            end do

            do i=1,N
             poly_x(i)=0
             poly_y(i)=0
            end do

c	    if(m .le. 2)then
c		total_area= 0
c		go to 120
c	     if(m==1)then
c	       total_area= pi*0.008*0.008
c	       total_area=0
c	       go to 120
c	     else if(m==2)then
c               twomot= (rel_motbx(1)-rel_motbx(2))**2+
c     1          (rel_motby(1)-rel_motby(2))**2+
c     1          (rel_motbz(1)-rel_motbz(2))**2
c	       tworut= sqrt(twomot)
c	       total_area= tworut*0.016
c	       go to 120
c	     end if
c	    end if

            temp_x= copy_beadmx(1); temp_y= copy_beadmy(1)
            do i=2,N                                                ! finding the point with smallest y-coordinate
             if(copy_beadmy(i) .lt. temp_y)then
               temp_y= copy_beadmy(i)
               temp_x= copy_beadmx(i)
             end if
            end do

            chosen_x= temp_x                                        !chosen point acc to smallest y-coordinate
            chosen_y= temp_y
            poly_x(1)=chosen_x                                      !initialize the point to the array of the polygon vertices
            poly_y(1)=chosen_y
            start_x= chosen_x                                       !store point in the start point variable
            start_y= chosen_y
            prev_x= 0
            prev_y= chosen_y

            im=1; ik=0; item=0

110         ik =ik +1
            do i=1,N
            if(copy_beadmx(i)==chosen_x .and.
     1          copy_beadmy(i)==chosen_y)then
              agl(i)= 1e-5
            else
              conss1= (chosen_y-prev_y)/(chosen_x-prev_x)                                                         !slope of line
              conss6=(copy_beadmy(i)-chosen_y)/(copy_beadmx(i)-chosen_x)
              conss7= (conss6-conss1)/(1+(conss6*conss1))
             if(conss7 .ge. 0)then
              agl(i)= atan(conss7)*v
             else
              agl(i)= 180+(atan(conss7)*v)
             end if
            end if

            do j=2,im
             if(copy_beadmx(i)==poly_x(j) .and.
     1          copy_beadmy(i)==poly_y(j))then
              agl(i)= 1e-5
             end if
            end do
            end do

            it= 0
            do i=1,N
            if(agl(i) == 1e-5)then
              it= it+1
            end if
            end do

            if(it == (N-1))then
             do i=1,N
              if(agl(i) .gt. 1e-5)then
                im= im+1
                poly_x(im)= copy_beadmx(i)
                poly_y(im)= copy_beadmy(i)
                prev_x= chosen_x
                prev_y= chosen_y
                chosen_x= copy_beadmx(i)
                chosen_y= copy_beadmy(i)
              end if
             end do
            go to 600
            end if

            do i=1,N
              if(agl(i) .gt. 1e-5)then
                temp_angle= agl(i)
                item= i
                exit
              end if
           end do

           do j=1,N
           if (agl(j) .gt. 1e-5 .and. agl(j) .lt. temp_angle)then
             temp_angle= agl(j)
             item = j
           end if
           end do
           im= im+1

           poly_x(im)= copy_beadmx(item)
           poly_y(im)= copy_beadmy(item)
           prev_x= chosen_x
           prev_y= chosen_y
           chosen_x= copy_beadmx(item)
           chosen_y= copy_beadmy(item)

           if(poly_x(im).eq. start_x .and. poly_y(im).eq. start_y)then
            go to 600
           end if

           it= 0
           do i=1,N
            if(agl(i) == 1e-5)then
             it= it+1
            end if
           end do
           if(it .lt. (N-1))then
                 go to 110
           end if

600	   continue

           sum_area1= poly_x(im-1)*poly_y(1)
           sum_area2= poly_x(1)*poly_y(im-1)
           do i=1,im-2
             sum_area1= sum_area1+ poly_x(i)*poly_y(i+1)
             sum_area2= sum_area2+ poly_x(i+1)*poly_y(i)
           end do
           total_area= (abs(sum_area1 - sum_area2))/2
           sm_area(itcount)= sm_area(itcount)+ total_area
           write(20,*)-F_ext,N,vmot,sigconf,time,total_area,sum_x
           end if

           counterprev= itime
           end if

           copy_cargox= final_cargox
           copy_cargoy= final_cargoy
           copy_cargoz= final_cargoz

	  end do
300       continue
	  print*, "significant conf=", sigconf
	  if (sigconf==maxconf)exit  
	 end do

	 avg_angle = avg_angle/sigconf
         avg_angle= mod(avg_angle,360.)	 
         avg_angle1 = avg_angle1/sigconf
	 avg_angle1= mod(avg_angle1,360.)
         avg_angle2 = avg_angle2/sigconf
         avg_angle2= mod(avg_angle2,360.)
      	 torq_avg= torq_avg/sigconf
         torq_avg1= torq_avg1/sigconf
         torq_avg2= torq_avg2/sigconf
	 print*, "avg torque=", torq_avg,"avg angle=",avg_angle
         print*, "avg torque=", torq_avg1,"avg angle=",avg_angle1
         print*, "avg torque=", torq_avg2,"avg angle=",avg_angle2
c        write(10,* ) "Clockwise rotn matrix, thermal F=0"
c        write(10,*) "Code-"
c      write(10,*)"cluster_rotate_lin_motion_classical_MC_motor_dynamics"
c       write(10,*) "Average Torque per unit time for F=0 for 1000 conf="
         write(10,*) iF,N,imtd,beadRadius,torq_avg,torq_avg1,torq_avg2
c        write(11,*) "Clockwise rotn, thermal F=0"
c        write(11,*) "Code-"
c      write(11,*)"cluster_rotate_lin_motion_classical_MC_motor_dynamics"
c        write(11,*) "avg angular displacement at F=0 for 1000 conf="
        write(11,*) iF,N,imtd,beadRadius,avg_angle,avg_angle1,avg_angle2

         RLavg= xsum/sigconf
         avgact= sum_act/sigconf
         avgt= sum_t/sigconf
         detavg= sum_det/(avgt*sigconf)
         attavg= sum_att/sigconf
         stpavg= sum_step/sigconf
         Favg= Favg/sigconf
         Fxavg= Fxavg/sigconf
         Fyavg= Fyavg/sigconf
         Fzavg= Fzavg/sigconf
         t_mot_avg= t_mot_avg/sigconf
         
         Xavgsq= xsqsum/sigconf
         SD= sqrt(Xavgsq-(RLavg**2))
         SEM= SD/sqrt(sigconf)

         Actsq_avg= sumsq_mot/sigconf
         active_SD= sqrt(Actsq_avg-(avgact**2))
         active_SEM= active_SD/sqrt(sigconf)

         tavg_sq= tsq_sum/sigconf
         tavg_SD= sqrt(tavg_sq-(avgt**2))
         tavg_SEM= tavg_SD/sqrt(sigconf)

         detsq_avg= detsq_sum/sigconf
         det_SD= sqrt(detsq_avg-(detavg**2))
         det_SEM= det_SD/sqrt(sigconf)

c         attsq_avg= attsq_sum/sigconf
c         att_SD= sqrt(attsq_avg-(attavg**2))
c         att_SEM= att_SD/sqrt(sigconf)

         stpsq_avg= stepsq_sum/sigconf
         step_SD= sqrt(stpsq_avg-(stpavg**2))
         step_SEM= step_SD/sqrt(sigconf)

         avgvel= vel_sum/sigconf
         vel_sqavg= vel_sqsum/sigconf
         vel_SD= sqrt(vel_sqavg- (avgvel**2))
         vel_SEM= vel_SD/sqrt(sigconf)

         print*, -F_ext, "detach= ",detavg," attach= ",attavg," step= ",
     1          stpavg, sigconf, " RL=",RLavg," active motors=", avgact
         print*, "Motor t=", t_mot_avg,"Cargo t=", avgt,"V=",avgvel

         write(22,*) -F_ext,Fxavg,N,imtd,beadRadius,detavg,det_SD,
     1          det_SEM,stpavg,step_SD,step_SEM,avgvel,vel_SD,vel_SEM

         write(23,*)-F_ext,Fxavg,N,imtd,beadRadius,RLavg,SD,SEM,avgt,
     1          tavg_SD,tavg_SEM,avgact,active_SD,active_SEM,
     1          t_mot_avg

	end do
	end do
	close(22)
	close(23)
	close(11)
	close(10)
	stop
	end
