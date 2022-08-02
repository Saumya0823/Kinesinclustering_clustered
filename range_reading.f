	program readtimetails
	parameter (N= 10000)
	dimension arr(N),templ_IMD(N),templ_area(N),tempup_IMD(N),
     1		tempup_area(N),arr_IMD(N),arr_ar(N),ran_IMD(N),
     1		ran_ar(N)
        open(unit=20, file="motor_tails_clusmotran_N5_D001_F0.dat")
	open(unit=22, file="time_points_clusmotran_N5_D001_F0.dat")
        open(unit=21, file="range_clusmotran_D001_N5_F0.dat")
	
	do i=1,188,1
	  arr(i)= 0
          templ_IMD(i)=0.5
          templ_area(i)=0.2
          tempup_IMD(i)=0
          tempup_area(i)=0	
	end do

	j=0
	do i=1,48265
         read(20,*)F_ext,M,vmot,sigconf,time,total_area,sum_x
	 if(time==0.200008288 )j=1
	 arr(j)= time
	 arr_IMD(j)= sum_x
	 arr_ar(j)= total_area
c	 print*, j, time, arr(j), arr_IMD(j), arr_ar(j)
	 if(tempup_IMD(j) .lt. arr_IMD(j))tempup_IMD(j)=arr_IMD(j)
         if(templ_IMD(j) .gt. arr_IMD(j))templ_IMD(j)=arr_IMD(j)
         if(tempup_area(j) .lt. arr_ar(j))tempup_area(j)=arr_ar(j)
         if(templ_area(j) .gt. arr_ar(j))templ_area(j)=arr_ar(j)
c         print*, j,time,arr(j),arr_IMD(j),templ_IMD(j),tempup_IMD(j)
c         print*, j,time,arr(j),arr_ar(j),templ_area(j),tempup_area(j)
	 j= j+1
	end do
	print*, j

	do i=1,188,1
	 read(22,*)Fn,vmot,M,time,ivar1,distIMD,Areaar
	 ran_IMD(i)= tempup_IMD(i)- templ_IMD(i)
	 ran_ar(i)= tempup_area(i)-templ_area(i)
	 print*, arr(i), ran_IMD(i), ran_ar(i)
         write(21,*) arr(i), distIMD,ran_IMD(i),Areaar,ran_ar(i)
c	 write(21,*) arr(i), ran_IMD(i), ran_ar(i),tempup_IMD(i),
c     1		templ_IMD(i),tempup_area(i),templ_area(i)
	end do

	close(21)
	close(20)
	close(22)
	stop
	end
