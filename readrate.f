	Program ratiorates
	open(unit=21, file= "Rate_V_clusmot_D001_N10.dat")
        open(unit=22, file= "Rate_V_clusmot_D01_N10.dat")
        open(unit=24, file= "Rate_V_clusmot_D1_N10.dat")
	open(unit=23, file= "ratios_rates_clusmot_V_N10_F0.dat")

	do i=1,33
	 if(i .le. 11)read(21,*)iF,N,imtd,beadRadius,detavg,
     1		det_SD,det_SEM,
     1          attavg,att_SD,att_SEM,stpavg,step_SD,step_SEM
         if(i .gt. 11 .and. i .le. 22)read(22,*)iF,N,imtd,
     1		beadRadius,detavg,det_SD,det_SEM,
     1          attavg,att_SD,att_SEM,stpavg,step_SD,step_SEM
         if(i .gt. 22 .and. i .le. 33)read(24,*)iF,N,imtd,
     1		beadRadius,detavg,det_SD,det_SEM,
     1          attavg,att_SD,att_SEM,stpavg,step_SD,step_SEM
	if(i==1 .or. i==12 .or. i==23)write(23,*)iF,N,imtd,
     1		beadRadius,detavg/attavg,detavg,det_SD,det_SEM,
     1          attavg,att_SD,att_SEM,stpavg,step_SD,step_SEM
	end do
	close(21)
	close(22)
	close(23)
	close(24)
	End Program
