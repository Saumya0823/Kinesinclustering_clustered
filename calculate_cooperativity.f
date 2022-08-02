	Program cooperativity index
	parameter (Vo=1000, Fs=5.,w=2)
	open(unit=23,file="RLV_clusmotran_N10_D001.dat")
	open(unit=24,file="V_clus_motran_D001_N10.dat")
	open(unit=25, file="coop_clus_motran_N10_D001.dat")

	do iF=0,10,1
	Fext=-float(iF)
	read(23,*)s,Fx,j,k,Rad,RL,RLSD,RLSEM,t,tSD,tSEM,
     1		atmot,atSD,atSEM,tmot
	read(24,*)l,m,n,Radii,V,VSD,VSEM,amot,aSD,aSEM

c	if(abs(Fext) .lt. Fs)then
c          th_V= (1-(abs(Fext/Fs)**w))
c	else
c	  var1= 1.07+ abs(0.186*Fext)
c	  th_V= 1/var1
c	end if

	tmot= tmot*1e5
	V_frac= V/1.00037169						!0.994736433		!/th_V
	variable= (1-(t/tmot))
	print*, V_frac, variable
	if(variable .lt. 1e-5)then
	  coop=0
	else
	  coop=variable*V_frac
	end if
	print*, -Fext,Rad,RL,V,th_V,V_frac,t,tmot,coop
	write(25,*)-Fext,j,k,Rad,RL,V,th_V,V_frac,t,tmot,coop
	end do
	close(23)
	close(24)
	close(25)
	stop
	end
