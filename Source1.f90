module CaviData
implicit none
	integer nx,ny
	integer case_study
	real (8) v_x(0:8),v_y(0:8)    
	real (8) wi(0:8)        
	logical, ALLOCATABLE  :: Solid(:,:)
	real (8), allocatable :: f(:,:,:),f_help(:,:,:), g(:,:,:),g_help(:,:,:)
	real (8), allocatable :: mixctr(:,:),ux_mid(:),ux_mid_his(:),T_mid(:),T_mid_his(:)
    integer BF_ctr,time
    real (8) pr,ra,da,a1,a2,a3,a4, density, teta, T0
    real (8) delta_x,delta_t,speed,Cs2,tau,tau_g,speed2,speed4
    real (8) nu_his,nu_change,ave_nuselt,stime,rtime
	real (8) Degree
    real (8) T_change,Ux_change
    real (8) Fraction
	real (8) e,visco,k1,c_0,c_1,u_x_0,u_y_0,u_x_1,u_y_1,cf,cf_1,L0
	real (8) rho_s,rho_f,Cp_s,Cp_f,sigmaa,m1,m2
    real (8) N_mac,e_mac,D_mac,delta,free_space,porous_space,block_size
    real (8) N_box,deltax,deltay,Lm,LOD,Lx_box,Ly_box,Lm_y
	real (8) tau_1
	real (8) pr_e,visco_f,K_f,k_s,tau_0,T_c,T_h,tau1,tau_loc
    character (len=300) :: FileNameSTR
contains
	function case_tag() result(tag)
	    implicit none
	    character (len=32) :: tag
	    character (len=16) :: case_number

	    write(case_number, '(I0.2)') case_study
	    tag='Case'//trim(case_number)//'_'
	end function

	subroutine ini_var(n_x,n_y)  ! allocation general variables
		implicit none
		integer , intent(in) :: n_x,n_y

		v_x(0)=0.d0; v_x(1)=1.d0; v_x(2)=0.d0; v_x(3)=-1.d0; v_x(4)=0.d0
		v_y(0)=0.d0; v_y(1)=0.d0; v_y(2)=1.d0; v_y(3)=0.d0 ; v_y(4)=-1.d0

		v_x(5)=1.d0; v_x(6)=-1.d0; v_x(7)=-1.d0; v_x(8)=1.d0
		v_y(5)=1.d0; v_y(6)=1.d0 ; v_y(7)=-1.d0; v_y(8)=-1.d0

		wi(0)=4.d0/9.d0
		wi(1)=1.d0/9.d0;  wi(2)=1.d0/9.d0;  wi(3)=1.d0/9.d0;  wi(4)=1.d0/9.d0
		wi(5)=1.d0/36.d0; wi(6)=1.d0/36.d0; wi(7)=1.d0/36.d0; wi(8)=1.d0/36.d0

		allocate (Solid(0:n_x-1,0:n_y-1),f(0:nx-1,0:ny-1,0:8),f_help(0:nx-1,0:ny-1,0:8),g(0:nx-1,0:ny-1,0:8),g_help(0:nx-1,0:ny-1,0:8),mixctr(0:nx-1,0:ny-1))
        allocate (ux_mid(0:n_y-1),ux_mid_his(0:n_y-1),T_mid(0:n_y-1),T_mid_his(0:n_y-1))
		solid(:,:)=.false.
        ux_mid_his(:)=1000.0;T_mid_his(:)=1000.0
		f(:,:,:)=0.0;f_help(:,:,:)=0.0;g(:,:,:)=0.0;g_help(:,:,:)=0.0;mixctr(:,:)=0.d0
		print*,'--> general variables are allocated'
	end subroutine
end module CaviData
!________________________________________________________________________________________________________________________
      program prog
      use CaviData
	  integer run_ctr, cmd_status, cmd_exit
      print*,' 1- natural convection in cavity'
      print*,' 2- natural convectin with temperature-dependent viscosity'
      print*;print*;
      read*,run_ctr
      call execute_command_line('mkdir -p result', cmdstat=cmd_status, exitstat=cmd_exit)
      if (cmd_status /= 0) error stop 'Unable to create result directory'
      if (cmd_exit /= 0) error stop 'Unable to create result directory'
      select case(run_ctr)
      case(1)
          call singlerun
     
      case(2)
          call singlerun2
      end select
     end program prog
!________________________________________________________________________________________________________________________      
!________________________________________________________________________________________________________________________
!________________________________________________________________________________________________________________________
subroutine singlerun
    use CaviData
    implicit none
    integer time_ctr
    Print*,'------------------------------------'
    print*,'          Single RUN                '
    Print*,'------------------------------------'
        Print*,' Input required information ...';        print*;
        print*,' nx, ny (for example 81,81):'
        read*,nx,ny;        print*;
        print*,' BF:    1) Natural    2) Natural+Porous   3)bidisperse porous mediume ' ;    print*;
        read*, Bf_ctr;    print*;
        call Ini_Calc
        call outputfiles
        call geo_maker
        call cpu_time(stime)
        time=0   
        open(77,file='result/'//trim(case_tag())//'Nu_Time.txt')
        do while ((T_change >= 1e-3) .OR. (Ux_change >= 1e-5) )
            time=time+1
            time_ctr=int(0.02/delta_t)
            if (time .ge. 10 .and. mod(time,time_ctr) .eq. 0) then  
                  call check_density
                  call converg_Check                  
	              call NU_calc
	              call write_results
                  write(77,*) time*delta_t,ave_nuselt
             end if 
            call stream
            call fluid_bc
            call temprature_bc
            call collision
        end do
        close(77)
        print*;print*;
        print*,'----- FINISHED !!!! -----'    
        write(20,*) '----------------------------'
        write(20,*) ' Final Nu(Ave) :', ave_nuselt
        close(20)
        pause
    end subroutine
	!_____________________________________________________________________________
	subroutine singlerun2
    use CaviData
    implicit none
    integer time_ctr
    Print*,'-----------------------------------------------------'
    print*,'     Single RUN,temperature-dependent viscosity                '
    Print*,'-----------------------------------------------------'
        Print*,' Input required information ...';        print*;
        print*,' nx, ny (for example 81,81):'
        read*,nx,ny;        print*;
        print*,' BF:    1) Natural    2) Natural+Porous  3)bidisperse porous mediume ' ;    print*;
        read*, Bf_ctr;    print*;
        call Ini_Calcvisco
        call outputfiles
        call geo_maker
        call cpu_time(stime)
        time=0   
        open(77,file='result/'//trim(case_tag())//'Nu_Time.txt')
        do while ((T_change >= 1e-5) .OR. (Ux_change >= 1e-5) )
            time=time+1
            time_ctr=int(0.02/delta_t)
            if (time .ge. 10 .and. mod(time,time_ctr) .eq. 0) then  
                  call check_density
                  call converg_Check                  
	              call NU_calc
	              call write_results
                  write(77,*) time*delta_t,ave_nuselt
             end if 
            call stream
            call fluid_bc
            call temprature_bc
            call collisionvisco
        end do
        close(77)
        print*;print*;
        print*,'----- FINISHED !!!! -----'    
        write(20,*) '----------------------------'
        write(20,*) ' Final Nu(Ave) :', ave_nuselt
        close(20)
    end subroutine
	!_____________________________________________________________________________

	subroutine Ini_calc
    use CaviData
    implicit none
    integer x,y,k
    call ini_var(nx,ny)
    Print*,'--> Reading Data from Input.txt'
    open(12,file='Input.txt')
	read (12,*)  !case_study
    read (12,*)  case_study
    read (12,*) !a1,a2
    read (12,*) a1,a2
    read (12,*) !a3,a4
    read (12,*) a3,a4
    read (12,*) ! pr
    read (12,*) pr
    read (12,*) ! Ra
    read (12,*) Ra
    read (12,*) ! Da
    read (12,*) Da
    read (12,*) ! micro porosity
	read (12,*) e
	read (12,*) ! degree of gravity
	read (12,*) degree
    close(12)
    print*,'--> Calculating other parameters'
    density=1.d0
    teta=(degree)*3.1415926535897932384626433832795/180.0   ! (Degree)
    T0=1.0    
    delta_x=1.0/real(ny-1)
    tau=0.833
    delta_t=(tau-0.5)*delta_x**2/3.0/Pr
    speed=delta_x/delta_t
    speed2=speed*speed
    speed4=speed2*speed2
    Cs2=speed*speed/3.0
    tau_g=3.0*delta_t/delta_x/delta_x+0.5
	!------------------------------------------------------------
	!e=0.5                                             !porosity
	rho_s=1
	rho_f=1
	Cp_s=1
	Cp_f=1
	sigmaa=e+(1-e)*rho_s*Cp_s/rho_f/Cp_f
	L0=1
	visco=Cs2*(tau-0.5)*delta_t
	Cf_1=150*(e**3)
    Cf=1.75/sqrt(Cf_1)
    K1=(L0**2)*Da
    C_0=0.5*(1.0+e*Delta_t*visco/2/K1)
    C_1=e*delta_t*Cf/2/sqrt(K1)
    !------------------------------------------------------------
   
	!------------------------------------------------------------
    nu_his=100.0
    nu_change=100.0
    T_change=100.0
    Ux_change=100.0
    Print*;print*;
    Print*,'               ----- '
    print*,' Nx, Ny        :',Nx,Ny
    Print*,' speed, teta   :',speed,teta
    print*,' tau, tau_g    :', tau,tau_g
    print*,' Pr            :',pr
    print*,' Da            :',Da
    print*,' Ra            :',Ra
    Print*,'               ----- ';  Print*;print*;
    print*,'--> Initialaize f,f_help,g,g_help'
    do x=0,nx-1
        do y=0,ny-1
            do k=0,8
                f(x,y,k)=wi(k)*density
                f_help(x,y,k)=wi(k)*density
                g(x,y,k)=wi(k)*0.0
                g_help(x,y,k)=wi(k)*0.0
            end do
        end do
    end do
    end subroutine
		
    !---------------------------------------------------------------------------
    
	subroutine outputfiles
    use CaviData
    implicit none
          character (len=300) :: filename1,pr_str,ra_str,Da_str,nx_str
          write (pr_str, "(A3,F4.2)") "_Pr", pr
          write (Ra_str, "(A3,F8.1)") "_Ra", Ra
          write (Da_str, "(A3,F8.6)") "_Da", Da
          write (Nx_str, "(A3,I3)") "_Nx", nx
          filename1=trim(case_tag())//'Results'// trim(adjustl(Nx_str))// trim(adjustl(Pr_str)) //  trim(adjustl(Da_str)) // trim(adjustl(Ra_str))//'.txt'
          open(20,file='result/'//trim(filename1))
          write(20,*) ' Run Information:'
          select case(bf_ctr)
            case(1)
                print*,'-->  BF is adjusted for pure natural convection'
                write(20,*) '-->  BF is adjusted for pure natural convection'                
            case(2)
                print*,'-->  BF is adjusted for natural convection in porous media'        
                write(20,*) '-->  BF is adjusted for natural convection in porous media'                  
            case(3)
                print*,'-->  BF is adjusted for natural convection in bidisperse porous media'  
                write(20,*) '-->  BF is adjusted for natural convection in bidisperse porous media'                 
          end select    
	      write(20,*) ' Nx,Ny   :',Nx,Ny 
	      write(20,*) ' a1 , a2 :',a1,a2
	      write(20,*) ' a3 , a4 :',a3,a4
          write(20,*) ' Tau_g   :',tau_g
          write(20,*) ' Tau     :',tau
          write(20,*) ' delta_t :',delta_t
          write(20,*) ' Pr      :',pr
          write(20,*) ' Da      :',da
          write(20,*) ' Ra      :',Ra
          print*,'-----------------------'
          !$OMP PARALLEL
          print*,'        Core ready ... '
          !$OMP END PARALLEL
          print*,'-----------------------'          
    end subroutine
	!______________________________________________________________________________
	
	!___________________________________________________________________________________________

	subroutine check_density
    use CaviData
    implicit none
    integer  x,y,k,hr,min,sec
    real (8) n_sum
    n_sum = 0.d0
    do x = 0, nx-1
        do y = 0, ny-1
            do k = 0, 8
                n_sum = n_sum + f_help(x,y,k)
            end do
        end do
    end do
    call cpu_time(rtime)
	hr=int((rtime-stime)/3600.d0)
	min=int(( (rtime-stime)-real(hr)*3600.d0 )/60.d0)
	sec=int( (rtime-stime)-real(hr)*3600.d0-real(min)*60.d0 ) 
    print*;print*;print*;
    Print*,'***********************', nx, '**********************'
	print*,' RunTime:'
	print*,hr,' Hours +',min,' Minutes +',sec,' Seconds'
    print*;
    print*, '-Iteration         :', time
    print*, '-Dimensionless time:', time*delta_t
    print*, '-Integral density :', n_sum
	print*, ' MLUPS            :',real(nx*ny)*real(time)/1.e6/(rtime-stime)         
    end subroutine
	!__________________________________________________________________________________
	subroutine Converg_Check
    use CaviData
    implicit none
    integer x,y,k,kk
    real (8) Density_loc, u_x,T_loc,T_err,Ux_err,u_y,bf(0:1),u_x_sum,u_y_sum
    x=int(real(nx-1)/2.d0)
    T_err=0.0
    UX_err=0.0
        do y=0,ny-1
				!-------------------------------- local calculations	
				select case(bf_ctr)
				case(1)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x_1=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y_1=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
				case(2)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
            do k=0,8
                density_loc=density_loc+f_help(x,y,k)
                T_loc=T_loc+g_help(x,y,k)
                u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
                u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
            end do
            Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
			Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))
				 
				case(3)!-------------bidisperse porous-------------------
				
				if ( (mixctr(x,y)==1)) then
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)


					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))
				else
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x_1=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y_1=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
				end if

				end select

				!--------------------------------------------------------
                ux_mid(y)=u_x_1
                T_mid(y)=T_loc
                T_err=T_err+abs(T_mid(y)-T_mid_his(y))
                Ux_err=Ux_err+abs(ux_mid(y)-ux_mid_his(y))              
                ux_mid_his(y)=ux_mid(y)
                T_mid_his(y)=T_mid(y)               
        end do
    Ux_change=(Ux_err/real(ny))       
    T_change=(T_err/real(ny))    
    print*,'-Ux_mid  change in % ',Ux_change*100.0
    print*,'-T_mid   change in % ',T_change*100.0
    end subroutine
	!________________________________________________________________________________________
	
	!____________________________________________________________________________________________________
	subroutine stream
    use CaviData
    implicit none
    integer x,y,x_next,y_next,k
    !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(x,y,x_next,y_next,k)
    do x=0,nx-1
        do y=0,ny-1
            do k=0,8        
                x_next=mod(x+int(v_x(k))+nx,nx)
                y_next=mod(y+int(v_y(k))+ny,ny)
                f_help(x_next,y_next,k)=f(x,y,k)
                g_help(x_next,y_next,k)=g(x,y,k)
			end do
        end do
    end do
    !$OMP END PARALLEL DO
    end subroutine
	!___________________________________________________________________________________________________
	
	!__________________________________________________________________________________________________-
	subroutine collision
    use CaviData
    implicit none
    integer x,y,k,kk
    real (8) Density_loc, u_x,u_y,T_loc,f_eq,g_eq,p,BF(0:1),GF,h(0:1),u_x_sum,u_y_sum!,u_x_0,u_y_0,u_x_1,u_y_1

    !============================== flow
    !$OMP PARALLEL
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(x,y,k,kk,Density_loc,u_x,u_y,T_loc,f_eq,g_eq,p,BF,GF,h,u_x_sum,u_y_sum,u_x_0,u_y_0,u_x_1,u_y_1)
    do x=0,nx-1
        do y=0,ny-1
           ! if (.not. solid(x,y)) then
				!-------------------------------- local calculations	
				select case(bf_ctr)
				case(1)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x+speed*v_y(k)*u_y
                    f_eq=wi(k)*density_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4-1.5*(u_x*u_x+u_y*u_y)/speed2)                  
					h(0)=(speed*v_x(k)-u_x)/cs2+p*speed*v_x(k)/cs2/cs2
					h(1)=(speed*v_y(k)-u_y)/cs2+p*speed*v_y(k)/cs2/cs2						
                    Gf=delta_t*(1.d0-1.d0/2.d0/tau)*wi(k)*(bf(0)*h(0)+bf(1)*h(1))
                    f(x,y,k) = f_help(x,y,k) + ( f_eq - f_help(x,y,k) ) / tau + Gf
                end do


				case(2)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)


					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))
				
		
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x_1+speed*v_y(k)*u_y_1
                    f_eq=wi(k)*density_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4/e-1.5*(u_x_1*u_x_1+u_y_1*u_y_1)/speed2/e) 
					Bf(0)=-(Pr/Da*u_x_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_x_1)*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-(Pr/Da*u_y_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_y_1)*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)               
					h(0)=(speed*v_x(k)-u_x_1)/cs2+p*speed*v_x(k)/cs2/cs2
					h(1)=(speed*v_y(k)-u_y_1)/cs2+p*speed*v_y(k)/cs2/cs2						
                    Gf=delta_t*(1.d0-1.d0/2.d0/tau)*wi(k)*(bf(0)*h(0)+bf(1)*h(1))
                    f(x,y,k) = f_help(x,y,k) + ( f_eq - f_help(x,y,k) ) / tau + Gf
                end do
				case(3)!-------------bidisperse porous-------------------
				
				if ( (mixctr(x,y)==1)) then
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)


					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))
				
		
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x_1+speed*v_y(k)*u_y_1
                    f_eq=wi(k)*density_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4/e-1.5*(u_x_1*u_x_1+u_y_1*u_y_1)/speed2/e) 
					Bf(0)=-(Pr/Da*u_x_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_x_1)*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-(Pr/Da*u_y_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_y_1)*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)               
					h(0)=(speed*v_x(k)-u_x_1)/cs2+p*speed*v_x(k)/cs2/cs2
					h(1)=(speed*v_y(k)-u_y_1)/cs2+p*speed*v_y(k)/cs2/cs2						
                    Gf=delta_t*(1.d0-1.d0/2.d0/tau)*wi(k)*(bf(0)*h(0)+bf(1)*h(1))
                    f(x,y,k) = f_help(x,y,k) + ( f_eq - f_help(x,y,k) ) / tau + Gf
                end do

				else


				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x+speed*v_y(k)*u_y
                    f_eq=wi(k)*density_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4-1.5*(u_x*u_x+u_y*u_y)/speed2)                  
					h(0)=(speed*v_x(k)-u_x)/cs2+p*speed*v_x(k)/cs2/cs2
					h(1)=(speed*v_y(k)-u_y)/cs2+p*speed*v_y(k)/cs2/cs2						
                    Gf=delta_t*(1.d0-1.d0/2.d0/tau)*wi(k)*(bf(0)*h(0)+bf(1)*h(1))
                    f(x,y,k) = f_help(x,y,k) + ( f_eq - f_help(x,y,k) ) / tau + Gf
                end do

            end if

		   end select

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    !============================== thermal
    !$OMP PARALLEL
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(x,y,k,kk,Density_loc,u_x,u_y,T_loc,f_eq,g_eq,p,BF,GF,h,u_x_sum,u_y_sum,u_x_0,u_y_0,u_x_1,u_y_1)
        do x=0,nx-1
        do y=0,ny-1
				!-------------------------------- local calculations	
				select case(bf_ctr)
				case(1)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x+speed*v_y(k)*u_y
                    g_eq=wi(k)*T_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4-1.5*(u_x*u_x+u_y*u_y)/speed2)                    
                    g(x,y,k) = g_help(x,y,k) + ( g_eq - g_help(x,y,k) ) / tau_g
                end do

				case(2)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))

				do k=0,8
                    p=speed*v_x(k)*u_x_1+speed*v_y(k)*u_y_1
                    g_eq=wi(k)*T_loc*(1.d0+3.d0*p/speed2)                    
                    g(x,y,k) = g_help(x,y,k) + ( g_eq - g_help(x,y,k) ) / tau_g
                end do
				case(3)!-------------bidisperse porous-------------------
				if ( (mixctr(x,y)==1)) then
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)


					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))
				else
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x_1=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y_1=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
				end if
				
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x_1+speed*v_y(k)*u_y_1
                    g_eq=wi(k)*T_loc*(1.d0+3.d0*p/speed2)                    
                    g(x,y,k) = g_help(x,y,k) + ( g_eq - g_help(x,y,k) ) / tau_g
                end do
				end select
        end do
        end do
    !$OMP END DO
    !$OMP END PARALLEL
    end subroutine
	!__________________________________________________________________________________________
	!------------------------------------------------------------------------------------------
	subroutine Ini_calcvisco
    use CaviData
    implicit none
    integer x,y,k
    call ini_var(nx,ny)
    Print*,'--> Reading Data from Input.txt'
    open(12,file='Input.txt')
    read (12,*) !case_study
    read (12,*) case_study
    read (12,*) !a1,a2
    read (12,*) a1,a2
    read (12,*) !a3,a4
    read (12,*) a3,a4
    read (12,*) ! pr
    read (12,*) pr
    read (12,*) ! Ra
    read (12,*) Ra
    read (12,*) ! Da
    read (12,*) Da
    read (12,*) ! micro porosity
    read (12,*) e
    read (12,*) ! degree of gravity
    read (12,*) degree
    close(12)
    print*,'--> Calculating other parameters'
    density=1.d0
    teta=(degree)*3.1415926535897932384626433832795/180.0   ! (Degree)
    Pr_e= 0.164
    T0=0.21
	!T_h=T0
	!T_c=0.0       
    delta_x=1.0/real(ny-1)
    tau=1.0
    delta_t=(tau-0.5)*delta_x**2/3.0/Pr_e
    speed=delta_x/delta_t
    speed2=speed*speed
    speed4=speed2*speed2
    Cs2=speed*speed/3.0
	!------------------------------------------------------------
	!------------------------------------------------------------
	Pr_e= 0.164
	!-------------------PAO(PolyAlphaolefin)---------------------
	!rho_f=768.5
	!visco_f=5.95*e-3
	!k_f=0.1424
	!Cp_f=1971.35
	!--------------Aluminum foam--------------------------------
	!rho_s=1
	!Cp_s=1
	!K_s=170.0
	!--------------------------tau_0=?---------------------------
	sigmaa=e+(1-e)*rho_s*Cp_s/rho_f/Cp_f
	tau_0=3*Pr_e*delta_t/delta_x/delta_x+0.5
	visco=Cs2*(tau_0-0.5)*delta_t
	!------------------------------------------------------------
	tau_g=(1/Pr_e)*(tau_0 - 0.5)+0.5
	!------------------------------------------------------------
	L0=1
	Cf_1=150*(e**3)
    Cf=1.75/sqrt(Cf_1)
    K1=(L0**2)*Da
    C_0=0.5*(1.0+e*Delta_t*visco/2/K1)
    C_1=e*delta_t*Cf/2/sqrt(K1)
	!------------------------------------------------------------
    nu_his=100.0
    nu_change=100.0
    T_change=100.0
    Ux_change=100.0
    Print*;print*;
    Print*,'               ----- '
    print*,' Nx, Ny        :',Nx,Ny
    Print*,' speed, teta   :',speed,teta
    print*,' tau, tau_g    :', tau,tau_g
    print*,' Pr            :',pr
    print*,' Da            :',Da
    print*,' Ra            :',Ra
    print*,' delta_t       :',delta_t
    print*,' tau_0         :',tau_0
    Print*,'               ----- ';  Print*;print*;
    print*,'--> Initialaize f,f_help,g,g_help'
    do x=0,nx-1
        do y=0,ny-1
            do k=0,8
                f(x,y,k)=wi(k)*density
                f_help(x,y,k)=wi(k)*density
                g(x,y,k)=wi(k)*0.5
                g_help(x,y,k)=wi(k)*0.5
            end do
        end do
    end do
    end subroutine

	!------------------------------------------------------------
	subroutine collisionvisco
    use CaviData
    implicit none
    integer x,y,k,kk
    real (8) Density_loc, u_x,u_y,T_loc,f_eq,g_eq,p,BF(0:1),GF,h(0:1),u_x_sum,u_y_sum !,u_x_0,u_y_0,u_x_1,u_y_1
    real (8) tau_T
    !============================== flow
    !$OMP PARALLEL
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(x,y,k,kk,Density_loc,u_x,u_y,T_loc,f_eq,g_eq,p,BF,GF,h,u_x_sum,u_y_sum,u_x_0,u_y_0,u_x_1,u_y_1,tau_loc)
    do x=0,nx-1
        do y=0,ny-1
           ! if (.not. solid(x,y)) then
				!-------------------------------- local calculations	
				select case(bf_ctr)
				case(1)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x+speed*v_y(k)*u_y
                    f_eq=wi(k)*density_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4-1.5*(u_x*u_x+u_y*u_y)/speed2)                  
					h(0)=(speed*v_x(k)-u_x)/cs2+p*speed*v_x(k)/cs2/cs2
					h(1)=(speed*v_y(k)-u_y)/cs2+p*speed*v_y(k)/cs2/cs2						
                    Gf=delta_t*(1.d0-1.d0/2.d0/tau)*wi(k)*(bf(0)*h(0)+bf(1)*h(1))
                    f(x,y,k) = f_help(x,y,k) + ( f_eq - f_help(x,y,k) ) / tau + Gf
                end do


				case(2)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)


					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed     
					              
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))
                !------------------------------------------------------
				!m1=(T_loc) /0.21
				!m2=m1**-1.0868
                !tau_1=3*Pr_e*delta_t/delta_x/delta_x
				!tau_loc=tau_1*m2+0.5 
				!tau_T=m2!*tau_g
                tau_loc=(m2)*(tau_0 - 0.5)+0.5   
           
				 !print*,   tau_loc 
                ! pause
                                                                      !!!!!!!!!!!!!!!!!relaxation time!!!!!!!!!!!!!!!!!!!
				                                                      !-------------------------------------------------!
				                                                      !-------------------------------------------------!
																	  !-------------------------------------------------!
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x_1+speed*v_y(k)*u_y_1
                    f_eq=wi(k)*density_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4/e-1.5*(u_x_1*u_x_1+u_y_1*u_y_1)/speed2/e) 
					Bf(0)=-(Pr/Da*u_x_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_x_1)*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-(Pr/Da*u_y_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_y_1)*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)               
					h(0)=(speed*v_x(k)-u_x_1)/cs2+p*speed*v_x(k)/cs2/cs2
					h(1)=(speed*v_y(k)-u_y_1)/cs2+p*speed*v_y(k)/cs2/cs2						
                    Gf=delta_t*(1.d0-1.d0/2.d0/tau)*wi(k)*(bf(0)*h(0)+bf(1)*h(1))
                    f(x,y,k) = f_help(x,y,k) + ( f_eq - f_help(x,y,k) ) / tau + Gf
                end do
				
				case(3)!-------------semi porous-------------------
				
				if ( (mixctr(x,y)==1)) then
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)


					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))

				!m1=(T_loc*30+6) /21.0
				!m2=m1**-1.0868
                !tau_1=3*Pr_e*delta_t/delta_x/delta_x
				!tau_loc=tau_1*m2+0.5 
				
		
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x_1+speed*v_y(k)*u_y_1
                    f_eq=wi(k)*density_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4/e-1.5*(u_x_1*u_x_1+u_y_1*u_y_1)/speed2/e) 
					Bf(0)=-(Pr/Da*u_x_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_x_1)*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-(Pr/Da*u_y_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_y_1)*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)               
					h(0)=(speed*v_x(k)-u_x_1)/cs2+p*speed*v_x(k)/cs2/cs2
					h(1)=(speed*v_y(k)-u_y_1)/cs2+p*speed*v_y(k)/cs2/cs2						
                    Gf=delta_t*(1.d0-1.d0/2.d0/tau)*wi(k)*(bf(0)*h(0)+bf(1)*h(1))
                    f(x,y,k) = f_help(x,y,k) + ( f_eq - f_help(x,y,k) ) / tau + Gf
                end do

				else


				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
			!	m1=(T_loc*30+6) /21.0
			!	m2=m1**-1.0868
             !   tau_1=3*Pr_e*delta_t/delta_x/delta_x
			!	tau_loc=tau_1*m2+0.5 
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x+speed*v_y(k)*u_y
                    f_eq=wi(k)*density_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4-1.5*(u_x*u_x+u_y*u_y)/speed2)                  
					h(0)=(speed*v_x(k)-u_x)/cs2+p*speed*v_x(k)/cs2/cs2
					h(1)=(speed*v_y(k)-u_y)/cs2+p*speed*v_y(k)/cs2/cs2						
                    Gf=delta_t*(1.d0-1.d0/2.d0/tau)*wi(k)*(bf(0)*h(0)+bf(1)*h(1))
                    f(x,y,k) = f_help(x,y,k) + ( f_eq - f_help(x,y,k) ) / tau + Gf
                end do

            end if

		   end select

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    !============================== thermal
    !$OMP PARALLEL
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(x,y,k,kk,Density_loc,u_x,u_y,T_loc,f_eq,g_eq,p,BF,GF,h,u_x_sum,u_y_sum,u_x_0,u_y_0,u_x_1,u_y_1,tau_loc)
        do x=0,nx-1
        do y=0,ny-1
				!-------------------------------- local calculations	
				select case(bf_ctr)
				case(1)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x+speed*v_y(k)*u_y
                    g_eq=wi(k)*T_loc*(1.d0+3.d0*p/speed2+4.5*p*p/speed4-1.5*(u_x*u_x+u_y*u_y)/speed2)                    
                    g(x,y,k) = g_help(x,y,k) + ( g_eq - g_help(x,y,k) ) / tau_g
                end do

				case(2)
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))

				do k=0,8
                    p=speed*v_x(k)*u_x_1+speed*v_y(k)*u_y_1
                    g_eq=wi(k)*T_loc*(1.d0+3.d0*p/speed2)                    
                    g(x,y,k) = g_help(x,y,k) + ( g_eq - g_help(x,y,k) ) / tau_g
                end do
				case(3)!-------------semi porous-------------------
				if ( (mixctr(x,y)==1)) then
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)


					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation		
					
				u_x_0=(u_x_sum+delta_t*bf(0)*density_loc*e/2.d0)/density_loc	! initial approximation	
				u_y_0=(u_y_sum+delta_t*bf(1)*density_loc*e/2.d0)/density_loc	! initial approximation	
				
				u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
                u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))
				else
				density_loc=0.d0;T_loc=0.d0;u_x_sum=0.d0;u_y_sum=0.d0
				do k=0,8
					density_loc=density_loc+f_help(x,y,k)
					T_loc=T_loc+g_help(x,y,k)
					u_x_sum=u_x_sum+f_help(x,y,k)*v_x(k)*speed
					u_y_sum=u_y_sum+f_help(x,y,k)*v_y(k)*speed                    
				end do
				Bf(0)=+Ra*Pr*T_loc*sin(teta)   ! initial approximation							
				Bf(1)=+Ra*Pr*T_loc*cos(teta)   ! initial approximation			
				u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc	! initial approximation	
				u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc	! initial approximation	
				do kk=1,10
					Bf(0)=-Pr/Da*u_x*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
					Bf(1)=-Pr/Da*u_y*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	
					u_x_1=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					u_y_1=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
				end do
				end if
				
				!--------------------------------------------------------
                do k=0,8
                    p=speed*v_x(k)*u_x_1+speed*v_y(k)*u_y_1
                    g_eq=wi(k)*T_loc*(1.d0+3.d0*p/speed2)                    
                    g(x,y,k) = g_help(x,y,k) + ( g_eq - g_help(x,y,k) ) / tau_g
                end do
				end select
        end do
        end do
    !$OMP END DO
    !$OMP END PARALLEL
    end subroutine