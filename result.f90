 subroutine NU_calc
    use CaviData
    implicit none
    integer x,y,i,y_start,y_stop
    real (8) T1,T2,T3,deff,nu_local(0:ny-1)
    x=0
    y_start=0
    y_stop=ny-1
	do y=y_start,y_stop
	     T1=T0
         T2=0.;T3=0.
         do i=0,8
              T2=T2+g_help(x+1,y,i)
			  T3=T3+g_help(x+2,y,i)
         end do
         deff=(T2-T1)/delta_x                 ! first order
         !deff=(4.0*T2-T3-3.0*T1)/2.0/delta_x  ! secind order
	     nu_local(y)=-deff/(real(nx-1)/real(ny-1))
    end do
    nu_local(0)=nu_local(1)
    nu_local(ny-1)=nu_local(ny-2)
	ave_nuselt=0.
	do y=y_start+1,y_stop-1
		ave_nuselt=ave_nuselt+nu_local(y)*delta_x
	end do
	ave_nuselt=ave_nuselt+0.5*delta_x*(nu_local(y_start)+nu_local(y_stop))
    Print*;print*,' -------------------------------       Nu_ave: ',ave_nuselt;    Print*;
    if (ave_nuselt>0) then
        nu_change=abs(ave_nuselt-nu_his)/nu_his
        print*,'-Nu_ave change in % ',nu_change*100.0
        nu_his=ave_nuselt
    end if
    end subroutine
	!-------------------------------------------------------------------------
	 subroutine  write_results
     use CaviData
      implicit none   
      integer  x,y,i,obsval,k,kk
      real (8) density_loc,u_x,u_y,T_loc,press,bf(0:1),u_x_sum,u_y_sum
    character (len=300) :: filename1,pr_str,ra_str,Da_str,nx_str
     write (pr_str, "(A3,F4.2)") "_Pr", pr
     write (Ra_str, "(A3,F8.1)") "_Ra", Ra
      write (Da_str, "(A3,F8.6)") "_Da", Da
      write (Nx_str, "(A3,I3)") "_Nx", nx
      filename1='Tecplot_'// trim(adjustl(Nx_str))// trim(adjustl(Pr_str)) //  trim(adjustl(Da_str)) // trim(adjustl(Ra_str))//'.DAT'
      open(11,file=filename1)
      !open(11,file='Data.dat')
      write(11,*) 'VARIABLES = X, Y, Ux, Uy, PRESS, Temprature, Porous, OBST' 
      write(11,*) 'ZONE I=', nx, ', J=', ny, ', F=POINT'
      do y = 0, ny-1
        do x = 0, nx-1
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
          if (solid(x,y)) then 
            obsval = 1
            press = density * cs2
          else
            press = density_loc * Cs2
            obsval = 0 
          end if
          write(11,*) x, y, u_x, u_y, press,T_loc, mixctr(x,y), obsval
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
					
			u_x_0=(u_x_sum+delta_t*(Bf(0))*density_loc*e/2.d0)/density_loc
            u_y_0=(u_y_sum+delta_t*(bf(1))*density_loc*e/2.d0)/density_loc

            u_x_1=u_x_0/(C_0+sqrt(C_0**2+C_1*abs(u_x_0)))
            u_y_1=u_y_0/(C_0+sqrt(C_0**2+C_1*abs(u_y_0)))

			!--------------------------------------------------------            
          if (solid(x,y)) then 
            obsval = 1
            press = density * cs2
          else
            press = density_loc * Cs2
            obsval = 0 
          end if
          write(11,*) x, y, u_x_1, u_y_1, press,T_loc, mixctr(x,y), obsval
          
            case(3)!------------------bidisperse porous----------------
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
          if (solid(x,y)) then 
            obsval = 1
            press = density * cs2
          else
            press = density_loc * Cs2
            obsval = 0 
          end if
          write(11,*) x, y, u_x_1, u_y_1, press,T_loc, mixctr(x,y), obsval
          
		  	  end select
		end do
	  end do


      close(11)
      return
    end subroutine

	!--------------------------------------------------------------------------
	subroutine  write_results_fr(fr)
     use CaviData
      implicit none   
      integer  x,y,i,obsval,k,kk
      real (8) density_loc,u_x,u_y,T_loc,press,fr,bf(0:1),u_x_sum,u_y_sum
    character (len=300) :: filename1,pr_str,ra_str,Da_str,nx_str,fr_str
     write (pr_str, "(A3,F4.2)") "_Pr", pr
    write (pr_str, "(A3,F5.3)") "_Fr", fr
     write (Ra_str, "(A3,F8.1)") "_Ra", Ra
      write (Da_str, "(A3,F8.6)") "_Da", Da
      write (Nx_str, "(A3,I3)") "_Nx", nx
      filename1='Tecplot_'// trim(adjustl(Nx_str))// trim(adjustl(fr_str))// trim(adjustl(Pr_str)) //  trim(adjustl(Da_str)) // trim(adjustl(Ra_str))//'.DAT'
      open(11,file=filename1)
      !open(11,file='Data.dat')
      write(11,*) 'VARIABLES = X, Y, Ux, Uy, PRESS, Temprature, Porous, OBST' 
      write(11,*) 'ZONE I=', nx, ', J=', ny, ', F=POINT'
      do y = 0, ny-1
        do x = 0, nx-1
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
          if (solid(x,y)) then 
            obsval = 1
            press = density * cs2
          else
            press = density_loc * Cs2
            obsval = 0 
          end if
          write(11,*) x, y, u_x, u_y, press,T_loc, mixctr(x,y), obsval
		  !--------------------------------------------

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

				!do kk=1,10

				!	Bf(0)=-(Pr/Da*u_x_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_x_1)*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
				!	Bf(1)=-(Pr/Da*u_y_1)-(Cf/sqrt(Da)*sqrt(u_x_1**2+u_y_1**2)*u_y_1)*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)	

					!u_x=(u_x_sum+delta_t*bf(0)/2.d0)/density_loc
					!u_y=(u_y_sum+delta_t*bf(1)/2.d0)/density_loc
			!	end do
			!--------------------------------------------------------            
          if (solid(x,y)) then 
            obsval = 1
            press = density * cs2
          else
            press = density_loc * Cs2
            obsval = 0 
          end if
          write(11,*) x, y, u_x_1, u_y_1, press,T_loc, mixctr(x,y), obsval
		  !--------------------------------------------------------------------------
            case(3)   !----------------bidisperse------------
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
                if (solid(x,y)) then 
            obsval = 1
            press = density * cs2
          else
            press = density_loc * Cs2
            obsval = 0 
          end if
          write(11,*) x, y, u_x_1, u_y_1, press,T_loc, mixctr(x,y), obsval
          !----------------------------------------------------------
		  end select
		end do
	  end do
	  
      close(11)
      return
    end subroutine