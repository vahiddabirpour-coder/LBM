 subroutine temprature_bc
    use CaviData
    implicit none
    integer x,y,k
    real (8) Th,Tc,Tup,Tdown
    	!------------------------------------------left wall heat transfer, Th
    Th=1.0
	x=0
	do y=0,ny-1
		g_help(x,y,1)=4.d0/6.d0*(Th-(g_help(x,y,0)+g_help(x,y,2)+g_help(x,y,3)+g_help(x,y,4)+g_help(x,y,6)+g_help(x,y,7)))
		g_help(x,y,5)=1.d0/6.d0*(Th-(g_help(x,y,0)+g_help(x,y,2)+g_help(x,y,3)+g_help(x,y,4)+g_help(x,y,6)+g_help(x,y,7)))
		g_help(x,y,8)=1.d0/6.d0*(Th-(g_help(x,y,0)+g_help(x,y,2)+g_help(x,y,3)+g_help(x,y,4)+g_help(x,y,6)+g_help(x,y,7)))
	end do

        !-----------------------------------------right wall heat transfer, Tc
        Tc=0.0
        x=nx-1
	do y=0,ny-1
		g_help(x,y,3)=4.d0/6.d0*(Tc-(g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,2)+g_help(x,y,4)+g_help(x,y,5)+g_help(x,y,8)))
		g_help(x,y,7)=1.d0/6.d0*(Tc-(g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,2)+g_help(x,y,4)+g_help(x,y,5)+g_help(x,y,8)))
		g_help(x,y,6)=1.d0/6.d0*(Tc-(g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,2)+g_help(x,y,4)+g_help(x,y,5)+g_help(x,y,8)))
	end do
	!------------------------------------------upper wall, adiabatic
         
        y=ny-1
	do x=1,nx-2
		Tup=0
		  do k=0,8
			Tup=Tup+g_help(x,y-1,k)
		  end do
		g_help(x,y,4)=4.d0/6.d0*(Tup-(g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,2)+g_help(x,y,3)+g_help(x,y,5)+g_help(x,y,6)))
		g_help(x,y,7)=1.d0/6.d0*(Tup-(g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,2)+g_help(x,y,3)+g_help(x,y,5)+g_help(x,y,6)))
		g_help(x,y,8)=1.d0/6.d0*(Tup-(g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,2)+g_help(x,y,3)+g_help(x,y,5)+g_help(x,y,6)))
 
	end do	
         !-----------------------------------------lower wall,adiabatic

       y=0
	do x=1,nx-2
        Tdown=0
		     do k=0,8
			  Tdown=Tdown+g_help(x,y+1,k)
			  end do       
        g_help(x,y,2)=4.d0/6.d0*(Tdown-(g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,4)+g_help(x,y,3)+g_help(x,y,7)+g_help(x,y,8)))
		g_help(x,y,5)=1.d0/6.d0*(Tdown-(g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,4)+g_help(x,y,3)+g_help(x,y,7)+g_help(x,y,8)))
		g_help(x,y,6)=1.d0/6.d0*(Tdown-(g_help(x,y,0)+g_help(x,y,1)+g_help(x,y,4)+g_help(x,y,3)+g_help(x,y,7)+g_help(x,y,8)))
    
          
      end do
    end subroutine
    !_______________________________________________________________________________________________________________
    
    subroutine fluid_bc
    use CaviData
    implicit none
    integer x,y,k
    real (8) ux,uy,rho,t_loc,bf(0:1) 
	! Upper wall
	ux=(0.d0)/speed
	uy=(0.d0)/speed
	y=ny-1
	do x=0,nx-1
		!-------------------------------- local calculations	
        T_loc=0.d0
        do k=0,8
            T_loc=T_loc+g_help(x,y-1,k)                 
        end do
        Bf(0)=-Pr/Da*ux*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
        Bf(1)=-Pr/Da*uy*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)			
		!----------------------------------------------------		
		!rho=(f_help(x,y,0)+f_help(x,y,1)+f_help(x,y,3)+2.d0*(f_help(x,y,2)+f_help(x,y,5)+f_help(x,y,6))+delta_t*bf(1)/2.d0/speed)
		f_help(x,y,4)=f_help(x,y,2)
		f_help(x,y,7)=f_help(x,y,5)+0.5*(f_help(x,y,1)-f_help(x,y,3))+delta_t*(bf(0)+bf(1))/4.d0/speed
		f_help(x,y,8)=f_help(x,y,6)-0.5*(f_help(x,y,1)-f_help(x,y,3))-delta_t*(bf(0)-bf(1))/4.d0/speed		
	end do
	! Lower wall
	ux=0.d0/speed
	uy=0.d0/speed
	y=0
	do x=0,nx-1
		!-------------------------------- local calculations	
        T_loc=0.d0
        do k=0,8
            T_loc=T_loc+g_help(x,y+1,k)                 
        end do
        Bf(0)=-Pr/Da*ux*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
        Bf(1)=-Pr/Da*uy*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)			
		!----------------------------------------------------	
		!rho=(f_help(x,y,0)+f_help(x,y,1)+f_help(x,y,3)+2.d0*(f_help(x,y,4)+f_help(x,y,7)+f_help(x,y,8))-delta_t*bf(1)/2.d0/speed)
		f_help(x,y,2)=f_help(x,y,4)
		f_help(x,y,5)=f_help(x,y,7)-0.5*(f_help(x,y,1)-f_help(x,y,3))-delta_t*(bf(0)+bf(1))/4.d0/speed
		f_help(x,y,6)=f_help(x,y,8)+0.5*(f_help(x,y,1)-f_help(x,y,3))+delta_t*(bf(0)-bf(1))/4.d0/speed
	end do      
	! Right wall
	ux=0.d0/speed
	uy=0.d0/speed
	x=nx-1
	do y=0,ny-1
		!-------------------------------- local calculations	
        T_loc=0.d0
        Bf(0)=-Pr/Da*ux*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
        Bf(1)=-Pr/Da*uy*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)			
		!----------------------------------------------------
		!rho=(f_help(x,y,0)+f_help(x,y,2)+f_help(x,y,4)+2.d0*(f_help(x,y,1)+f_help(x,y,5)+f_help(x,y,8))+delta_t*bf(0)/2.d0/speed)
		f_help(x,y,3)=f_help(x,y,1)
		f_help(x,y,7)=f_help(x,y,5)+0.5*(f_help(x,y,2)-f_help(x,y,4))+delta_t*(bf(0)+bf(1))/4.d0/speed
		f_help(x,y,6)=f_help(x,y,8)-0.5*(f_help(x,y,2)-f_help(x,y,4))+delta_t*(bf(0)-bf(1))/4.d0/speed
	end do  	
	! Left wall
	ux=0.d0/speed
	uy=0.d0/speed
	x=0
	do y=0,ny-1
		!-------------------------------- local calculations	
        T_loc=T0
        Bf(0)=-Pr/Da*ux*mixctr(x,y)+Ra*Pr*T_loc*sin(teta)							
        Bf(1)=-Pr/Da*uy*mixctr(x,y)+Ra*Pr*T_loc*cos(teta)			
		!----------------------------------------------------
		!rho=(f_help(x,y,0)+f_help(x,y,2)+f_help(x,y,4)+2.d0*(f_help(x,y,7)+f_help(x,y,3)+f_help(x,y,6))-delta_t*bf(0)/2.d0/speed)
		f_help(x,y,1)=f_help(x,y,3)
		f_help(x,y,5)=f_help(x,y,7)-0.5*(f_help(x,y,2)-f_help(x,y,4))-delta_t*(bf(0)+bf(1))/4.d0/speed
		f_help(x,y,8)=f_help(x,y,6)+0.5*(f_help(x,y,2)-f_help(x,y,4))-delta_t*(bf(0)-bf(1))/4.d0/speed
	end do 
    end subroutine
    