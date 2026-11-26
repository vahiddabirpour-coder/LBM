subroutine geo_maker
    use CaviData
    implicit none
    integer x,y,o1,o2
	integer blocknumber
    real (8) ofr
    ofr=0.1
    o1=int(real(ny-1)*(1.d0-ofr)/2.d0)
    o2=o1+int(ofr*real(ny-1))
    do x=0,nx-1
	    solid(x,0) = .true.
	    solid(x,ny-1) = .true.
    end do
    do y=0,ny-1
	    solid(0,y) = .true.
	    solid(nx-1,y) = .true.
    end do
    solid(0,0)=.true.
    solid(0,ny-1)=.true.
    solid(nx-1,0)=.true.
    solid(nx-1,ny-1)=.true.
    select case(bf_ctr)
    case(1)
        print*,'-->  BF is adjusted for pure natural convection'
        mixctr(:,:)=0.d0
    case(2)
        print*,'-->  BF is adjusted for natural convection in porous media'        
        mixctr(:,:)=1.d0
    case(3)
        print*,'-->  BF is adjusted for natural convection in bidisperse porous media'     
		print*,'---------------------------------------------------------------------'
		print*,'number of blocks : 1) 9   2) 16   3) 36   4) 64'
		read*,blocknumber
		print*,'insert macroscopic porosity (for example 0.5)'
		read*,e_mac
        print*,'insert midel lineX  (0 < LM < 0.5)'
        read*,Lm
        print*,'insert midel lineY (0 < LM_y < 0.5)'
        read*,Lm_y
		select case(blocknumber)
   case(1)
    N_mac=9.00
    !e_mac=0.9          !free space in porous media
    N_box=sqrt(N_mac)
    D_mac=sqrt((1-e_mac)/N_mac)
    LOD=1.0-Lm    ! long of macroscopic domain*
    deltax=(LOD-(sqrt(N_mac))*D_mac)/2/sqrt(N_mac)
    deltay=(1.0-(sqrt(N_mac))*D_mac)/2/sqrt(N_mac)
    Lx_box=(LOD/N_box)
    Ly_box=1.0/N_box
    free_space=((4*(delta*delta))+(4*(delta*D_mac)))*9
    porous_space=(D_mac*D_mac)*9
    block_size=D_mac+(2*delta)
    !-----------------------------------------------------------
    print*,'                                '
    print*,'-----macroscopic porosity------'
	print*,'number of block=',N_mac
    print*,'number of box=' ,N_box
    print*,'long of macroscopic domain=' ,LOD
    print*,'long of box(x,y)=',Lx_box,Ly_box
    print*,'macroscopic porosity',e_mac
    print*,'D_mac=',D_mac
    print*,'delta=',delta
    print*,'free_space=',free_space
    print*,'porous_space=',porous_space
    print*,'block size=',block_size
    print*,'-------------------------------'
    print*,'                               '
	!------------------A-------------------
	do x= (Lm+deltax)*int(nx-1),((Lx_box+Lm)-deltax)*int(nx-1)
	  do y=((2*Ly_box)+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
	    mixctr(x,y)=1.d0
            end do
        end do
     !---------------B---------------------
	 do x=((Lx_box+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
	 do y=(0.66666+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
	    mixctr(x,y)=1.d0
            end do
        end do
     !-------------------C------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(nx-1)-deltax*int(nx-1)
	 do y=(0.66666+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
	    mixctr(x,y)=1.d0
            end do
        end do

    !________________________________________________________
	!-------------------D-----------------------------------
	do x= (Lm+deltax)*int(nx-1),((Lx_box+Lm)-deltax)*int(nx-1)
	  do y=(0.333333+deltay)*int(ny-1),(0.66666-deltay)*int(ny-1)
	    mixctr(x,y)=1.d0
            end do
        end do
     !---------------E---------------------
	 do x=((Lx_box+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
	 do y=(0.333333+deltay)*int(ny-1),(0.66666-deltay)*int(ny-1)
	    mixctr(x,y)=1.d0
            end do
        end do
     !-------------------F------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(nx-1)-deltax*int(nx-1)
	 do y=(0.333333+deltay)*int(ny-1),(0.66666-deltay)*int(ny-1)
	    mixctr(x,y)=1.d0
            end do
        end do
	!_______________________________________________________
	!-------------------G-----------------------------------
	do x= (Lm+deltax)*int(nx-1),((Lx_box+Lm)-deltax)*int(nx-1)
	  do y=deltay*int(ny-1),(0.333333-deltay)*int(ny-1)
	    mixctr(x,y)=1.d0
            end do
        end do
     !---------------H---------------------
	 do x=((Lx_box+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
	 do y=deltay*int(ny-1),(0.333333-deltay)*int(ny-1)
	    mixctr(x,y)=1.d0
            end do
        end do
     !-------------------I------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(nx-1)-deltax*int(nx-1)
	 do y=deltay*int(ny-1),(0.333333-deltay)*int(ny-1)
	    mixctr(x,y)=1.d0
            end do
        end do
	!______________________________________________________________________________________________________________________________

   case(2)
    N_mac=16.00
    !e_mac=0.9          !free space in porous media
    N_box=sqrt(N_mac)
    D_mac=sqrt((1-e_mac)/N_mac)
    LOD=1.0-Lm    ! long of macroscopic domain*
    deltax=(LOD-(sqrt(N_mac))*D_mac)/2/sqrt(N_mac)
    deltay=(1.0-(sqrt(N_mac))*D_mac)/2/sqrt(N_mac)
    Lx_box=(LOD/N_box)
    Ly_box=1.0/N_box
    free_space=((4*(delta*delta))+(4*(delta*D_mac)))*9
    porous_space=(D_mac*D_mac)*9
    block_size=D_mac+(2*delta)
    !-----------------------------------------------------------
    print*,'                                '
    print*,'-----macroscopic porosity------'
	print*,'number of block=',N_mac
    print*,'macroscopic porosity',e_mac
    print*,'D_mac=',D_mac
    print*,'delta=',delta
    print*,'free_space=',free_space
    print*,'porous_space=',porous_space
    print*,'block size=',block_size
    print*,'-------------------------------'
    print*,'                               '
        !-----------------A------------------
        do x=(Lm+deltax)*int(nx-1),((Lx_box+Lm)-deltax)*int(nx-1)
            do y=(0.75+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
        end do
        !----------------B-------------------
         do x=((Lx_box+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
            do y=(0.75+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
         end do
         !---------------c-------------------
          do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(ny-1)
            do y=(0.75+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
          end do
          !--------------D--------------------
           do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(nx-1)-deltax*int(nx-1)
            do y=(0.75+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
           end do
           !_________________________________________________________
           !----------------------E--------------------------
            do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
               do y=(0.5+deltay)*int(ny-1),(0.75-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
            end do
            !-----------------------F------------------------
             do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
               do y=(0.5+deltay)*int(ny-1),(0.75-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
             end do
             !-------------------G---------------------------
             do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(ny-1)
               do y=(0.5+deltay)*int(ny-1),(0.75-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
             end do
             !---------------H---------------------------
             do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(nx-1)-deltax*int(nx-1)
               do y=(0.5+deltay)*int(ny-1),(0.75-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
             end do
             
            !_____________________________________________________________
             !---------------I---------------------------
             do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
               do y=(0.25+deltay)*int(ny-1),(0.5-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
             end do
             !---------------J--------------------------
             do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
               do y=(0.25+deltay)*int(ny-1),(0.5-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
             end do
             !---------------K---------------------------
             do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(ny-1)
               do y=(0.25+deltay)*int(ny-1),(0.5-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
             end do
             !--------------L----------------------------
             do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(nx-1)-deltax*int(nx-1)
               do y=(0.25+deltay)*int(ny-1),(0.5-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
             end do
             !___________________________________________________________
             !------------M------------------------------
             do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
               do y=deltay*int(ny-1),(0.25-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
             end do
             !----------N-------------------------------
              do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
               do y=deltay*int(ny-1),(0.25-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
              end do
              !--------------O------------------------
              do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(ny-1)
               do y=deltay*int(ny-1),(0.25-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
              end do
              !-------------P------------------------
             do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(nx-1)-deltax*int(nx-1)
               do y=deltay*int(ny-1),(0.25-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
               end do
             end do
             !__________________________________________________________________________
			  case(3)
    N_mac=36.00
    !e_mac=0.9          !free space in porous media
    N_box=sqrt(N_mac)
    D_mac=sqrt((1-e_mac)/N_mac)
    LOD=1.0-Lm    ! long of macroscopic domain*
    deltax=(LOD-(sqrt(N_mac))*D_mac)/2/sqrt(N_mac)
    deltay=(1.0-(sqrt(N_mac))*D_mac)/2/sqrt(N_mac)
    Lx_box=(LOD/N_box)
    Ly_box=1.0/N_box
    free_space=((4*(delta*delta))+(4*(delta*D_mac)))*9
    porous_space=(D_mac*D_mac)*9
    block_size=D_mac+(2*delta)
    !-----------------------------------------------------------
    print*,'                                '
    print*,'-----macroscopic porosity------'
	print*,'number of block=',N_mac
    print*,'macroscopic porosity',e_mac
    print*,'D_mac=',D_mac
    print*,'delta=',delta
    print*,'free_space=',free_space
    print*,'porous_space=',porous_space
    print*,'block size=',block_size
    print*,'-------------------------------'
    print*,'                               '
	!______________________1_____________________________________
	  !-------------------1-------------------------
	  do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.166666+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------2--------------------------------
	   do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.166666+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------3------------------------------------
	    do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.166666+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !------------4-----------------------------------
	   do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.166666+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------5-----------------------------------
	   do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.166666+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------6---------------------------------------
	   do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.166666+deltay)*int(ny-1),(ny-1)-deltay*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !____________________2_________________________________
	   !-------------------1-------------------------
	  do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.166666+deltay)*int(ny-1),(5*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------2--------------------------------
	   do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.166666+deltay)*int(ny-1),(5*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------3------------------------------------
	    do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.166666+deltay)*int(ny-1),(5*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !------------4-----------------------------------
	   do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.166666+deltay)*int(ny-1),(5*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------5-----------------------------------
	   do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.166666+deltay)*int(ny-1),(5*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------6---------------------------------------
	   do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.166666+deltay)*int(ny-1),(5*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
    !_________________________3________________________________
	 !-------------------1-------------------------
	  do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.166666+deltay)*int(ny-1),(4*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------2--------------------------------
	   do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.166666+deltay)*int(ny-1),(4*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------3------------------------------------
	    do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.166666+deltay)*int(ny-1),(4*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !------------4-----------------------------------
	   do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.166666+deltay)*int(ny-1),(4*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------5-----------------------------------
	   do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.166666+deltay)*int(ny-1),(4*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------6---------------------------------------
	   do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.166666+deltay)*int(ny-1),(4*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !_______________________4________________________________

	   !-------------------1-------------------------
	  do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.166666+deltay)*int(ny-1),(3*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------2--------------------------------
	   do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.166666+deltay)*int(ny-1),(3*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------3------------------------------------
	    do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.166666+deltay)*int(ny-1),(3*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !------------4-----------------------------------
	   do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.166666+deltay)*int(ny-1),(3*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------5-----------------------------------
	   do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.166666+deltay)*int(ny-1),(3*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------6---------------------------------------
	   do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.166666+deltay)*int(ny-1),(3*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !_______________________________5________________________
	   !-------------------1-------------------------
	  do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.166666+deltay)*int(ny-1),(2*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------2--------------------------------
	   do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.166666+deltay)*int(ny-1),(2*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------3------------------------------------
	    do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.166666+deltay)*int(ny-1),(2*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !------------4-----------------------------------
	   do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.166666+deltay)*int(ny-1),(2*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------5-----------------------------------
	   do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.166666+deltay)*int(ny-1),(2*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------6---------------------------------------
	   do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.166666+deltay)*int(ny-1),(2*0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !__________________________6_____________________________
	   !-------------------1-------------------------
	  do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=deltay*int(ny-1),(0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------2--------------------------------
	   do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=deltay*int(ny-1),(0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------3------------------------------------
	    do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=deltay*int(ny-1),(0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !------------4-----------------------------------
	   do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=deltay*int(ny-1),(0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !----------------5-----------------------------------
	   do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=deltay*int(ny-1),(0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !--------------6---------------------------------------
	   do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=deltay*int(ny-1),(0.166666-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
            end do
       end do
	   !_________________________________________________________
	   case(4)
	   N_mac=64.00
    !e_mac=0.9          !free space in porous media
     N_box=sqrt(N_mac)
    D_mac=sqrt((1-e_mac)/N_mac)
    LOD=1.0-Lm    ! long of macroscopic domain*
    deltax=(LOD-(sqrt(N_mac))*D_mac)/2/sqrt(N_mac)
    deltay=(1.0-(sqrt(N_mac))*D_mac)/2/sqrt(N_mac)
    Lx_box=(LOD/N_box)
    Ly_box=1.0/N_box
    free_space=((4*(delta*delta))+(4*(delta*D_mac)))*9
    porous_space=(D_mac*D_mac)*9
    block_size=D_mac+(2*delta)
    !-----------------------------------------------------------
    print*,'                                '
    print*,'-----macroscopic porosity------'
	print*,'number of block=',N_mac
    print*,'macroscopic porosity',e_mac
    print*,'D_mac=',D_mac
    print*,'delta=',delta
    print*,'free_space=',free_space
    print*,'porous_space=',porous_space
    print*,'block size=',block_size
    print*,'-------------------------------'
    print*,'                               '
	!___________________________________________________________
	!____________________________1______________________________
	!---------------------------1-------------------------------
	do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(7*0.125+deltay)*int(ny-1),(8*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------2----------------------------------
	 do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(7*0.125+deltay)*int(ny-1),(8*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-----------------------3------------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(7*0.125+deltay)*int(ny-1),(8*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------4-----------------------------------
	 do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(7*0.125+deltay)*int(ny-1),(8*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !----------------------5-------------------------------------
	 do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(7*0.125+deltay)*int(ny-1),(8*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !---------------------6--------------------------------------
	 do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(7*0.125+deltay)*int(ny-1),(8*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !--------------------7---------------------------------------
	 do x=(((6*Lx_box)+Lm)+deltax)*int(nx-1),(((7*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(7*0.125+deltay)*int(ny-1),(8*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-------------------8----------------------------------------
	 do x=(((7*Lx_box)+Lm)+deltax)*int(nx-1),(((8*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(7*0.125+deltay)*int(ny-1),(8*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do

    !______________________________2____________________________
	!---------------------------1-------------------------------
	do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(6*0.125+deltay)*int(ny-1),(7*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------2----------------------------------
	 do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(6*0.125+deltay)*int(ny-1),(7*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-----------------------3------------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(6*0.125+deltay)*int(ny-1),(7*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------4-----------------------------------
	 do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(6*0.125+deltay)*int(ny-1),(7*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !----------------------5-------------------------------------
	 do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(6*0.125+deltay)*int(ny-1),(7*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !---------------------6--------------------------------------
	 do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(6*0.125+deltay)*int(ny-1),(7*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !--------------------7---------------------------------------
	 do x=(((6*Lx_box)+Lm)+deltax)*int(nx-1),(((7*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(6*0.125+deltay)*int(ny-1),(7*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-------------------8----------------------------------------
	 do x=(((7*Lx_box)+Lm)+deltax)*int(nx-1),(((8*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(6*0.125+deltay)*int(ny-1),(7*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !__________________________3_________________________________
	 !---------------------------1-------------------------------
	do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.125+deltay)*int(ny-1),(6*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------2----------------------------------
	 do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.125+deltay)*int(ny-1),(6*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-----------------------3------------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.125+deltay)*int(ny-1),(6*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------4-----------------------------------
	 do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.125+deltay)*int(ny-1),(6*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !----------------------5-------------------------------------
	 do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.125+deltay)*int(ny-1),(6*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !---------------------6--------------------------------------
	 do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.125+deltay)*int(ny-1),(6*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !--------------------7---------------------------------------
	 do x=(((6*Lx_box)+Lm)+deltax)*int(nx-1),(((7*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.125+deltay)*int(ny-1),(6*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-------------------8----------------------------------------
	 do x=(((7*Lx_box)+Lm)+deltax)*int(nx-1),(((8*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(5*0.125+deltay)*int(ny-1),(6*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !____________________________4_______________________________
	 !---------------------------1-------------------------------
	do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.125+deltay)*int(ny-1),(5*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------2----------------------------------
	 do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.125+deltay)*int(ny-1),(5*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-----------------------3------------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.125+deltay)*int(ny-1),(5*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------4-----------------------------------
	 do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.125+deltay)*int(ny-1),(5*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !----------------------5-------------------------------------
	 do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.125+deltay)*int(ny-1),(5*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !---------------------6--------------------------------------
	 do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.125+deltay)*int(ny-1),(5*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !--------------------7---------------------------------------
	 do x=(((6*Lx_box)+Lm)+deltax)*int(nx-1),(((7*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.125+deltay)*int(ny-1),(5*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-------------------8----------------------------------------
	 do x=(((7*Lx_box)+Lm)+deltax)*int(nx-1),(((8*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(4*0.125+deltay)*int(ny-1),(5*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !__________________________5_______________________________
	 !---------------------------1-------------------------------
	do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.125+deltay)*int(ny-1),(4*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------2----------------------------------
	 do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.125+deltay)*int(ny-1),(4*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-----------------------3------------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.125+deltay)*int(ny-1),(4*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------4-----------------------------------
	 do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.125+deltay)*int(ny-1),(4*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !----------------------5-------------------------------------
	 do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.125+deltay)*int(ny-1),(4*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !---------------------6--------------------------------------
	 do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.125+deltay)*int(ny-1),(4*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !--------------------7---------------------------------------
	 do x=(((6*Lx_box)+Lm)+deltax)*int(nx-1),(((7*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.125+deltay)*int(ny-1),(4*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-------------------8----------------------------------------
	 do x=(((7*Lx_box)+Lm)+deltax)*int(nx-1),(((8*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(3*0.125+deltay)*int(ny-1),(4*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !_____________________6______________________________________
	 !---------------------------1-------------------------------
	do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.125+deltay)*int(ny-1),(3*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------2----------------------------------
	 do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.125+deltay)*int(ny-1),(3*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-----------------------3------------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.125+deltay)*int(ny-1),(3*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------4-----------------------------------
	 do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.125+deltay)*int(ny-1),(3*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !----------------------5-------------------------------------
	 do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.125+deltay)*int(ny-1),(3*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !---------------------6--------------------------------------
	 do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.125+deltay)*int(ny-1),(3*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !--------------------7---------------------------------------
	 do x=(((6*Lx_box)+Lm)+deltax)*int(nx-1),(((7*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.125+deltay)*int(ny-1),(3*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-------------------8----------------------------------------
	 do x=(((7*Lx_box)+Lm)+deltax)*int(nx-1),(((8*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(2*0.125+deltay)*int(ny-1),(3*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !___________________________7______________________________
	  !---------------------------1-------------------------------
	do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.125+deltay)*int(ny-1),(2*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------2----------------------------------
	 do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.125+deltay)*int(ny-1),(2*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-----------------------3------------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.125+deltay)*int(ny-1),(2*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------4-----------------------------------
	 do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.125+deltay)*int(ny-1),(2*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !----------------------5-------------------------------------
	 do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.125+deltay)*int(ny-1),(2*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !---------------------6--------------------------------------
	 do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.125+deltay)*int(ny-1),(2*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !--------------------7---------------------------------------
	 do x=(((6*Lx_box)+Lm)+deltax)*int(nx-1),(((7*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.125+deltay)*int(ny-1),(2*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-------------------8----------------------------------------
	 do x=(((7*Lx_box)+Lm)+deltax)*int(nx-1),(((8*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(0.125+deltay)*int(ny-1),(2*0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !____________________8___________________________________
	 !---------------------------1-------------------------------
	do x=(Lm+deltax)*int(nx-1),(((1*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(deltay)*int(ny-1),(0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------2----------------------------------
	 do x=(((1*Lx_box)+Lm)+deltax)*int(nx-1),(((2*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(deltay)*int(ny-1),(0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-----------------------3------------------------------------
	 do x=(((2*Lx_box)+Lm)+deltax)*int(nx-1),(((3*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(deltay)*int(ny-1),(0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !------------------------4-----------------------------------
	 do x=(((3*Lx_box)+Lm)+deltax)*int(nx-1),(((4*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(deltay)*int(ny-1),(0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !----------------------5-------------------------------------
	 do x=(((4*Lx_box)+Lm)+deltax)*int(nx-1),(((5*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(deltay)*int(ny-1),(0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !---------------------6--------------------------------------
	 do x=(((5*Lx_box)+Lm)+deltax)*int(nx-1),(((6*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(deltay)*int(ny-1),(0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !--------------------7---------------------------------------
	 do x=(((6*Lx_box)+Lm)+deltax)*int(nx-1),(((7*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(deltay)*int(ny-1),(0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
     end do
	 !-------------------8----------------------------------------
	 do x=(((7*Lx_box)+Lm)+deltax)*int(nx-1),(((8*Lx_box)+Lm)-deltax)*int(nx-1)
           do y=(deltay)*int(ny-1),(0.125-deltay)*int(ny-1)
                mixctr(x,y)=1.d0
           end do
	 end do

    !___________________________________________________________        
             end select
         
    end select
    print*;
    end subroutine
	!-----------------------------------------------------------------------------

	