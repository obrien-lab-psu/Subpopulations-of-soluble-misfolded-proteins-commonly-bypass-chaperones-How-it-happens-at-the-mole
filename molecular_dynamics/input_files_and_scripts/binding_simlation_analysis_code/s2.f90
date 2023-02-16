       program Ct_St_TCF

       implicit none
       real(kind=8),allocatable :: rprot(:,:),rw(:,:,:),dist(:),rprot1(:,:)
       real(kind=8),allocatable :: counter(:),correlation(:), con(:,:)
       real(kind=8):: boxl(3),boxl1(3),boxl2(3),boxl3(3),pi,dr(3),dr1(3),cutoff
       real(kind=8):: inter_oxy_dist,dist11,dist12,dist21,dist22, d
       real(kind=8):: temp_array(4),angle_conv,cond_3,ohdistsq,theta, costheta
       real(kind=8):: u1(3),u2(3),min_angle,oodistsq,min_dist,sum2, width = 2.0
       real(kind=8):: corr,sum1,count1,count2,dt=0.004d0,count3,p,q, hist(90)

       integer,allocatable :: gwi(:),rgwi(:)
       integer(kind=1),allocatable :: location(:,:)
       integer,allocatable :: ht_val(:,:)
       integer,allocatable :: info(:,:,:)!,hb_info(:,:)
       integer :: gw, rgw,jmol_water,kmol_water,corr_frame,ip,j,k,mi, maxbin=90
       integer :: natom_prot,iatom_prot,jatom_prot,natom_eth,iatom_eth,natom_prot1
       integer :: nmol_eth,imol_eth,iframe,nframe,i,ii,pair
       integer :: hbpair,kk,ih,flag, bin, total_contact
      
       character(len=4),allocatable :: prot_atom(:),acceptor(:)
       character(len=4) :: junk2
       character(len=6) :: junk1
       character(len=100) :: prot_gro,sol_gro, prot1_gro,traj4_gro,traj6_gro
      
      traj4_gro='groel.gro'
      traj6_gro='pro.gro'
       open(unit=1,file=traj4_gro,action='read')
       open(unit=3,file=traj6_gro, action='read')
       !open(unit=3,file='1hew.gro',action='read')
       ! System Information
       ! Number of Protein Atoms (Lysozyme)
       natom_prot = 1309
       natom_prot1=384
       write(567,*)'Protein Atom',natom_prot
       ! Number of Water molecules
       !nmol_eth=4199;natom_eth=1
       !write(567,*)'Water Molecule',nmol_eth
       ! Total Number of Snapshot
       nframe = 19680
       write(568,*)'Total Frame',nframe
       ! Constants
       pi=4.0d0*atan(1.0d0)
       angle_conv=180.0d0/pi
       cutoff=4.25d0
       oodistsq=3.5d0*3.5d0
       ohdistsq=2.45d0*2.45d0
       cond_3=30.0d0
       ! Allocation
       allocate(rprot(natom_prot,3))
       allocate(rprot1(natom_prot1,3))
       !allocate(prot_atom(natom_prot))
       !allocate(dist(natom_prot))
       allocate(con(natom_prot,natom_prot1))
       !allocate(gwi(nmol_eth),rgwi(nmol_eth))
       !rw=0.0d0; rprot=0.0d0; count1=0.0d0; hist=0.0d0
       ! Reading the Protein atom name from 1hew.gro
       !read(3,*)
       !read(3,*)
       !do iatom_prot = 1 , natom_prot
        !read(3,*) junk1, prot_atom(iatom_prot)
       !end do
       ! ==========================
         do iatom_prot = 1 , natom_prot
         do jatom_prot = 1, natom_prot1
            con(iatom_prot,jatom_prot)= 0
         enddo
         enddo
       	do iframe = 1 , nframe
       !gwi=0;rgwi=0
        !gw=0
        ! Reading Protein Coordinates
        read(1,*)
        read(1,*)
        do iatom_prot = 1 , natom_prot
         read(1,*) junk1, junk2, i, rprot(iatom_prot,:)
        end do
        read(1,*)boxl1(:)
        ! Converting into \AA
        rprot=rprot*10.0d0;boxl1=boxl1*10.0d0
       ! write(21,*) rprot
        read(3,*)
        read(3,*)
        do jatom_prot = 1 , natom_prot1
        read(3,*) junk1, junk2, i, rprot1(jatom_prot,:)
        end do
        read(3,*)boxl3(:)
        ! Converting into \AA
        rprot1=rprot1*10.0d0;boxl3=boxl3*10.0d0
        write(111,*)rprot1
        ! Reading Water coordinates
       ! read(2,*)
        !read(2,*)
        !do imol_eth = 1 , nmol_eth
         !do iatom_eth = 1 , natom_eth
!          if(imol_water <= 3333) then
          ! read(2,*) junk1, junk2, i, rw(imol_eth,iatom_eth,:)
!          else
!           read(2,*) junk1, junk2,    rw(imol_water,iatom_water,:)
!          end if
       ! write(100,*) rw(imol_eth,iatom_eth,:)
        ! end do
        !end do
        !read(2,*)boxl2(:)
        ! Converting into \AA
        !rw=rw*10.0d0;boxl2=boxl2*10.d0
        !if(mod(iframe,1000)==0)then 
         !write(567,*)'Frame No.',iframe
         !write(567,*)boxl1
         !write(567,*)boxl2
        !end if
        ! Searching the water 
         !gw = 0
       ! gwi=0
        !do imol_eth = 1 , nmol_eth
         !dist=999.99d0
         do iatom_prot = 1 , natom_prot
          do jatom_prot = 1, natom_prot1
        dr(:)=rprot1(jatom_prot,:)-rprot(iatom_prot,:) 
       ! write(402,*) dr       
!dr(:)=rprot(iatom_prot,:)-rprot(iatom_prot,:)
        ! dr(:)= dr(:)-boxl1(:)*anint(dr(:)/boxl1(:))
        d = sqrt(dot_product(dr,dr))
          !write(401,*) d
         ! write(400,*) dr
          gw=0 ; rgw=0
      if (d <= 8.0) then
         con(iatom_prot,jatom_prot)=con(iatom_prot,jatom_prot) + 1  
          gw=gw+1 
          endif
          ! write(200,*) con(iatom_prot,jatom_prot)
           ! write(240,*) iatom_prot, jatom_prot, gw

           !write(300,*) gw
            end do
             enddo

        total_contact=0
        do iatom_prot = 1 , natom_prot
          do jatom_prot = 1, natom_prot1
        dr(:)=rprot1(jatom_prot,:)-rprot(iatom_prot,:)
       ! dr(:)= dr(:)-boxl1(:)*anint(dr(:)/boxl1(:))
        d = sqrt(dot_product(dr,dr))
       if (d <= 8.0) then
        total_contact=total_contact+1   
         endif
          end do
           enddo
      if (total_contact>=140) then
        write(1000,*) iframe, 1
        write(1100,*) iframe, 1
         else
         write(1000,*)iframe,0
          endif
         
            enddo
      do iatom_prot = 1 , natom_prot
          do jatom_prot = 1, natom_prot1
          con(iatom_prot,jatom_prot)=con(iatom_prot,jatom_prot)/nframe
      ! write(100,*) iatom_prot, jatom_prot, con(iatom_prot,jatom_prot)
       !write(101,*) con(iatom_prot,jatom_prot)
       enddo  
         enddo
       end program Ct_St_TCF  
