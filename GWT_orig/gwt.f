***********************************************************************
*                                                                     *
*                            >>  PROGRAM GWT  <<                      *
*                                                                     *
*---------------------------------------------------------------------*
*  purpose : main program of GWT                                      *
*---------------------------------------------------------------------*
*  subprograms    : read_dim, read_grid, read_matr, read_riv,         *
*                   read_hydro, read_rech, lin_equ, lin_solv,         *
*                   write_inter, write_heads                          *
*---------------------------------------------------------------------*
*  Last update    : 16-OCT-2010 by Toon Leijnse                       *
***********************************************************************
      program gwt
*----------------------------------------------------------------------
      implicit none
	double precision delt
	integer itim, nx, ny, nrbl, ngrid, ndim, neq
	logical active
	allocatable active(:)
	integer ieq, ixr, iyr
	allocatable ieq(:), ixr(:), iyr(:)
	double precision totar, cumbal, delx, dely, xc, yc
	dimension cumbal(7)
	allocatable delx(:), dely(:), xc(:), yc(:)
	double precision h, hold, bot, permx, permy, por, rech
	allocatable h(:), hold(:), bot(:), permx(:), permy(:), por(:), 
     &            rech(:)
	double precision rar, res, rl
      allocatable rar(:), res(:), rl(:)
	integer ir, ic
	double precision am, bm, sol, sc, p, ap, r
	allocatable ir(:), ic(:), am(:), bm(:), sol(:), sc(:), p(:), 
     &            ap(:), r(:)
*----------------------------------------------------------------------
	integer i, itmax, it, ios
	double precision eps
	logical fler
*----------------------------------------------------------------------
*     Open log file
*
      open (15,file='gwt_logfile',position='append')
*----------------------------------------------------------------------
*     Open timestep input file and read timestep number and 
*     timestepsize
*
      open (10, file='gwt_timestepinput', status='old',iostat=ios)
	if (ios.ne.0) then
	   write (15,'(a)') 'Error opening file gwt_timestepinput'
	   ngrid = 10
	   allocate (h(ngrid))
	   h = -999.d0
	   goto 10
	endif
	read (10,*,iostat=ios) itim, delt
	if (ios.ne.0) then
	   write (15,'(a,a)') 'Error reading timestep number from ',
     &                      'gwt_timestepinput'
	   goto 10
	endif
	cumbal = 0.d0
*----------------------------------------------------------------------
*     Read dimensioning data
*
      call read_dim (itim, nx, ny, nrbl, ngrid, ndim, neq, totar, 
     &               cumbal, fler)
	if (fler) goto 10
*----------------------------------------------------------------------
*     Allocate space
*
      allocate (active(ngrid), ieq(ngrid), ixr(max(nrbl,1)), 
     &          iyr(max(nrbl,1)))
      allocate (delx(nx), dely(ny), xc(nx), yc(ny))
	allocate (h(ngrid), hold(ngrid), bot(ngrid), permx(ngrid), 
     &          permy(ngrid), por(ngrid), rech(ngrid), 
     &          rar(max(nrbl,1)), res(max(nrbl,1)), rl(max(nrbl,1)))
	allocate (ir(neq+1), ic(ndim), am(ndim), bm(neq), sol(neq),
     &          sc(neq), p(neq), ap(neq), r(neq))
*----------------------------------------------------------------------
*     Read grid definition
*
      call read_grid (itim, nx, ny, ngrid, neq, delx, dely, xc, yc,
     &                active, ieq, totar, fler)
	if (fler) goto 10
*----------------------------------------------------------------------
*     Read or determine matrix structure
*
      call read_matr (itim, nx, ny, ngrid, active, ir, ic, ndim, neq,
     &                ieq, fler)
	if (fler) goto 10
*----------------------------------------------------------------------
*     Read river data
*
      if (nrbl.gt.0) call read_riv (itim, nx, ny, nrbl, ixr, iyr, rar, 
     &                              res, fler)
	if (fler) goto 10
*----------------------------------------------------------------------
*     Read hydrological data
*
      call read_hydro (itim, ngrid, h, bot, permx, permy, por, fler)
	if (fler) goto 10
	hold = h
*----------------------------------------------------------------------
*     Read recharge and river levels
*
      call read_rech (ngrid, nrbl, rech, rl, fler)
	if (fler) goto 10
*----------------------------------------------------------------------
*     Set up linear equations and solve
*      
      call lin_equ (nx, ny, ngrid, neq, ndim, am, bm, ir, ic, active,
     &              ieq, h, bot, delx, dely, xc, yc, permx, permy, 
     &              por, delt, rech, nrbl, ixr, iyr, rar, res, rl)
*
*     Convergence criteria and starting values for the solution
*
	eps = 1.d-10
	itmax = max(100,ngrid)
	do i = 1,ngrid
	   if (active(i)) sol(ieq(i)) = h(i)
	enddo
	call lin_solv (neq, ndim, am, bm, sol, ir, ic, eps, itmax, it, 
     &               sc, r, p, ap)
	fler = it.gt.itmax
*----------------------------------------------------------------------
*     Write results
* 
   10 h = -999.d0
      write (15,'(50(''*''))')
	write (15,'(/,a,i5,10x,a,1pe12.4)') 'Timestep number =',itim,
     &            'Timestep size = ',delt
      if (.not.fler) then
         do i = 1,ngrid
	      if (active(i)) h(i) = sol(ieq(i))
	   enddo
	   write (15,'(10x,a,i5,a)') 'Linear equations converged after',
     &                            it,' iterations'
*----------------------------------------------------------------------
*     Write water balance
*
         call write_balance (nx, ny, nrbl, delx, dely, active, ixr,
     &                       iyr, rar, res, rl, h, hold, rech, por,
     &                       delt, totar, cumbal, ngrid) 
*----------------------------------------------------------------------
*     Write intermediate file
* 
         call write_inter (nx, ny, nrbl, ndim, neq, totar, delx, dely, 
     &                     xc, yc, active, ieq, ir, ic, ixr, iyr, rar,
     &                     res, h, bot, permx, permy, por, cumbal, 
     &                     ngrid)
      else
	    if (it.gt.itmax) write (15,'(10x,a)') 
     &                          'No convergence for linear equations'
      endif
*----------------------------------------------------------------------
*     Write heads to output file
*
      call write_heads (ngrid,h) 
*----------------------------------------------------------------------
	stop
      end