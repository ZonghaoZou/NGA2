! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   Refactored by Federico Perini, 10/6/2022
!   Based on the netlib library by Paul Dierckx
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************
module fitpack_core
    use iso_c_binding, only: c_double,c_int32_t,c_bool
    implicit none
    private

    ! Precision and array size
    integer, parameter, public :: FP_REAL = c_double
    integer, parameter, public :: FP_SIZE = c_int32_t
    integer, parameter, public :: FP_FLAG = c_int32_t
    integer, parameter, public :: FP_BOOL = c_bool

    ! Curve fitting routines
    public :: curfit ! * General curve fitting

    ! Curve approximation routines
    public :: splev  ! * Evaluation of a spline function

    ! Spline behavior for points not in the support
    integer(FP_FLAG), parameter,  public :: OUTSIDE_EXTRAPOLATE = 0 ! extrapolated from the end spans
    integer(FP_FLAG), parameter,  public :: OUTSIDE_ZERO        = 1 ! spline evaluates to zero
    integer(FP_FLAG), parameter,  public :: OUTSIDE_NOT_ALLOWED = 2 ! an error flag is returned
    integer(FP_FLAG), parameter,  public :: OUTSIDE_NEAREST_BND = 3 ! evaluate to value of nearest boundary point

    ! Spline degrees
    integer(FP_SIZE), parameter, public :: MAX_ORDER = 19    ! Max spline order (for array allocation)

    integer(FP_FLAG), parameter,  public :: FITPACK_OK                   = 0  ! ok for spline, abs(fp-s)/s <= tol=0.001
    integer(FP_FLAG), parameter,  public :: FITPACK_INTERPOLATING_OK     = -1 ! ok for interpolating spline, fp=0
    integer(FP_FLAG), parameter,  public :: FITPACK_LEASTSQUARES_OK      = -2 ! ok for weighted least-squares polynomial of degree k.
    integer(FP_FLAG), parameter,  public :: FITPACK_INSUFFICIENT_STORAGE = 1
    integer(FP_FLAG), parameter,  public :: FITPACK_S_TOO_SMALL          = 2
    integer(FP_FLAG), parameter,  public :: FITPACK_MAXIT                = 3
    integer(FP_FLAG), parameter,  public :: FITPACK_INVALID_RANGE        = 4
    integer(FP_FLAG), parameter,  public :: FITPACK_INPUT_ERROR          = 5

    ! Internal Parameters
    logical(FP_BOOL), parameter, public :: FP_FALSE = .false._FP_BOOL

    integer(FP_SIZE), parameter :: IONE   = 1_FP_SIZE

    real(FP_REAL), parameter, public :: one     = 1.0_FP_REAL
    real(FP_REAL), parameter, public :: zero    = 0.0_FP_REAL
    real(FP_REAL), parameter, public :: half    = 0.5_FP_REAL
    real(FP_REAL), parameter, public :: three   = 3.0_FP_REAL
    real(FP_REAL), parameter, public :: smallnum03 = 1.0e-03_FP_REAL

    contains
      pure subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)

      !  given the set of data points (x(i),y(i)) and the set of positive
      !  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
      !  approximation of degree k on the interval xb <= x <= xe.
      !  if iopt=-1 curfit calculates the weighted least-squares spline
      !  according to a given set of knots.
      !  if iopt>=0 the number of knots of the spline s(x) and the position
      !  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
      !  ness of s(x) is then achieved by minimalizing the discontinuity
      !  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
      !  n-k-1. the amount of smoothness is determined by the condition that
      !  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
      !  negative constant, called the smoothing factor.
      !  the fit s(x) is given in the b-spline representation (b-spline coef-
      !  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
      !  subroutine splev.
      !
      !  calling sequence:
      !     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
      !    * lwrk,iwrk,ier)
      !
      !  parameters:
      !   iopt  : integer flag. on entry iopt must specify whether a weighted
      !           least-squares spline (iopt=-1) or a smoothing spline (iopt=
      !           0 or 1) must be determined. if iopt=0 the routine will start
      !           with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
      !           k+1. if iopt=1 the routine will continue with the knots
      !           found at the last call of the routine.
      !           attention: a call with iopt=1 must always be immediately
      !           preceded by another call with iopt=1 or iopt=zero
      !           unchanged on exit.
      !   m     : integer. on entry m must specify the number of data points.
      !           m > k. unchanged on exit.
      !   x     : real array of dimension at least (m). before entry, x(i)
      !           must be set to the i-th value of the independent variable x,
      !           for i=1,2,...,m. these values must be supplied in strictly
      !           ascending order. unchanged on exit.
      !   y     : real array of dimension at least (m). before entry, y(i)
      !           must be set to the i-th value of the dependent variable y,
      !           for i=1,2,...,m. unchanged on exit.
      !   w     : real array of dimension at least (m). before entry, w(i)
      !           must be set to the i-th value in the set of weights. the
      !           w(i) must be strictly positive. unchanged on exit.
      !           see also further comments.
      !   xb,xe : real values. on entry xb and xe must specify the boundaries
      !           of the approximation interval. xb<=x(1), xe>=x(m).
      !           unchanged on exit.
      !   k     : integer. on entry k must specify the degree of the spline.
      !           1<=k<=5. it is recommended to use cubic splines (k=3).
      !           the user is strongly dissuaded from choosing k even,together
      !           with a small s-value. unchanged on exit.
      !   s     : real.on entry (in case iopt>=0) s must specify the smoothing
      !           factor. s >=zero unchanged on exit.
      !           for advice on the choice of s see further comments.
      !   nest  : integer. on entry nest must contain an over-estimate of the
      !           total number of knots of the spline returned, to indicate
      !           the storage space available to the routine. nest >=2*k+2.
      !           in most practical situation nest=m/2 will be sufficient.
      !           always large enough is  nest=m+k+1, the number of knots
      !           needed for interpolation (s=0). unchanged on exit.
      !   n     : integer.
      !           unless ier =10 (in case iopt >=0), n will contain the
      !           total number of knots of the spline approximation returned.
      !           if the computation mode iopt=1 is used this value of n
      !           should be left unchanged between subsequent calls.
      !           in case iopt=-1, the value of n must be specified on entry.
      !   t     : real array of dimension at least (nest).
      !           on successful exit, this array will contain the knots of the
      !           spline,i.e. the position of the interior knots t(k+2),t(k+3)
      !           ...,t(n-k-1) as well as the position of the additional knots
      !           t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
      !           the b-spline representation.
      !           if the computation mode iopt=1 is used, the values of t(1),
      !           t(2),...,t(n) should be left unchanged between subsequent
      !           calls. if the computation mode iopt=-1 is used, the values
      !           t(k+2),...,t(n-k-1) must be supplied by the user, before
      !           entry. see also the restrictions (ier=10).
      !   c     : real array of dimension at least (nest).
      !           on successful exit, this array will contain the coefficients
      !           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
      !   fp    : real. unless ier=10, fp contains the weighted sum of
      !           squared residuals of the spline approximation returned.
      !   wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
      !           used as working space. if the computation mode iopt=1 is
      !           used, the values wrk(1),...,wrk(n) should be left unchanged
      !           between subsequent calls.
      !   lwrk  : integer. on entry,lwrk must specify the actual dimension of
      !           the array wrk as declared in the calling (sub)program.lwrk
      !           must not be too small (see wrk). unchanged on exit.
      !   iwrk  : integer array of dimension at least (nest).
      !           used as working space. if the computation mode iopt=1 is
      !           used,the values iwrk(1),...,iwrk(n) should be left unchanged
      !           between subsequent calls.
      !   ier   : integer. unless the routine detects an error, ier contains a
      !           non-positive value on exit, i.e.
      !    ier=0  : normal return. the spline returned has a residual sum of
      !             squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !             ive tolerance set to 0.001 by the program.
      !    ier=-1 : normal return. the spline returned is an interpolating
      !             spline (fp=0).
      !    ier=-2 : normal return. the spline returned is the weighted least-
      !             squares polynomial of degree k. in this extreme case fp
      !             gives the upper bound fp0 for the smoothing factor s.
      !    ier=1  : error. the required storage space exceeds the available
      !             storage space, as specified by the parameter nest.
      !             probably causes : nest too small. if nest is already
      !             large (say nest > m/2), it may also indicate that s is
      !             too small
      !             the approximation returned is the weighted least-squares
      !             spline according to the knots t(1),t(2),...,t(n). (n=nest)
      !             the parameter fp gives the corresponding weighted sum of
      !             squared residuals (fp>s).
      !    ier=2  : error. a theoretically impossible result was found during
      !             the iteration process for finding a smoothing spline with
      !             fp = s. probably causes : s too small.
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=3  : error. the maximal number of iterations maxit (set to 20
      !             by the program) allowed for finding a smoothing spline
      !             with fp=s has been reached. probably causes : s too small
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=10 : error. on entry, the input data are controlled on validity
      !             the following restrictions must be satisfied.
      !             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
      !             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
      !             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
      !                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
      !                       the schoenberg-whitney conditions, i.e. there
      !                       must be a subset of data points xx(j) such that
      !                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
      !             if iopt>=0: s>=0
      !                         if s=0 : nest >= m+k+1
      !             if one of these conditions is found to be violated,control
      !             is immediately repassed to the calling program. in that
      !             case there is no approximation returned.
      !
      !  further comments:
      !   by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the spline will be too smooth and signal will be
      !   lost ; if s is too small the spline will pick up too much noise. in
      !   the extreme cases the program will return an interpolating spline if
      !   s=0 and the weighted least-squares polynomial of degree k if s is
      !   very large. between these extremes, a properly chosen s will result
      !   in a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the weights w(i). if these are
      !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
      !   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
      !   sqrt(2*m)). if nothing is known about the statistical error in y(i)
      !   each w(i) can be set equal to one and s determined by trial and
      !   error, taking account of the comments above. the best is then to
      !   start with a very large value of s ( to determine the least-squares
      !   polynomial and the corresponding upper bound fp0 for s) and then to
      !   progressively decrease the value of s ( say by a factor 10 in the
      !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
      !   approximation shows more detail) to obtain closer fits.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt=zero
      !   if iopt=1 the program will continue with the set of knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if curfit is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   function underlying the data. but, if the computation mode iopt=1
      !   is used, the knots returned may also depend on the s-values at
      !   previous calls (if these were smaller). therefore, if after a number
      !   of trials with different s-values and iopt=1, the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   curfit once more with the selected value for s but now with iopt=0.
      !   indeed, curfit may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !
      !  other subroutines required:
      !    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
      !
      !  references:
      !   dierckx p. : an algorithm for smoothing, differentiation and integ-
      !                ration of experimental data using spline functions,
      !                j.comp.appl.maths 1 (1975) 165-184.
      !   dierckx p. : a fast algorithm for smoothing data on a rectangular
      !                grid while using spline functions, siam j.numer.anal.
      !                19 (1982) 1286-1304.
      !   dierckx p. : an improved algorithm for curve fitting with spline
      !                functions, report tw54, dept. computer science,k.u.
      !                leuven, 1981.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : may 1979
      !  latest update : march 1987
      !
      !  ..
      !  ..scalar arguments..
      real(FP_REAL),    intent(in)    :: xb,xe,s
      real(FP_REAL),    intent(inout) :: fp
      integer(FP_SIZE), intent(in)    :: iopt,m,k,nest,lwrk
      integer(FP_FLAG), intent(out)   :: ier
      integer(FP_SIZE), intent(inout) :: n
      !  ..array arguments..
      real(FP_REAL),    intent(in)    :: x(m),y(m),w(m)
      real(FP_REAL),    intent(inout) :: t(nest),c(nest),wrk(lwrk)
      integer(FP_SIZE), intent(inout) :: iwrk(nest)
      !  ..local scalars..
      integer(FP_SIZE) :: i,ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,nmin
      !  ..
      !  we set up the parameters tol and maxit
      real(FP_REAL), parameter :: tol = smallnum03
      integer(FP_SIZE), parameter :: maxit = 20

      k1   = k+1
      k2   = k1+1
      nmin = 2*k1

      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if (k<=0 .or. k>5)         return
      if (iopt<(-1) .or. iopt>1) return
      if (m<k1 .or. nest<nmin)   return
      lwest = m*k1+nest*(7+3*k)
      if (lwrk<lwest)             return
      if (xb>x(1) .or. xe<x(m))  return
      if (any(x(1:m-1)>x(2:m)))  return

      if (iopt>=0) then
          if (s<zero .or. (equal(s,zero) .and. nest<(m+k1))) return
      else
          if (n<nmin .or. n>nest) return
          j = n
          do i=1,k1
             t(i) = xb
             t(j) = xe
             j = j-1
          end do
          ier = fpchec(x,m,t,n,k); if (ier/=0) return
      endif

      ier = FITPACK_OK

      ! we partition the working space and determine the spline approximation.
      ifp = 1
      iz = ifp+nest
      ia = iz+nest
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
      call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp, &
                  wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)

      end subroutine curfit

      pure subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol, &
                             maxit,k1,k2,n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)

      !  ..
      !  ..scalar arguments..
      real(FP_REAL),    intent(in)    :: xb,xe,s,tol
      real(FP_REAL),    intent(out)   :: fp
      integer(FP_SIZE), intent(in)    :: iopt,m,k,nest,maxit,k1,k2
      integer(FP_SIZE), intent(inout) :: n
      integer(FP_FLAG), intent(out)   :: ier

      !  ..array arguments..
      real(FP_REAL), intent(in)    :: x(m),y(m),w(m)
      real(FP_REAL), intent(inout) :: t(nest),c(nest),fpint(nest),z(nest),a(nest,k1),b(nest,k2),&
                                    g(nest,k2),q(m,k1)
      integer(FP_SIZE), intent(inout) :: nrdata(nest)
      !  ..local scalars..
      real(FP_REAL) :: acc,cos,fpart,fpms,fpold,fp0,f1,f2,f3,p,pinv,piv,p1,p2,p3,rn,sin,store,&
                     term,wi,xi,yi
      integer(FP_SIZE) :: i,it,iter,i2,j,k3,l,l0,mk1,nk1,nmax,nmin,nplus,npl1,nrint,n8
      !  ..local arrays..
      real(FP_REAL) :: h(MAX_ORDER+1)
      logical(FP_BOOL) :: new,check1,check3,success

      fpold = zero
      fp0   = zero
      nplus = 0
      fpms  = huge(zero)
      nk1   = n-k1

      ! calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc  = tol*s

      ! determine nmax, the number of knots for spline interpolation.
      nmax = m+k1

      ! *****
      !  part 1: determination of the number of knots and their position
      ! *****
      !  given a set of knots we compute the least-squares spline sinf(x), and the corresponding sum
      !  of squared residuals fp=f(p=inf).
      !  if iopt=-1 sinf(x) is the requested approximation.
      !  if iopt=0 or iopt=1 we check whether we can accept the knots:
      !    if fp <=s we will continue with the current set of knots.
      !    if fp > s we will increase the number of knots and compute the corresponding least-
      !    squares spline until finally fp<=s.
      !  the initial choice of knots depends on the value of s and iopt.
      !    if s=0 we have spline interpolation; in that case the number of knots equals nmax = m+k+1.
      !    if (s>0 and iopt=0) we first compute the least-squares polynomial curve of degree k;
      !      n = nmin = 2*k+2
      !    iopt=1 we start with the set of knots found at the last call of the routine, except for
      !    the case that s > fp0; then we compute directly the least-squares polynomial of degree k.
      ! *****

      !  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      bootstrap: if (iopt>=0) then

          interpolating: if (s<=zero) then

              !  if s=0, s(x) is an interpolating spline.
              !  test whether the required storage space exceeds the available one.
              n = nmax
              if (nmax>nest) then
                 ier = FITPACK_INSUFFICIENT_STORAGE
                 return
              end if

              !  find the position of the interior knots in case of interpolation.
              mk1 = m-k1
              if (mk1/=0) then
                  k3 = k/2
                  i  = k2
                  j  = k3+2
                  do l=1,mk1
                    t(i) = merge( x(j) , (x(j)+x(j-1))*half , k3*2/=k)
                    i = i+1
                    j = j+1
                  end do
              endif

          else interpolating

              !  if s>0 our initial choice of knots depends on the value of iopt.
              !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
              !  polynomial of degree k which is a spline without interior knots.
              !  if iopt=1 and fp0>s we start computing the least squares spline
              !  according to the set of knots found at the last call of the routine.
              use_last_call: if (iopt/=0 .and. n/=nmin) then
                 fp0   = fpint(n)
                 fpold = fpint(n-1)
                 nplus = nrdata(n)
                 if (fp0<=s) then
                     n         = nmin
                     fpold     = zero
                     nplus     = 0
                     nrdata(1) = m-2
                 endif
              else use_last_call
                 n         = nmin
                 fpold     = zero
                 nplus     = 0
                 nrdata(1) = m-2
              endif use_last_call
          endif interpolating

      endif bootstrap

      !  main loop for the different sets of knots. m is a save upper bound
      !  for the number of trials.
      iter = 0
      main_loop: do while (iter<=m)

        iter = iter+1

        if (n==nmin) ier = FITPACK_LEASTSQUARES_OK

        ! find nrint, tne number of knot intervals.
        nrint = n-nmin+1

        ! find the position of the additional knots which are needed for
        ! the b-spline representation of s(x).
        nk1        = n-k1
        t(1:k1)    = xb
        t(nk1+1:n) = xe

        ! compute the b-spline coefficients of the least-squares spline
        ! sinf(x). the observation matrix a is built up row by row and
        ! reduced to upper triangular form by givens transformations.
        ! at the same time fp=f(p=inf) is computed.
        fp = zero

        ! initialize the observation matrix a.
        z(1:nk1)      = zero
        a(1:nk1,1:k1) = zero
        l = k1
        coefs: do it=1,m

            ! fetch the current data point x(it),y(it).
            xi = x(it)
            wi = w(it)
            yi = y(it)*wi

            ! search for knot interval t(l) <= xi < t(l+1).
            do while (xi>=t(l+1) .and. l/=nk1)
                l = l+1
            end do

            ! evaluate the (k+1) non-zero b-splines at xi and store them in q.
            h = fpbspl(t,n,k,xi,l)

            q(it,1:k1) = h(1:k1)
            h(:k1) = wi*h(:k1)

            ! rotate the new row of the observation matrix into triangle.
            j = l-k1
            rotate_row: do i=1,k1

                j = j+1

                piv = h(i); if (equal(piv,zero)) cycle rotate_row

                ! calculate the parameters of the givens transformation.
                call fpgivs(piv,a(j,1),cos,sin)

                ! transformations to right hand side.
                call fprota(cos,sin,yi,z(j))

                ! transformations to left hand side.
                if (i<k1) call fprota(cos,sin,h(i+1:k1),a(j,2:k1-i+1))

            end do rotate_row

            !  add contribution of this row to the sum of squares of residual
            !  right hand sides.
            fp = fp+yi*yi

        end do coefs

        if (ier==FITPACK_LEASTSQUARES_OK) fp0 = fp
        fpint(n-1:n) = [fpold,fp0]
        nrdata(n) = nplus

        ! backward substitution to obtain the b-spline coefficients.
        c(:nk1) = fpback(a,z,nk1,k1,nest)

        ! test whether the approximation sinf(x) is an acceptable solution.
        if (iopt<0) return

        fpms = fp-s; if(abs(fpms)<acc) return

        ! if f(p=inf) < s accept the choice of knots.
        if (fpms<zero) exit main_loop

        ! if n = nmax, sinf(x) is an interpolating spline.
        if (n==nmax) then
            ier = FITPACK_INTERPOLATING_OK
            return
        end if

        ! increase the number of knots.
        ! if n=nest we cannot increase the number of knots because of the storage capacity limitation.
        if (n==nest) then
            ier = FITPACK_INSUFFICIENT_STORAGE
            return
        end if

        ! determine the number of knots nplus we are going to add.
        if (ier==FITPACK_OK) then
            npl1 = nplus*2
            rn = nplus
            if (fpold-fp>acc) npl1 = int(rn*fpms/(fpold-fp))
            nplus = min(nplus*2,max(npl1,nplus/2,1))
        else
            nplus = 1
            ier   = FITPACK_OK
        endif

        ! Initialize update
        fpold = fp

        ! compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
        ! t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = zero
        i = 1
        l = k2
        new = .false.
        square_residuals: do it=1,m
            if (x(it)>=t(l) .and. l<=nk1) then
              new = .true.
              l = l+1
            endif
            l0   = l-k2
            term = dot_product(c(l0+1:l0+k1),q(it,1:k1))

            term = (w(it)*(term-y(it)))**2
            fpart = fpart+term

            if (new) then
                store = term*half
                fpint(i) = fpart-store
                i = i+1
                fpart = store
                new = .false.
            endif
        end do square_residuals

        fpint(nrint) = fpart

        add_new_knots: do l=1,nplus

            ! add a new knot.
            call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,IONE)

            ! if n=nmax we locate the knots as for interpolation.
            if (n==nmax) then
                mk1 = m-k1
                if (mk1/=0) then
                    k3 = k/2
                    i  = k2
                    j  = k3+2
                    do l0=1,mk1
                      t(i) = merge( x(j) , (x(j)+x(j-1))*half , k3*2/=k)
                      i = i+1
                      j = j+1
                    end do
                endif

                ! Restart main loop
                iter = 0
                cycle main_loop

            end if

            ! test whether we cannot further increase the number of knots.
            if (n==nest) exit add_new_knots

        end do add_new_knots
      !  restart the computations with the new set of knots.
      end do main_loop

      !  test whether the least-squares kth degree polynomial is a solution
      !  of our approximation problem.
      if (ier==FITPACK_LEASTSQUARES_OK) return

      ! *****
      !  part 2: determination of the smoothing spline sp(x).
      ! *****
      !  we have determined the number of knots and their position.
      !  we now compute the b-spline coefficients of the smoothing spline sp(x). the observation matrix a
      !  is extended by the rows of matrix b expressing that the kth derivative discontinuities of sp(x)
      !  at the interior knots t(k+2),...t(n-k-1) must be zero. the corresponding weights of these
      !  additional rows are set to 1/p.
      !  iteratively we then have to determine the value of p such that f(p), the sum of squared
      !  residuals be = s. we already know that the least squares kth degree polynomial corresponds
      !  to p=0, and that the least-squares spline corresponds to p=infinity. the iteration process
      !  which is proposed here, makes use of rational interpolation. since f(p) is a convex and strictly
      !  decreasing function of p, it can be approximated by a rational function r(p) = (u*p+v)/(p+w).
      !  three values of p(p1,p2,p3) with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)
      !  are used to calculate the new value of p such that r(p)=s. convergence is guaranteed by taking
      !  f1>0 and f3<zero
      ! *****

      !  evaluate the discontinuity jump of the kth derivative of the
      !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)

      !  initial value for p.
      p1 = zero
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p  = sum(a(1:nk1,1))
      rn = nk1
      p  = rn/p
      check1 = FP_FALSE
      check3 = FP_FALSE
      n8 = n-nmin

      !  iteration process to find the root of f(p) = s.
      iter = 0
      find_root: do while (iter<maxit)

          iter = iter+1

          !  the rows of matrix b with weight 1/p are rotated into the
          !  triangularised observation matrix a which is stored in g.
          pinv = one/p
          c(1:nk1) = z(1:nk1)
          g(1:nk1,1:k1) = a(1:nk1,1:k1)
          g(1:nk1,  k2) = zero
          b_rows: do it=1,n8
              ! the row of matrix b is rotated into triangle by givens transformation
              h(1:k2) = b(it,1:k2)*pinv
              yi = zero
              b_cols: do j=it,nk1
                  piv = h(1)

                  ! calculate the parameters of the givens transformation.
                  call fpgivs(piv,g(j,1),cos,sin)

                  ! transformations to right hand side.
                  call fprota(cos,sin,yi,c(j))
                  if (j==nk1) cycle b_rows

                  ! transformations to left hand side.
                  i2 = merge(nk1-j,k1,j>n8)+1
                  call fprota(cos,sin,h(2:i2),g(j,2:i2))
                  h(1:i2) = [h(2:i2),zero]
              end do b_cols
          end do b_rows

          ! backward substitution to obtain the b-spline coefficients.
          c(:nk1) = fpback(g,c,nk1,k2,nest)

          ! computation of f(p).
          fp = zero
          l = k2

          get_fp: do it=1,m
              if (x(it)>=t(l) .and. l<=nk1) l = l+1
              l0   = l-k2
              term = dot_product(c(l0+1:l0+k1),q(it,:))
              fp   = fp+(w(it)*(term-y(it)))**2
          end do get_fp

          ! SUCCESS! the approximation sp(x) is an acceptable solution.
          fpms = fp-s; if (abs(fpms)<acc) return

          ! find the new value of p and carry out one more step.
          call root_finding_iterate(p1,f1,p2,f2,p3,f3,p,fpms,acc,check1,check3,success)
          if (.not.success) then
             ier = FITPACK_S_TOO_SMALL
             return
          end if

      end do find_root

      ! Maximum number of iterations reached
      ier = FITPACK_MAXIT
      return

      end subroutine fpcurf

      !  subroutine fpdisc calculates the discontinuity jumps of the kth
      !  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
      pure subroutine fpdisc(t,n,k2,b,nest)

      !  ..scalar arguments..
      integer(FP_SIZE), intent(in) :: n,k2,nest
      !  ..array arguments..
      real(FP_REAL), intent(in) :: t(n)
      real(FP_REAL), intent(inout) :: b(nest,k2)
      !  ..local scalars..
      real(FP_REAL) :: an,fac,prod
      integer(FP_SIZE) :: i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
      !  ..local array..
      real(FP_REAL) :: h(12)
      !  ..
      k1    = k2-1
      k     = k1-1
      nk1   = n-k1
      nrint = nk1-k
      an    = nrint
      fac   = an/(t(nk1+1)-t(k1))
      do l=k2,nk1
         lmk = l-k1
         do j=1,k1
            ik = j+k1
            lj = l+j
            lk = lj-k2
            h(j)  = t(l)-t(lk)
            h(ik) = t(l)-t(lj)
         end do
         lp = lmk
         do j=1,k2
           jk = j
           prod = h(j)
           do i=1,k
             jk = jk+1
             prod = prod*h(jk)*fac
           end do
           lk = lp+k1
           b(lmk,j) = (t(lk)-t(lp))/prod
           lp = lp+1
         end do
      end do
      return
      end subroutine fpdisc

      !  subroutine fpgivs calculates the parameters of a givens transformation .
      elemental subroutine fpgivs(piv,ww,cos,sin)
          real(FP_REAL), intent(in)    :: piv
          real(FP_REAL), intent(inout) :: ww
          real(FP_REAL), intent(out)   :: cos,sin
          !  ..local scalars..
          real(FP_REAL) :: dd,store

          store = abs(piv)
          dd  = merge(store*sqrt(one+(ww/piv)**2), &
                      ww   *sqrt(one+(piv/ww)**2), store>=ww)
          cos = ww/dd
          sin = piv/dd
          ww  = dd
          return
      end subroutine fpgivs

      !  subroutine fpknot locates an additional knot for a spline of degree k and adjusts the
      !  corresponding parameters,i.e.
      !    t     : the position of the knots.
      !    n     : the number of knots.
      !    nrint : the number of knot intervals.
      !    fpint : the sum of squares of residual right hand sides
      !            for each knot interval.
      !    nrdata: the number of data points inside each knot interval.
      !  istart indicates that the smallest data point at which the new knot may be added is x(istart+1)
      pure subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)

      !  ..scalar arguments..
      integer(FP_SIZE), intent(in)    :: m,nest,istart
      integer(FP_SIZE), intent(inout) :: n,nrint
      !  ..array arguments..
      real(FP_REAL), intent(in)    :: x(m)
      real(FP_REAL), intent(inout) :: t(nest)
      real(FP_REAL), intent(inout) :: fpint(nest)
      integer(FP_SIZE),  intent(inout) :: nrdata(nest)

      !  ..local scalars..
      real(FP_REAL) :: an,am,fpmax
      integer(FP_SIZE) :: ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt,next,nrx,number
      !  ..
      number = 0
      maxpt  = 0
      maxbeg = 0
      k      = (n-nrint-1)/2
      !  search for knot interval t(number+k) <= x <= t(number+k+1) where fpint(number) is maximal on the
      !  condition that nrdata(number)/=0 .
      fpmax  = zero
      jbegin = istart      
      do j=1,nrint
        jpoint = nrdata(j)

        if (fpmax<fpint(j) .and. jpoint/=0) then
           fpmax = fpint(j)
           number = j
           maxpt = jpoint
           maxbeg = jbegin
        endif

        jbegin = jbegin+jpoint+1
      end do
      
      !  let coincide the new knot t(number+k+1) with a data point x(nrx)
      !  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx   = maxbeg+ihalf
      next  = number+1

      !  adjust the different parameters.
      if (next<=nrint) then
         do j=next,nrint
            jj = next+nrint-j
            fpint(jj+1) = fpint(jj)
            nrdata(jj+1) = nrdata(jj)
            jk = jj+k
            t(jk+1) = t(jk)
         end do
      endif
      
      if (number>0) then 
          nrdata(number) = ihalf-1
          nrdata(next)   = maxpt-ihalf
          am = maxpt
          an = nrdata(number)
          fpint(number) = fpmax*an/am
      endif
      
      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)
      
      n     = n+1
      nrint = nrint+1

      end subroutine fpknot

      ! three values of p (p1,p2,p3) with corresponding values of
      ! f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value of p
      ! such that r(p)=s. convergence is guaranteed by taking f1>0,f3<zero
      elemental subroutine root_finding_iterate(p1,f1,p2,f2,p3,f3,p,fpms,acc,check1,check3,success)
          real(FP_REAL), intent(inout) :: p1,f1,p2,f2,p3,f3,p
          real(FP_REAL), intent(in)    :: fpms,acc
          logical(FP_BOOL), intent(inout) :: check1,check3
          logical(FP_BOOL), intent(out) :: success

          !  set constants
          real(FP_REAL), parameter :: con1 = 0.1e0_FP_REAL
          real(FP_REAL), parameter :: con9 = 0.9e0_FP_REAL
          real(FP_REAL), parameter :: con4 = 0.4e-01_FP_REAL

          success = .true.

          p2 = p
          f2 = fpms
          if (.not.check3) then
             if ((f2-f3)>acc) then
                check3=f2<zero
             else
                ! our initial choice of p is too large.
                p3 = p2
                f3 = f2
                p  = p*con4
                if (p<=p1) p=p1*con9 + p2*con1
                return
             endif
          endif

          if (.not.check1) then
             if ((f1-f2)>acc) then
                 check1 = f2>zero
             else
                 ! our initial choice of p is too small
                 p1 = p2
                 f1 = f2
                 p  = p/con4
                 if (p3>=zero .and. p>=p3) p = p2*con1 + p3*con9
                 return
             endif
          endif

          ! test whether the iteration process proceeds as theoretically expected.
          if (f2>=f1 .or. f2<=f3) then
             success = .false.
             return
          else
             ! find the new value of p.
             call fprati(p1,f1,p2,f2,p3,f3,p)
          endif

      end subroutine root_finding_iterate

      ! subroutine fprota applies a givens rotation to a and b.
      elemental subroutine fprota(cos,sin,a,b)

          !  ..scalar arguments..
          real(FP_REAL), intent(in)    :: cos,sin
          real(FP_REAL), intent(inout) :: a,b

          ! ..local scalars..
          real(FP_REAL) :: stor1,stor2

          !  ..
          stor1 = a
          stor2 = b
          b = cos*stor2+sin*stor1
          a = cos*stor1-sin*stor2
          return

      end subroutine fprota

      ! subroutine splev evaluates in a number of points x(i),i=1,2,...,m a spline s(x) of degree k,
      ! given in its b-spline representation.
      pure subroutine splev(t,n,c,k,x,y,m,e,ier)

      !  calling sequence:
      !     call splev(t,n,c,k,x,y,m,e,ier)
      !
      !  input parameters:
      !    t    : array,length n, which contains the position of the knots.
      !    n    : integer(FP_SIZE), giving the total number of knots of s(x).
      !    c    : array,length n, which contains the b-spline coefficients.
      !    k    : integer(FP_SIZE), giving the degree of s(x).
      !    x    : array,length m, which contains the points where s(x) must be evaluated.
      !    m    : integer(FP_SIZE), giving the number of points where s(x) must be evaluated.
      !    e    : integer(FP_SIZE), boundary condition for points outside the support
      !           0 = the spline is extrapolated from the end spans
      !           1 = the spline evaluates to zero for those points,
      !           2 = extrapolation not allowed, ier is set to 1 and the subroutine returns,
      !           3 = the spline evaluates to the value of the nearest boundary point.
      !
      !  output parameter:
      !    y    : array,length m, giving the value of s(x) at the different points.
      !    ier  : error flag
      !
      !  restrictions:
      !    m >= 1
      !
      !  other subroutines required: fpbspl.
      !
      !  references :
      !    de boor c  : on calculating with b-splines, j. approximation theory 6 (1972) 50-62.
      !    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths applics 10 (1972) 134-149.
      !    dierckx p. : curve and surface fitting with splines, monographs on numerical analysis, oxford
      !                 university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  ..scalar arguments..
      integer(FP_SIZE), intent(in)  :: n, k, m
      integer(FP_FLAG), intent(in)  :: e
      integer(FP_FLAG), intent(out) :: ier
      !  ..array arguments..
      real(FP_REAL), intent(in)  :: t(n), c(n), x(m)
      real(FP_REAL), intent(out) :: y(m)

      !  ..local scalars..
      integer(FP_SIZE) :: i, k1, l, l1, nk1,k2
      real(FP_REAL) :: arg, tb, te
      !  ..local array..
      real(FP_REAL) :: h(MAX_ORDER+1)
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if (m<1) return

      ier = FITPACK_OK
      !  fetch tb and te, the boundaries of the approximation interval.
      k1  = k  + 1
      k2  = k1 + 1
      nk1 = n - k1
      tb  = t(k1)
      te  = t(nk1 + 1)
      l   = k1
      l1  = l + 1

      !  main loop for the different points.
      user_points: do i = 1, m

        ! fetch a new x-value arg.
        arg = x(i)

        ! check if arg is in the support
        outside: if (arg<tb .or. arg>te) then
            select case (e)
               case (OUTSIDE_EXTRAPOLATE)
                ! Continue normally
               case (OUTSIDE_ZERO)
                  y(i) = zero
                  cycle user_points
               case (OUTSIDE_NOT_ALLOWED)
                  ier = FITPACK_INVALID_RANGE
                  return
               case (OUTSIDE_NEAREST_BND)
                  arg = max(min(arg,te),tb)
            end select
        endif outside

        ! search for knot interval t(l) <= arg < t(l+1)
        do while (arg<t(l) .and. l1/=k2)
          l1 = l
          l  = l - 1
        end do
        do while (arg>=t(l1) .and. l/=nk1)
          l  = l1
          l1 = l + 1
        end do

        ! evaluate the non-zero b-splines at arg.
        h = fpbspl(t, n, k, arg, l)

        ! find the value of s(x) at x=arg.
        y(i) = dot_product(c(l-k:l),h(1:k1))
      end do user_points

      end subroutine splev

      !  function fpback calculates the solution of the system of equations a*c = z with
      !  a a n x n upper triangular matrix of bandwidth k.
      pure function fpback(a,z,n,k,nest) result(c)

      !  ..scalar arguments..
      integer(FP_SIZE), intent(in) :: n,k,nest
      !  ..array arguments..
      real(FP_REAL), intent(in)  :: a(nest,k),z(n)
      real(FP_REAL)              :: c(n)
      !  ..local scalars..
      real(FP_REAL) :: store
      integer(FP_SIZE) :: i,i1,j,k1,l,m,jm1
      !  ..
      k1   = k-1
      c(n) = z(n)/a(n,1)
      i    = n-1
      if (i==0) return

      jm1 = 1
      rows: do j=2,n
        store = z(i)
        i1 = merge(jm1,k1,j<=k1)
        m = i
        do l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
        end do
        c(i) = store/a(i,1)
        i = i-1
        jm1 = j
      end do rows

      end function fpback

      !  subroutine fpchec verifies the number and the position of the knots t(j),j=1,2,...,n of a spline
      !  of degree k, in relation to the number and the position of the data points x(i),i=1,2,...,m.
      !  If all of the following conditions are fulfilled, the error parameter ier is set to zero. if one
      !  of the conditions is violated, an error flag is returned.
      pure integer(FP_FLAG) function fpchec(x,m,t,n,k) result(ier)
         integer(FP_SIZE), intent(in)  :: m,n,k
         real(FP_REAL), intent(in) :: x(m),t(n)

         ! Local variables
         integer(FP_SIZE) :: i,j,k1,k2,l,nk1,nk2,nk3
         real(FP_REAL) :: tj,tl

         ! Init sizes
         k1 = k+1
         k2 = k1+1
         nk1 = n-k1
         nk2 = nk1+1

         ier = FITPACK_INPUT_ERROR

         ! 1) k+1 <= n-k-1 <= m
         if(nk1<k1 .or. nk1>m) return

         ! 2) monotonicity
         !    t(1) <= t(2) <= ... <= t(k+1)
         !    t(n-k) <= t(n-k+1) <= ... <= t(n)

         j = n
         monotonic: do i=1,k
            if(t(i)>t(i+1)) return
            if(t(j)<t(j-1)) return
            j = j-1
         end do monotonic

         ! 3) t(k+1) < t(k+2) < ... < t(n-k)
         do i=k2,nk2
            if(t(i)<=t(i-1)) return
         end do

         ! 4) t(k+1) <= x(i) <= t(n-k)
         ! 5) schoenberg and whitney conditions: they must hold for at least one subset of data points, i.e.
         !    there must be a subset of data points y(j) such that
         !         t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1

         if(x(1)<t(k1) .or. x(m)>t(nk2)) return
         if(x(1)>=t(k2) .or. x(m)<=t(nk1)) return

         i   = 1
         l   = k2
         nk3 = nk1-1

         if (nk3>=2) then
             do j=2,nk3
                tj = t(j)
                l  = l+1
                tl = t(l)
                do while (i<m .and. x(i)<=tj)
                    i = i+1
                    if(i>=m) return
                end do
                if (x(i)>=tl) return
             end do
         endif

         ! All checks passed
         ier = FITPACK_OK

      end function fpchec

      !  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati  gives the value of p such
      !  that the rational interpolating function of the form r(p) = (u*p+v)/(p+w) equals zero at p.
      elemental subroutine fprati(p1,f1,p2,f2,p3,f3,p)
      !  ..scalar arguments..
      real(FP_REAL), intent(inout) :: p1,f1,p3,f3
      real(FP_REAL), intent(in)    :: p2,f2
      real(FP_REAL), intent(out)   :: p

      !  ..local scalars..
      real(FP_REAL) :: h1,h2,h3
      !  ..
      if (p3>zero) then
         h1 = f1*(f2-f3)
         h2 = f2*(f3-f1)
         h3 = f3*(f1-f2)
         p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
      else
         !  value of p in case p3 = infinity.
         p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      end if

      !  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
      if (f2>=zero) then
        p1 = p2
        f1 = f2
      else
        p3 = p2
        f3 = f2
      endif

      return
      end subroutine fprati

      !  function fpbspl evaluates the (k+1) non-zero b-splines of degree k at t(l) <= x < t(l+1) using
      !  the stable recurrence relation of de boor and cox.
      !  Travis Oliphant 2007 changed so that weighting of 0 is used when knots with multiplicity are present.
      !   Also, notice that l+k <= n and 1 <= l+1-k or else the routine will be accessing memory outside t
      !   Thus it is imperative that that k <= l <= n-k but this is not checked.
      pure function fpbspl(t,n,k,x,l) result(h)
         integer(FP_SIZE), intent(in)  :: n,k,l
         real(FP_REAL), intent(in)  :: x,t(n)
         real(FP_REAL) :: h(MAX_ORDER+1)

         ! Local variables
         real(FP_REAL) :: f,hh(MAX_ORDER+1)
         integer(FP_SIZE) :: i,j,li,lj

         h(1) = one
         do j=1,k
           hh(1:j) = h(1:j)
           h(1) = zero
           do i=1,j
             li = l+i
             lj = li-j
             if (not_equal(t(li),t(lj))) then
                f = hh(i)/(t(li)-t(lj))
                h(i)   = h(i)+f*(t(li)-x)
                h(i+1) = f*(x-t(lj))
             else
                h(i+1) = zero
             endif
           end do
         end do

      end function fpbspl

      ! Test if two reals have the same representation at the current precision
      elemental logical(FP_BOOL) function equal(a,b)
         real(FP_REAL), intent(in) :: a,b
         equal = abs(a-b)<spacing(merge(a,b,abs(a)<abs(b)))
      end function equal
      elemental logical(FP_BOOL) function not_equal(a,b)
         real(FP_REAL), intent(in) :: a,b
         not_equal = .not.equal(a,b)
      end function not_equal
end module fitpack_core