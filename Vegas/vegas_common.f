      integer mxdim
      parameter(mxdim=25)
      integer ndim,ncall,itmx,nprn
      double precision xl(mxdim),xu(mxdim),acc
      double precision xi(50,mxdim),si,si2,swgt,schi
      integer ndo,it,idum
      logical readin,writeout
      logical warmup,stopvegas
      character*(72) ingridfile,outgridfile
      common/gridinfo_logic/readin,writeout,warmup,stopvegas
      common/gridinfo_char/ingridfile,outgridfile
      common/bveg1/xl,xu,acc,ndim,ncall,itmx,nprn
      common/bveg2/xi,si,si2,swgt,schi,ndo,it
      common/ranno/idum
