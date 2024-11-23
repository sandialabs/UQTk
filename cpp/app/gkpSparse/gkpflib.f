c$$$=====================================================================================
c$$$
c$$$                      The UQ Toolkit (UQTk) version 3.1.5
c$$$                          Copyright (2024) NTESS
c$$$                        https://www.sandia.gov/UQToolkit/
c$$$                        https://github.com/sandialabs/UQTk
c$$$
c$$$     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
c$$$     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
c$$$     retains certain rights in this software.
c$$$
c$$$     This file is part of The UQ Toolkit (UQTk)
c$$$
c$$$     UQTk is open source software: you can redistribute it and/or modify
c$$$     it under the terms of BSD 3-Clause License
c$$$
c$$$     UQTk is distributed in the hope that it will be useful,
c$$$     but WITHOUT ANY WARRANTY; without even the implied warranty of
c$$$     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c$$$     BSD 3 Clause License for more details.
c$$$
c$$$     You should have received a copy of the BSD 3 Clause License
c$$$     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
c$$$
c$$$     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
c$$$     Sandia National Laboratories, Livermore, CA, USA
c$$$=====================================================================================
      subroutine heap_ext(n,isgn,i,j,index)

       implicit none
       integer n,isgn,i,j,index

       integer l,l1,n1

       common /hpwrk/ l,l1,n1

       if (index) 90,10,80
 10    n1 = n
       l  = 1+n/2
 20    l  = l-1
 30    l1 = l
 40    i  = l1+l1
       if (i-n1) 50,60,70
 50    j  = i+1
       index = -2
       return
 60    j  = l1
       l1 = i
       index = -1
       return
 70    if (l.gt.1) goto 20
       if (n1.eq.1) goto 110
       i  = n1
       n1 = n1-1
       j  = 1
       index = 1
       return
 80    if (index-1) 30,30,40
 90    if (index.eq.-1) goto 100
       if (isgn.lt.0) i=i+1
       goto 60
 100   if (isgn.le.0) goto 70
       index = 2
       return
 110   index = 0
       return
       end
